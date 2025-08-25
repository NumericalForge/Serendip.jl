# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export DynAnalysis

# DynAnalysis_params = [
#     FunInfo(:DynAnalysis, "Dynamic mechanical analyses."),
#     KwArgInfo(:stress_state, "Stress model", :d3, values=(:plane_stress, :plane_strain, :axisymmetric, :d3)),
#     KwArgInfo(:thickness, "Thickness for 2d analyses", 1.0, cond=:(thickness>0)),
#     KwArgInfo(:g, "Gravity acceleration", 0.0, cond=:(g>=0))
# ]
# @doc docstring(DynAnalysis_params) DynAnalysis()

# mutable struct DynAnalysisProps<:TransientAnalysis
#     stress_state::Symbol # plane stress, plane strain, etc.
#     thickness::Float64  # thickness for 2d analyses
#     g::Float64 # gravity acceleration

#     function DynAnalysisProps(; kwargs...)
#         args = checkargs(kwargs, DynAnalysis_params)
#         this = new(args.stress_state, args.thickness, args.g)
#         return this
#     end
# end

# DynAnalysis = DynAnalysisProps


# DynamicAnalysis_params = [
#     FunInfo(:DynAnalysis, "Dynamic analysis"),
#     ArgInfo(:model, "Finite element model"),
# ]
# @doc docstring(DynamicAnalysis_params) DynamicAnalysis


mutable struct DynamicAnalysis<:Analysis
    model ::FEModel
    ctx   ::Context
    # sctx  ::SolverContext

    stages  ::Vector{Stage}
    loggers ::Vector{Logger}
    monitors::Vector{Monitor}

    function DynamicAnalysis(model::FEModel; outdir=".", outkey="out")
        this = new(model, model.ctx)
        this.stages = []
        this.loggers = []
        this.monitors = []
        this.sctx = SolverContext()

        this.sctx.outkey = outkey
        this.sctx.outdir = rstrip(outdir, ['/', '\\'])
        isdir(this.sctx.outdir) || mkdir(this.sctx.outdir) # create output directory if it does not exist

        model.ctx.thickness = model.thickness
        if model.ctx.stress_state==:none
            if model.ctx.ndim==2
                model.ctx.stress_state = :plane_strain
            else
                model.ctx.stress_state = :d3
            end
        end

        return this
    end
end

const DynAnalysis = DynamicAnalysis


# Assemble the global mass matrix
function mount_M(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            Me, rmap, cmap = elem_mass(elem)

            nr, nc = size(Me)
            for i in 1:nr
                for j in 1:nc
                    val = Me[i,j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    local M
    try
        M = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show err
    end

    yield()

    return M
end


function sismic_force(model::FEModel, bcs, M::SparseMatrixCSC{Float64, Int}, F::Vect, AS::Array{Float64,2}, keysis::Symbol, tfia::Float64, tds::Float64)

    dofs, nu = configure_dofs(model, bcs)
    ndofs = length(dofs)

    ndat = length(AS)    # quantity of aceleration data
    AS2  = zeros(ndat+1) # sismic aceleration data how vector
    c = 0
    for i in 1:size(AS,1)
        for j in 1:size(AS,2)
            c += 1
            AS2[c+1] = AS[i,j]
        end
    end

    vts = zeros(ndat+1) # time vetor correspond to acelerations

    for i in 1:ndat+1
        vts[i] = (i-1)*tds/ndat
    end

    FAS = hcat(vts,AS2) # Function of aceleration

    # Interpolation of the aceleration value

    inic = 0
    fin = 0

    for i in 1:ndat
        if FAS[i,1]<=tfia<=FAS[i+1,1]
            inic = i
            fin = i+1
        end
        if inic!=0 & fin!=0; break end
    end

    m = (FAS[fin,2]-FAS[inic,2])/(FAS[fin,1]-FAS[inic,1])
    acel = FAS[inic,2] + m*(tfia - FAS[inic,1])

    #Dof sismic aceleration

    VAS  = zeros(ndofs) #Sismic aceleration vector accord dof

    for node in model.nodes
        dof = node.dofdict[keysis]
        VAS[dof.eq_id] += acel
    end

    #Dof sismic force
    FS = M*VAS
    F += FS

    return F

end


# # function solve!(model::FEModel, ana::DynAnalysis; args...)
# #     name = "Solver for dynamic analyses"
# #     status = stage_iterator(name, dyn_stage_solver!, model; args...)
# #     return status
# # end

# dyn_solver_params = [
#     FunInfo( :solve!, "Solves a dynamic mechanical analysis."),
#     KwArgInfo( :tol, "Force tolerance", 0.01, cond=:(tol>0)),
#     # ArgInfo( :rtol, "Relative tolerance in terms of displacements", 0.01, cond=:(rtol>0)),
#     KwArgInfo( :dTmin, "Relative minimum increment size", 1e-7, cond=:(0<dTmin<1) ),
#     KwArgInfo( :dTmax, "Relative maximum increment size", 0.1, cond=:(0<dTmax<1) ),
#     KwArgInfo( :rspan, "Relative span to residue reapplication", 0.01, cond=:(0<rspan<1) ),
#     KwArgInfo( :scheme, "Global solving scheme", :FE, values=(:FE, :ME, :BE, :Ralston) ),
#     KwArgInfo( :maxits, "Maximum number of NR iterations", 5, cond=:(1<=maxits<=10)),
#     KwArgInfo( :autoinc, "Flag to set auto-increments", false),
#     KwArgInfo( :quiet, "Flat to set silent mode", false),
#     KwArgInfo( :alpha, "", 0.0),
#     KwArgInfo( :beta, "", 0.0),
#     KwArgInfo( :sism, "", false),
#     KwArgInfo( :tss, "", 0.0),
#     KwArgInfo( :tds, "", 0.0),
#     KwArgInfo( :sism_file, "", ""),
#     KwArgInfo( :sism_dir, "", "fx"),
# ]
# @doc docstring(dyn_solver_params) solve!(::DynAnalysis; args...)


function solve!(ana::DynamicAnalysis; args...)
    args = checkargs(args, dyn_solver_params)
    if !args.quiet
        printstyled("Solver for mechanical dynamic analyses", "\n", bold=true, color=:cyan)
        println("  stress model: ", ana.ctx.stress_state)
    end

    status = stage_iterator(dyn_stage_solver!, ana; args...)
    return status
end


function dyn_stage_solver!(ana::DynamicAnalysis, stage::Stage; args...)
    args = NamedTuple(args)


    tol     = args.tol
    ΔTmin   = args.dTmin
    ΔTmax   = args.dTmax
    rspan   = args.rspan
    scheme  = args.scheme
    maxits  = args.maxits
    autoinc = args.autoinc
    quiet   = args.quiet
    sism    = args.sism
    alpha   = args.alpha
    beta    = args.beta


    sctx = ana.sctx
    model = ana.model
    ctx  = model.ctx

    println(sctx.log, "Dynamic FE analysis: Stage $(stage.id)")

    solstatus = success()
    # scheme in ("FE", "ME", "BE", "Ralston") || error("solve! : invalid scheme \"$(scheme)\"")

    nincs    = stage.nincs
    nouts    = stage.nouts
    bcs      = stage.bcs
    tspan    = stage.tspan
    saveouts = stage.nouts > 0

    stress_state = ctx.stress_state
    ctx.ndim==3 && @check stress_state==:d3

    # Get active elements
    for elem in stage.toactivate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs(model, stage.bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    model.ndofs = length(dofs)
    println(sctx.log, "unknown dofs: $nu")
    println(sctx.info, "unknown dofs: $nu")

    # Dictionary of data keys related with a dof
    components_dict = Dict(:ux => (:ux, :fx, :vx, :ax),
                           :uy => (:uy, :fy, :vy, :ay),
                           :uz => (:uz, :fz, :vz, :az),
                           :rx => (:rx, :mx, :vrx, :arx),
                           :ry => (:ry, :my, :vry, :ary),
                           :rz => (:rz, :mz, :vrz, :arz))

    # Setup quantities at dofs
    if stage.id == 1
        for dof in dofs
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[us] = 0.0
            dof.vals[fs] = 0.0
            dof.vals[vs] = 0.0
            dof.vals[as] = 0.0
        end

        #If the problem is sismic, read the sismic acelerations asking to user the file's name AS:SeismicAcelerations
        if sism
            AS = readdlm(sism_file)
            AS= 9.81*AS
            keysis = Symbol(sism_dir)
        end

        # Initial accelerations
        K = mount_K(model.elems, ndofs)
        M = mount_M(model.elems, ndofs)
        A = zeros(ndofs)
        V = zeros(ndofs)
        Uex, Fex = get_bc_vals(model, bcs) # get values at time t

        # if the problem has a sism, add the sismic force
        if sism && tss<=0 # tss:time when seismic activity starts tds: time of seismic duration 0:current time =0s
            Fex = sismic_force(model, bcs, M, Fex, AS, keysis, 0.0, tds)
        end
        solve_system!(M, A, Fex, nu)

        # Initial values at nodes
        for (i,dof) in enumerate(dofs)
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[vs] = V[i]
            dof.vals[as] = A[i]
            dof.vals[fs] = Fex[i]
        end

        # Save initial file and loggers
        update_records!(ana, force=true)
    end

    # Get the domain current state and backup
    State = [ ip.state for elem in active_elems for ip in elem.ips ]
    StateBk = copy.(State)

    # Incremental analysis
    ΔTbk    = 0.0
    ΔTcheck = saveouts ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT, ΔTmax, ΔTcheck))

    t  = sctx.t

    inc  = 0             # increment counter
    iout = sctx.out       # file output counter
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    Fin  = zeros(ndofs)  # total internal force
    ΔFin = zeros(ndofs)  # internal forces for current increment
    ΔUa  = zeros(ndofs)  # essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # essential values for current iteration obtained from NR algorithm
    Fina = zeros(ndofs)  # current internal forces
    TFin = zeros(ndofs)
    Aa   = zeros(ndofs)
    Va   = zeros(ndofs)
    Rc   = zeros(ndofs)  # vector of cumulated residues

    sysstatus = ReturnStatus()

    # while t < tend - ttol
    while T < 1.0-ΔTmin
        sctx.ΔT = ΔT
        Δt = tspan*ΔT
        sctx.t = t

        # Update counters
        inc += 1
        sctx.inc = inc

        println(sctx.log, "  inc $inc")

        Uex, Fex = get_bc_vals(model, bcs, t+Δt) # get values at time t+Δt

        # If the problem has a sism, the force sismic is added
        if sism && tss<=t+Δt && tds+tss>=t+Δt
            M = mount_M(model.elems, ndofs)
            Fex = sismic_force(model, bcs, M, Fex, AS, keysis, t+Δt, tds)
        end

        #R   .= FexN - F    # residual
        Fex_Fin = Fex-Fina    # residual
        ΔUa .= 0.0
        ΔUi .= Uex    # essential values at iteration i # TODO: check for prescribed disps

        # Newton Rapshon iterations
        residue   = 0.0
        maxfails  = 3    # maximum number of it. fails with residual change less than 90%
        nfails    = 0    # counter for iteration fails
        nits      = 0
        residue1  = 0.0
        converged = false

        for it=1:maxits
            yield()
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = residue # residue from last iteration

            # Try FE step
            K = mount_K(model.elems, ndofs)
            M = mount_M(model.elems, ndofs)

            C   = alpha*M + beta*K # Damping matrix
            Kp  = K + (4/(Δt^2))*M + (2/Δt)*C # pseudo-stiffness matrix
            ΔFp = Fex_Fin + M*(A + 4*V/Δt - 4*ΔUa/Δt^2) + C*(V - 2*ΔUa/Δt)

            # Solve
            solve_system!(Kp, ΔUi, ΔFp, nu)

            # Restore the state to last converged increment
            copyto!.(State, StateBk)

            # Get internal forces and update data at integration points (update ΔFin)
            ΔFin .= 0.0
            ΔUt   = ΔUa + ΔUi
            ΔFin, sysstatus = update_state(model.elems, ΔUt, Δt)

            Fina = Fin + ΔFin

            # Update V and A
            # Optional (V and A may be used from the interval beginning)
            Va = -V + 2*(ΔUt)/Δt
            Aa = -A + 4*((ΔUt) - V*Δt)/(Δt^2)


            TFin = Fina + C*Va + M*Aa  # Internal force including dynamic effects

            residue = maximum(abs, (Fex-TFin)[umap] )

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            Fex_Fin .= Fex .- Fina  # Check this variable, it is not the residue actually
            Fex_Fin[pmap] .= 0.0  # Zero at prescribed positions

            @printf(sctx.log, "    it %d  residue: %-10.4e\n", it, residue)

            it==1 && (residue1=residue)
            residue > tol && (Fina -= ΔFin)
            residue < tol  && (converged=true; break)
            isnan(residue) && break
            it>maxits      && break
            it>1 && residue>lastres && break
            residue>0.9*lastres && (nfails+=1)
            nfails==maxfails    && break
        end

        if converged
            # Update forces and displacement for the current stage
            Fin = Fina
            U .+= ΔUa

            # Backup converged state at ips
            copyto!.(StateBk, State)

            # Update vectors for velocity and acceleration
            A = Aa
            V = Va

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                us, fs, vs, as = components_dict[dof.name]
                dof.vals[us] = U[i]
                dof.vals[fs] = TFin[i]
                dof.vals[vs] = V[i]
                dof.vals[as] = A[i]
            end

            # update_single_loggers!(model)

            # Update time t and Δt
            T += ΔT
            sctx.T = T
            t += Δt
            sctx.t = t
            sctx.residue = residue

            # Check for saving output file
            checkpoint = T>Tcheck-ΔTmin
            if checkpoint
                sctx.out += 1
                # update_embedded_disps!(active_elems, model.node_data["U"])
                Tcheck += ΔTcheck # find the next output time
            end

            # Check for saving output file
            # if T>Tcheck-ΔTmin && saveouts
            #     sctx.out += 1
            #     iout = sctx.out

            #     rm.(glob("*conflicted*.dat", "$outdir/"), force=true)

            #     update_output_data!(model)
            #     update_embedded_disps!(active_elems, model.node_data["U"])

            #     update_multiloggers!(model)
            #     save(model, "$outdir/$outkey-$iout.vtu", quiet=true) #!

            #     Tcheck += ΔTcheck # find the next output time
            # end

            # update_single_loggers!(model)
            # update_monitors!(model)
            # flush(sctx.log)

            rstatus = update_records!(ana, checkpoint=checkpoint)
            if failed(rstatus)
                println(sctx.alerts, rstatus.message)
                return rstatus
            end

            if autoinc
                if ΔTbk>0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    if nits==1
                        q = (1+tanh(log10(tol/residue1)))^1
                    else
                        q = 1.0
                    end

                    ΔTtr = min(q*ΔT, 1/nincs, 1-T)
                    if T+ΔTtr>Tcheck-ΔTmin
                        ΔTbk = ΔT
                        ΔT = Tcheck-T
                        @assert ΔT>=0.0
                    else
                        ΔT = ΔTtr
                        ΔTbk = 0.0
                    end
                end
            end
        else
            # Restore counters
            inc -= 1
            sctx.inc -= 1

            copyto!.(State, StateBk)

            if autoinc
                println(sctx.log, "      increment failed")
                q = (1+tanh(log10(tol/residue1)))
                q = clamp(q, 0.2, 0.9)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < ΔTmin
                    solstatus = failure("solver did not converge")
                    break
                end
            else
                solstatus = failure("solver did not converge")
                break
            end
        end
    end

    failed(solstatus) && update_records!(ana, force=true)

    return solstatus
end