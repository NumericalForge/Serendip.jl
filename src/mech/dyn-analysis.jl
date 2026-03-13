# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export DynamicAnalysis


"""
    DynamicAnalysis(model::FEModel; outdir="./output", outkey="out")

Create a dynamic mechanical analysis for the given finite element model.

# Arguments
- `model::FEModel`: Finite element model definition.
- `outdir::String`: Directory for output files. Defaults to `"./output"`.
- `outkey::String`: Key prefix for output files. Defaults to `"out"`.

# Behavior
- Initializes an `AnalysisData` instance with the given output settings.
- If `model.ctx.stress_state == :auto` and `model.ctx.ndim == 2`, the stress state is set to `:plane_strain`.
- Marks the analysis as transient (`this.data.transient = true`).

# Example
```julia
model = FEModel(mapper)
analysis = DynamicAnalysis(model; outdir="results", outkey="run1")
```
"""
mutable struct DynamicAnalysis<:Analysis
    name::String
    model ::FEModel
    data::AnalysisData

    function DynamicAnalysis(model::FEModel; outdir="./output", outkey="out")
        name = "Dynamic mechanical analysis"
        data = AnalysisData(outdir=outdir, outkey=outkey)
        this = new(name, model, data)

        if model.ctx.stress_state==:auto && model.ctx.ndim==2
            model.ctx.stress_state = :plane_strain
        end
        this.data.transient = true

        return this
    end
end

get_natural_keys(::DynamicAnalysis) = [:fx, :fy, :fz]
get_essential_keys(::DynamicAnalysis) = [:ux, :uy, :uz]

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



function stage_solver(ana::DynamicAnalysis, stage::Stage, solver_settings::SolverSettings; quiet=quiet)

    tol     = solver_settings.tol
    ΔT0     = solver_settings.dT0
    ΔTmin   = solver_settings.dTmin
    ΔTmax   = solver_settings.dTmax
    maxits  = solver_settings.maxits
    autoinc = solver_settings.autoinc
    alpha   = solver_settings.alpha
    beta    = solver_settings.beta

    model = ana.model
    data  = ana.data
    println(data.log, "Dynamic mechanical FE analysis: Stage $(stage.id)")

    solstatus = success()

    nincs    = stage.nincs
    nouts    = stage.nouts
    bcs      = stage.bcs
    ctx      = model.ctx
    tspan    = stage.tspan
    saveouts = stage.nouts > 0
    ftol     = tol

    @check tspan > 0.0 "stage_solver: stage.tspan must be positive for dynamic analysis"

    stress_state = ctx.stress_state
    ctx.ndim == 3 && @assert(stress_state == :auto)

    # Get active elements
    for elem in stage.activate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs(model, bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    data.nu = nu

    println(data.log, "unknown dofs: $nu")
    quiet || nu == ndofs && println(data.alerts, "No essential boundary conditions")

    # Dictionary of data keys related with a dof
    components_dict = Dict(:ux => (:ux, :fx, :vx, :ax),
                           :uy => (:uy, :fy, :vy, :ay),
                           :uz => (:uz, :fz, :vz, :az),
                           :rx => (:rx, :mx, :vrx, :arx),
                           :ry => (:ry, :my, :vry, :ary),
                           :rz => (:rz, :mz, :vrz, :arz))

    if stage.id == 1
        commit_state(active_elems) # commits changes to states from elem_init

        # Setup quantities at dofs
        for dof in dofs
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[us] = 0.0
            dof.vals[fs] = 0.0
            dof.vals[vs] = 0.0
            dof.vals[as] = 0.0
        end

        # Initial accelerations
        A0 = zeros(ndofs)
        Uex = zeros(ndofs)
        Fex = zeros(ndofs)
        for bc in bcs
            compute_bc_values(ana, bc, data.t, Uex, Fex)
        end

        M0 = mount_M(active_elems, ndofs)
        sysstatus = solve_system!(M0, A0, Fex, nu)
        failed(sysstatus) && return sysstatus

        # Initial values at dofs
        for (i,dof) in enumerate(dofs)
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[us] = 0.0
            dof.vals[fs] = Fex[i]
            dof.vals[vs] = 0.0
            dof.vals[as] = A0[i]
        end

        update_records!(ana, force=true)
    end

    # Incremental analysis
    ΔTbk    = 0.0
    ΔTcheck = saveouts ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT, ΔTmax, ΔTcheck, ΔT0))

    t = data.t

    inc     = 0             # increment counter
    U       = zeros(ndofs)  # total displacements for current stage
    Fin     = zeros(ndofs)  # total internal force
    A       = zeros(ndofs)  # current acceleration
    V       = zeros(ndofs)  # current velocity
    ΔFin    = zeros(ndofs)  # internal forces for current increment
    ΔUa     = zeros(ndofs)  # essential values for this increment
    ΔUi     = zeros(ndofs)  # essential values for current iteration
    Fina    = zeros(ndofs)  # trial internal force
    TFin    = zeros(ndofs)  # internal force including dynamic effects
    Aa      = zeros(ndofs)  # trial acceleration
    Va      = zeros(ndofs)  # trial velocity
    R       = zeros(ndofs)  # residual force vector
    Fex_Fin = zeros(ndofs)

    for (i,dof) in enumerate(dofs)
        us, fs, vs, as = components_dict[dof.name]
        U[i] = dof.vals[us]
        V[i] = dof.vals[vs]
        A[i] = dof.vals[as]
    end

    # Rebuild internal force from committed element states before applying increments.
    Fin .= 0.0
    for elem in active_elems
        Fe, map, status = elem_internal_forces(elem)
        failed(status) && return status
        Fin[map] .+= Fe
    end
    Fina .= Fin

    sysstatus = ReturnStatus()

    while T < 1.0-ΔTmin
        data.ΔT = ΔT
        Δt = tspan*ΔT
        data.t = t

        # Update counters
        inc += 1
        data.inc = inc

        println(data.log, "  inc $inc")

        # Get boundary conditions
        Uex = zeros(ndofs)  # essential displacements
        Fex = zeros(ndofs)  # external forces
        for bc in bcs
            compute_bc_values(ana, bc, t+Δt, Uex, Fex) # get values at time t+Δt
        end

        Fex_Fin .= Fex .- Fina
        ΔUa .= 0.0
        ΔUi .= 0.0
        # Prescribed values in the linear system are displacement increments.
        ΔUi[pmap] .= Uex[pmap] .- U[pmap]

        # Newton-Raphson iterations
        residue   = 0.0
        residue1  = 0.0
        nits      = 0
        maxfails  = 3    # maximum number of it. fails with residual change less than 90%
        nfails    = 0    # counter for iteration fails
        converged = false
        syserror  = false

        for it in 1:maxits
            yield()

            nits += 1
            it > 1 && (ΔUi .= 0.0) # essential values are applied only at first iteration
            lastres = residue

            K = mount_K(active_elems, ndofs)
            M = mount_M(active_elems, ndofs)

            C   = alpha*M + beta*K # damping matrix
            Kp  = K + (4/(Δt^2))*M + (2/Δt)*C # pseudo-stiffness matrix
            ΔFp = Fex_Fin + M*(A + 4*V/Δt - 4*ΔUa/Δt^2) + C*(V - 2*ΔUa/Δt)

            sysstatus = solve_system!(Kp, ΔUi, ΔFp, nu)
            failed(sysstatus) && (syserror=true; break)

            ΔUt = ΔUa + ΔUi
            ΔFin, sysstatus = update_state(active_elems, ΔUt, Δt)
            failed(sysstatus) && (syserror=true; break)

            Fina .= Fin .+ ΔFin

            Va .= -V + 2*ΔUt/Δt
            Aa .= -A + 4*(ΔUt - V*Δt)/(Δt^2)

            TFin .= Fina + C*Va + M*Aa

            R .= Fex .- TFin
            residue = norm(R[umap], Inf)

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            Fex_Fin .= Fex .- Fina
            Fex_Fin[pmap] .= 0.0

            @printf(data.log, "    it %d  residue: %-10.4e\n", it, residue)

            it == 1 && (residue1 = residue)
            residue < ftol  && (converged=true; break)
            isnan(residue)  && break
            it > 1 && residue > lastres && break
            residue > 0.9*lastres && (nfails += 1)
            nfails == maxfails && break
        end

        if syserror
            println(data.log, sysstatus.message)
            quiet || sysstatus.message != "" && println(data.alerts, sysstatus.message)
            converged = false
        end

        if converged
            # Update force/displacement and kinematic vectors for current stage
            Fin .= Fina
            U .+= ΔUa
            A .= Aa
            V .= Va

            # Backup converged state at ips
            commit_state(active_elems)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                us, fs, vs, as = components_dict[dof.name]
                dof.vals[us] = U[i]
                dof.vals[fs] = TFin[i]
                dof.vals[vs] = V[i]
                dof.vals[as] = A[i]
            end

            # Update time
            T += ΔT
            t += Δt
            data.T = T
            data.t = t
            data.residue = residue

            # Check for saving output file
            checkpoint = T > Tcheck - ΔTmin
            if checkpoint
                data.out += 1
                Tcheck += ΔTcheck # find the next output time
            end

            rstatus = update_records!(ana, checkpoint=checkpoint)
            if failed(rstatus)
                println(data.alerts, rstatus.message)
                return rstatus
            end

            if autoinc
                if ΔTbk > 0.0
                    ΔT = min(ΔTbk, Tcheck-T)
                    ΔTbk = 0.0
                else
                    if nits == 1
                        q = 1 + tanh(log10(ftol/(residue1+eps())))
                    else
                        q = 1.0
                    end

                    ΔTtr = min(q*ΔT, ΔTmax, 1-T)
                    if T + ΔTtr > Tcheck - ΔTmin
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
            data.residue = -1.0

            # Restore counters and state from last converged increment
            inc -= 1
            data.inc -= 1
            reset_state(active_elems)

            if autoinc
                println(data.log, "      increment failed")
                q = 1+tanh(log10(ftol/(residue1+eps())))
                q = clamp(q, 0.2, 0.9)
                syserror && (q = 0.7)
                ΔT = q*ΔT
                ΔT = round(ΔT, sigdigits=3)  # round to 3 significant digits
                if ΔT < ΔTmin
                    solstatus = failure("Solver did not converge.")
                    break
                end
            else
                solstatus = failure("Solver did not converge. Try `autoinc=true`. ")
                break
            end
        end

        data.ΔT = ΔT
    end

    failed(solstatus) && update_records!(ana, force=true)

    return solstatus
end
