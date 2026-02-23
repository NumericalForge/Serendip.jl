# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechAnalysis

"""
    MechAnalysis(model::FEModel; outdir="", outkey="")

Create a static mechanical analysis for the given finite element model.

# Arguments
- `model::FEModel`: Finite element model definition.
- `outdir::String`: Directory for output files. Defaults to `"./output"`.
- `outkey::String`: Key prefix for output files. Defaults to `"out"`.

# Behavior
- Initializes an `AnalysisData` instance with the given output settings.
- If `model.ctx.stress_state == :auto` and `model.ctx.ndim == 2`, the stress state is set to `:plane_strain`.

# Example
```julia
model = FEModel(mapper) # mapper is a RegionMapper with defined mappings
analysis = MechAnalysis(model; outdir="results", outkey="job1")
```
"""
mutable struct MechAnalysis<:Analysis
    name::String
    model ::FEModel
    data::AnalysisData

    function MechAnalysis(model::FEModel; outdir="./output", outkey="out")
        name = "Mechanical analysis"
        data = AnalysisData(outdir=outdir, outkey=outkey)
        this = new(name, model, data)

        if model.ctx.stress_state==:auto && model.ctx.ndim==2
            model.ctx.stress_state = :plane_strain
        end
        this.data.transient = false

        return this
    end
end

get_natural_keys(::MechAnalysis) = [:fx, :fy, :fz]
get_essential_keys(::MechAnalysis) = [:ux, :uy, :uz]


# Assemble the global stiffness matrix
function mount_K(elems::Array{<:Element,1}, ndofs::Int )

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            Ke, rmap, cmap = elem_stiffness(elem)

            nr, nc = size(Ke)
            for i in 1:nr
                for j in 1:nc
                    val = Ke[i,j]
                    abs(val) < eps() && continue
                    # @redude R = append!(Int64[], rmap[i])
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    local K
    try
        K = sparse(R, C, V, ndofs, ndofs)
    catch err
        @show err
    end

    yield()

    return K
end


function mount_RHS(model::FEModel, ndofs::Int64, Δt::Float64)
    RHS = zeros(ndofs)
    for elem in active_elems
        F, map = elem_RHS(elem) # RHS[map] = elem_RHS(elem, Δt::Float64)
        RHS[map] = F
    end
    return RHS
end


# Get internal forces and update data at integration points
function update_state(elems::Vector{<:Element}, ΔUt::Vector{Float64}, t::Float64)
    ndofs = length(ΔUt)
    @withthreads begin
        ΔFin     = zeros(ndofs)
        statuses = ReturnStatus[]
        for elem in elems
            map = dof_map(elem)
            ΔUe = isempty(map) ? Float64[] : ΔUt[map]
            ΔF, map, status = elem_internal_forces(elem, ΔUe, t)
            if failed(status)
                push!(statuses, status)
                break
            end
            ΔFin[map] .+= ΔF
        end
    end

    yield()

    length(statuses)>0 && return ΔFin, statuses[1]
    any(isnan.(ΔFin)) && return ΔFin, failure("solve_system!: NaN values in internal forces vector")
    return ΔFin, success()
end


"""
    project_activation_equilibrium!(
        activated_cohesives,
        affected_idxs,
        lambda=0.0,
    )

Draft helper for cohesive activation transfer (not wired into `stage_solver` yet).

Builds a local patch around mobilized cohesive elements, computes force mismatch on
split-affected dofs, solves a regularized local correction problem, and applies an
incremental cohesive displacement correction through `elem_internal_forces`.

Returns `ReturnStatus`.
"""
function project_activation_equilibrium!(
    activated_cohesives::Vector{<:Element{MechCohesive}},
    affected_idxs::Vector{Int},
    lambda::Float64=0.0,
)
    lambda_factor = 1e-8
    isempty(activated_cohesives) && return success()

    # 1) Patch dofs: all displacement dofs connected to mobilized cohesive elements.
    #    Force fitting is restricted to `affected_idxs` later.
    patch_dofs = Int[]
    for elem in activated_cohesives
        append!(patch_dofs, dof_map(elem))
    end
    affected_set = Set(affected_idxs)
    patch_dofs = sort(unique(patch_dofs))
    nloc = length(patch_dofs)
    nloc == 0 && return failure("project_activation_equilibrium!: Empty local patch.")

    l2g = copy(patch_dofs)
    g2l = Dict{Int,Int}(gid => i for (i, gid) in enumerate(l2g))

    # 2) Build force target from neighboring elements.
    #    Include bulk/special elements and previously mobilized cohesive interfaces
    #    attached to activated-node stars; exclude only the cohesive elements that
    #    are being initialized in this call.
    ftarget = zeros(nloc)
    activated_set = Set(activated_cohesives)
    neigh_elems = Set{Element}()
    for cohe in activated_cohesives
        for node in cohe.nodes
            for elem in node.elems
                elem in activated_set && continue
                push!(neigh_elems, elem)
            end
        end
    end

    for elem in neigh_elems
        Fe, map, status = elem_internal_forces(elem)
        failed(status) && return status
        for (a, gid) in enumerate(map)
            gid in affected_set || continue
            iloc = get(g2l, gid, 0)
            iloc == 0 && continue
            ftarget[iloc] -= Fe[a]
        end
    end

    # 3) Assemble cohesive local stiffness Kpatch on local indexing.
    R = Int[]
    C = Int[]
    V = Float64[]
    for elem in activated_cohesives
        Ke, rmap, cmap = elem_stiffness(elem)
        nr, nc = size(Ke)
        for i in 1:nr
            iloc = get(g2l, rmap[i], 0)
            iloc == 0 && continue
            for j in 1:nc
                jloc = get(g2l, cmap[j], 0)
                jloc == 0 && continue
                val = Ke[i, j]
                abs(val) < eps() && continue
                push!(R, iloc)
                push!(C, jloc)
                push!(V, val)
            end
        end
    end
    Kpatch = sparse(R, C, V, nloc, nloc)

    # 4) Solve regularized least-squares: min ||K u - f||^2 + λ||u||^2
    KtK = Matrix(Kpatch' * Kpatch)
    rhs = Kpatch' * ftarget
    λ = lambda
    if λ <= 0.0
        trKtK = sum(diag(KtK))
        scale = trKtK > 0 ? trKtK/nloc : 1.0
        λ = lambda_factor*scale
    end
    
    Kreg = KtK + λ*Matrix(I, nloc, nloc)
    u_patch = Kreg \ rhs

    # 5) Initialize activated cohesive states from solved total local displacement.
    #    We zero stress values (other values are already zero), then update using elem_internal_forces.
    for elem in activated_cohesives
        for ip in elem.ips
            if hasfield(typeof(ip.state), :σ)
                ip.state.σ = zeros(Vec3)
                ip.cstate.σ = zeros(Vec3)
            end
            if hasfield(typeof(ip.state), :w)
                ip.state.w = zeros(Vec3)
                ip.cstate.w = zeros(Vec3)
            end
        end

        # Local element displacement vector extracted from patch solution.
        emap = dof_map(elem)
        ue = zeros(length(emap))
        for (a, gid) in enumerate(emap)
            @assert haskey(g2l, gid) "project_activation_equilibrium!: Missing gid in local map: $gid"
            ue[a] = u_patch[g2l[gid]]
        end

        _, _, status = elem_internal_forces(elem, ue, 0.0)
        
        failed(status) && return status
        
    end

    return success()
end


function stage_solver(ana::MechAnalysis, stage::Stage, solver_settings::SolverSettings; quiet=quiet)

    tol     = solver_settings.tol
    rtol    = solver_settings.rtol
    ΔT0     = solver_settings.dT0
    ΔTmin   = solver_settings.dTmin
    ΔTmax   = solver_settings.dTmax
    rspan   = solver_settings.rspan
    tangent_scheme  = solver_settings.tangent_scheme
    maxits  = solver_settings.maxits
    autoinc = solver_settings.autoinc

    model = ana.model
    data  = ana.data
    println(data.log, "Mechanical FE analysis: Stage $(stage.id)")

    solstatus = success()

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    ctx       = model.ctx
    saveouts  = stage.nouts > 0
    ftol      = tol

    stress_state = ctx.stress_state
    ctx.ndim==3 && @assert(stress_state==:auto)
    
    # Get active elements
    for elem in stage.activate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)
    
    # Check for cohesive elements
    has_cohesive_elems = any( elem -> elem isa Element{MechCohesive}, active_elems )
    cohesive_schemes = (:backward_euler, :heun, :ralston)
    has_cohesive_elems && !(tangent_scheme in cohesive_schemes) && alert("Using cohesive elements requires an implicit tangent scheme for stability $(repr(cohesive_schemes)). Current scheme: $(repr(tangent_scheme))")

    # Get dofs organized according to boundary conditions
    dofs, nu = configure_dofs(model, bcs) # unknown dofs first
    ndofs = length(dofs)
    umap  = 1:nu         # map for unknown displacements
    pmap  = nu+1:ndofs   # map for prescribed displacements
    data.nu  = nu

    println(data.log, "unknown dofs: $nu")

    quiet || nu==ndofs && println(data.alerts, "No essential boundary conditions")

    if stage.id == 1
        commit_state(active_elems) # commits changes to states from elem_init
        # Setup quantities at dofs
        for dof in dofs
            dof.vals[dof.name]    = 0.0
            dof.vals[dof.natname] = 0.0
        end

        update_records!(ana, force=true)
    end

    # Get the domain current state and backup

    # Incremental analysis
    ΔTbk    = 0.0
    ΔTcheck = saveouts ? 1/nouts : 1.0
    Tcheck  = ΔTcheck

    T  = 0.0
    ΔT = 1.0/nincs       # initial ΔT value
    autoinc && (ΔT=min(ΔT, ΔTmax, ΔTcheck, ΔT0))

    inc  = 1             # counter for the processing increment
    F    = zeros(ndofs)  # total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration
    Rc   = zeros(ndofs)  # vector of cumulated residues
    linear_domain = true
    sysstatus = ReturnStatus()

    # Get boundary conditions
    Uex = zeros(ndofs)  # essential displacements
    Fex = zeros(ndofs)  # external forces

    for bc in stage.bcs
        compute_bc_values(ana, bc, 0.0, Uex, Fex)
    end

    # Add unbalanced forces from all active elements
    for elem in active_elems
        Fe, map, status = elem_internal_forces(elem)
        failed(status) && return status
        Fex[map] .-= Fe
    end

    if tangent_scheme==:forward_euler
        p1=1.0; q11=1.0
    elseif tangent_scheme==:backward_euler
        p1=1.0; q11=1.0; a1=0.0; a2=1.0
    elseif tangent_scheme==:heun
        p1=1.0; q11=1.0; a1=0.5; a2=0.5
    elseif tangent_scheme==:ralston
        p1=2/3; q11=2/3; a1=1/4; a2=3/4
    end

    local K::SparseMatrixCSC{Float64,Int64}
    data.ΔT = ΔT

    while T < 1.0-ΔTmin

        println(data.log, "inc $(inc)  ΔT=$ΔT")

        ΔUex, ΔFex = ΔT*Uex, ΔT*Fex     # increment of external vectors

        αcr  = min(ΔT/rspan, 1.0)    # fraction of cumulated residues to apply
        ΔFex .+= αcr.*Rc # addition of residuals

        R   .= ΔFex
        ΔUa .= 0.0
        ΔUi .= ΔUex  # essential values at iteration i

        # Newton Rapshon iterations
        nits      = 0
        # err       = 0.0
        res       = 0.0
        res1      = 0.0
        converged = false
        syserror  = false

        for it in 1:maxits
            yield()

            nits += 1
            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = res # residue from last iteration

            # Predictor step
            K = mount_K(active_elems, ndofs)

            ΔUitr = p1*ΔUi
            Rtr   = q11*R
            sysstatus = solve_system!(K, ΔUitr, Rtr, nu)   # Changes unknown positions in ΔUi and R
            failed(sysstatus) && (syserror=true; break)
            # reset_state(active_elems)
            ΔUt   = ΔUa .+ ΔUitr
            ΔFin, sysstatus = update_state(active_elems, ΔUt, 0.0)
            failed(sysstatus) && (syserror=true; break)

            # Corrector step
            if tangent_scheme==:forward_euler
                ΔUi = ΔUitr
            else
                K2 = mount_K(active_elems, ndofs)
                if tangent_scheme==:backward_euler
                    K = K2
                else
                    K = a1*K + a2*K2
                end
                sysstatus = solve_system!(K, ΔUi, R, nu)   # Changes unknown positions in ΔUi and R
                failed(sysstatus) && (syserror=true; break)
                # reset_state(active_elems)
                ΔUt   = ΔUa .+ ΔUi
                ΔFin, sysstatus = update_state(active_elems, ΔUt, 0.0)
                failed(sysstatus) && (syserror=true; break)
            end

            # Update accumulated displacement
            ΔUa .+= ΔUi

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions
            res = norm(R, Inf)
            @printf(data.log, "    it %d  residue: %-10.4e\n", it, res)
            
            err = norm(ΔUi, Inf)/norm(ΔUa, Inf)
            if it==1
                res1 = res
                # res > 0.5*norm(ΔFin[umap], Inf) && (converged=false; break)
                # res1>5*ftol && (converged=false; break)
            end
            it>1  && (linear_domain=false)
            
            res<ftol && (converged=true; break)
            # err<rtol && (converged=true; break)

            isnan(res) && break
            it>maxits  && break

            it>1 && res>lastres && break
        end

        q = 0.0 # increment size factor for autoinc
        
        if syserror
            println(data.alerts, sysstatus.message)
            println(data.log, sysstatus.message)
            converged = false
        end

        if converged
            data.residue = res

            inc += 1
            data.inc = inc

            # ❱❱❱ Hook: Topology update for cohesive elements
            dof_map, new_idxs, affected_idxs, _nu, mobilized_elems = update_model_cohesive_elems(model, dofs, bcs)

            if length(affected_idxs) > 0
                nu      = _nu
                data.nu = nu
                ndofs   = length(dofs)
                umap    = 1:nu
                pmap    = nu+1:ndofs

                # Update displacements
                U   = U[dof_map]
                ΔUa = ΔUa[dof_map] 
                ΔUi = ΔUi[dof_map]
                Uex = Uex[dof_map] 

                # Re-evaluate boundary conditions
                Uex = zeros(ndofs)
                Fex = zeros(ndofs)
                for bc in stage.bcs
                    compute_bc_values(ana, bc, T, Uex, Fex)
                end
                
                # Resize F
                F = F[dof_map]
                F[new_idxs] .= 0.0 # Removes forces at new dofs. Does not perform redistribution!

                ΔFin = ΔFin[dof_map]
                ΔFin[new_idxs] .= 0.0 # Assuming that those positions are in equilibrium

                # Re-evaluate residuals
                # R  = ΔT.*Fex .- ΔFin
                # R[pmap] .= 0.0  # zero at prescribed positions
                
                R = R[dof_map]
                Rc = Rc[dof_map]
                R[new_idxs] .= 0.0
                Rc[new_idxs] .= 0.0
            end

            if length(affected_idxs) > 0 && length(mobilized_elems) > 0
                status = project_activation_equilibrium!(mobilized_elems, affected_idxs)
                if failed(status)
                    sysstatus = status
                    syserror = true
                    converged = false
                    break
                end
            end
            
            # Update forces and displacement for the current stage
            U .+= ΔUa
            F .+= ΔFin
            Rc .= (1.0-αcr).*Rc .+ R  # update cumulated residue

            # Backup converged state at ips
            commit_state(active_elems)

            # Update nodal variables at dofs
            for (i,dof) in enumerate(dofs)
                dof.vals[dof.name]    += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
            end

            # Update time
            T += ΔT
            data.T = T

            # Check for saving output file
            checkpoint = T > Tcheck - ΔTmin
            if checkpoint
                data.out += 1
                Tcheck += ΔTcheck # find the next output time
            end

            # Update analysis records
            solstatus = update_records!(ana, checkpoint=checkpoint)
            if failed(solstatus)
                println(data.alerts, solstatus.message)
                break
            end

            # Adjust increment size for next increment
            if autoinc
                if ΔTbk > 0.0
                    ΔT = min(ΔTbk/2, Tcheck-T)
                    if linear_domain
                        ΔT = min(1.5*ΔT, ΔTmax, Tcheck-T)
                    end
                    ΔTbk = 0.0
                else
                    q = 1+tanh(log10(ftol/(res1+eps())))
                    q = max(q, 1.1)

                    if linear_domain
                        q = 2
                    end
                    
                    ΔTtr = min(q*ΔT, ΔTmax, 1.0 - T)
                    
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

            reset_state(active_elems)

            if autoinc
                println(data.log, "      increment failed")

                q = 1+tanh(log10(ftol/(res1+eps())))
                q = clamp(q, 0.2, 0.9)
                syserror && (q=0.7)
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
