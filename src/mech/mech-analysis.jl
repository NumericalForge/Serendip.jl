# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechAnalysis

"""
    MechAnalysis(model::FEModel; outdir="./output", outkey="out")

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


get_natural_keys(::MechAnalysis)   = [:fx, :fy, :fz, :mx, :my, :mz]
get_essential_keys(::MechAnalysis) = [:ux, :uy, :uz, :rx, :ry, :rz]


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
    any(isnan.(ΔFin)) && return ΔFin, failure("update_state: NaN values in internal forces vector")
    return ΔFin, success()
end


const _ROTATIONAL_DOF_NAMES = (:rx, :ry, :rz)


function _mechanical_force_scale(elems::AbstractVector{<:Element}, dofs::Vector{Dof})
    ndofs = length(dofs)
    force_scale = ones(ndofs)
    rotation_length = fill(NaN, ndofs)
    is_rotation = falses(ndofs)

    for dof in dofs
        dof.name in _ROTATIONAL_DOF_NAMES && (is_rotation[dof.eq_id] = true)
    end
    any(is_rotation) || return force_scale

    for elem in elems
        map = dof_map(elem)
        any(id -> is_rotation[id], map) || continue

        diameter = cell_diameter(elem)
        isfinite(diameter) && diameter > 0.0 || throw(SerendipException(
            "stage_solver: Element $(elem.id) has rotational dofs but no positive diameter."
        ))

        for id in map
            is_rotation[id] || continue
            rotation_length[id] = isnan(rotation_length[id]) ?
                                  diameter : max(rotation_length[id], diameter)
        end
    end

    for dof in dofs
        dof.name in _ROTATIONAL_DOF_NAMES || continue
        id = dof.eq_id
        diameter = rotation_length[id]
        isfinite(diameter) && diameter > 0.0 || throw(SerendipException(
            "stage_solver: Rotational dof $(dof.name) with equation id $id has no associated element diameter."
        ))
        force_scale[id] = 1.0/diameter
    end

    return force_scale
end


function _scaled_inf_norm(values::AbstractVector, scales::AbstractVector, ids=eachindex(values))
    result = 0.0
    @inbounds for id in ids
        result = max(result, abs(scales[id]*values[id]))
    end
    return result
end


function _scaled_two_norm(values::AbstractVector, scales::AbstractVector)
    result = 0.0
    @inbounds for id in eachindex(values, scales)
        result += abs2(scales[id]*values[id])
    end
    return sqrt(result)
end


function _force_reference(Fext::AbstractVector, Fint::AbstractVector,
                          force_scale::AbstractVector, prescribed_ids)
    external_ref = _scaled_inf_norm(Fext, force_scale)
    internal_ref = _scaled_inf_norm(Fint, force_scale)
    reaction_ref = 0.0
    @inbounds for id in prescribed_ids
        reaction_ref = max(reaction_ref, abs(force_scale[id]*(Fint[id] - Fext[id])))
    end
    return max(external_ref, internal_ref, reaction_ref)
end


function _iteration_step_factor(iterations::Int, target_iterations::Int)
    iterations > 0 || throw(ArgumentError("iterations must be positive"))
    target_iterations > 0 || throw(ArgumentError("target_iterations must be positive"))
    return clamp(sqrt(target_iterations/iterations), 0.5, 1.5)
end


function stage_solver(ana::MechAnalysis, stage::Stage, solver_settings::SolverSettings; quiet=quiet)

    tol     = solver_settings.tol
    rtol    = solver_settings.rtol
    utol    = solver_settings.utol
    ΔT0     = solver_settings.dT0
    ΔTmin   = solver_settings.dTmin
    ΔTmax   = solver_settings.dTmax
    rspan   = solver_settings.rspan
    maxits  = solver_settings.maxits
    target_its = min(solver_settings.target_its, maxits)
    autoinc = solver_settings.autoinc
    lin_sol = solver_settings.linear_solver

    model = ana.model
    data  = ana.data

    solstatus = success()

    nincs     = stage.nincs
    nouts     = stage.nouts
    bcs       = stage.bcs
    constraints = stage.constraints
    ctx       = model.ctx
    saveouts  = stage.nouts > 0
    stress_state = ctx.stress_state
    ctx.ndim==3 && @assert(stress_state==:auto)
    
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
    A1, b = mount_constraint_matrices(constraints, nu, ndofs)
    has_constraints = !isempty(b)
    data.nu     = nu
    data.nfails = 0
    force_scale = _mechanical_force_scale(active_elems, dofs)

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
    Fext_acc = zeros(ndofs) # committed external force accumulated in current stage
    Fext_tr  = zeros(ndofs) # trial total external force for current stage
    Fint_tr  = zeros(ndofs) # trial total internal force for current stage
    U    = zeros(ndofs)  # total displacements for current stage
    R    = zeros(ndofs)  # vector for residuals of natural values
    ΔFin = zeros(ndofs)  # vector of internal natural values for current increment
    ΔUa  = zeros(ndofs)  # vector of essential values (e.g. displacements) for this increment
    ΔUi  = zeros(ndofs)  # vector of essential values for current iteration
    Rc   = zeros(ndofs)  # vector of cumulated residues
    Rcon = zeros(length(b))
    λ    = zeros(size(A1, 1))
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

    local K::SparseMatrixCSC{Float64,Int64}
    data.ΔT = ΔT

    while T < 1.0-ΔTmin

        println(data.log, "inc $(inc)  output $(data.out)  T=$(round(T, digits=4))   ΔT=$(round(ΔT, sigdigits=4))")

        ΔUex, ΔFex = ΔT*Uex, ΔT*Fex     # increment of external vectors

        η    = 1.0
        αcap = η*ΔT*_scaled_two_norm(F, force_scale)/(T+ΔT)/
               (_scaled_two_norm(Rc, force_scale)+eps())
        αcr  = min(ΔT/rspan, αcap, 1.0)    # fraction of cumulated residues to apply
        # αcr  = min(ΔT/rspan, 1.0)    # fraction of cumulated residues to apply

        ΔFex .+= αcr.*Rc # addition of residuals
        Fext_tr .= Fext_acc .+ ΔFex

        R   .= ΔFex
        ΔUa .= 0.0
        ΔUi .= ΔUex  # essential values at iteration i
        if has_constraints
            Rcon .= ΔT .* b
        end

        # Newton Rapshon iterations
        res   = 0.0 # residual for current iteration
        res_c = 0.0 # constraint residual
        force_ref = 0.0
        force_tol = tol
        nit = 0

        converged = false
        syserror  = false

        for it in 1:maxits
            nit = it
            yield()

            it>1 && (ΔUi.=0.0) # essential values are applied only at first iteration
            lastres = res # residue from last iteration

            K = mount_K(active_elems, ndofs)
            if !has_constraints
                sysstatus = solve_system(K, ΔUi, R, nu, lin_sol)   # Changes unknown positions in ΔUi and R
            else
                sysstatus = solve_system(K, ΔUi, R, nu, A1, Rcon, λ, lin_sol)
            end
            failed(sysstatus) && (syserror=true; break)
            ΔUt   = ΔUa .+ ΔUi
            ΔFin, sysstatus = update_state(active_elems, ΔUt, 0.0)
            failed(sysstatus) && (syserror=true; break)

            # Update accumulated displacement
            ΔUa .+= ΔUi
            Fint_tr .= F .+ ΔFin
            if it == 1
                force_ref = _force_reference(Fext_tr, Fint_tr, force_scale, pmap)
                force_tol = tol + rtol*force_ref
            end

            # Residual vector for next iteration
            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0  # zero at prescribed positions
            if has_constraints
                R[umap] .-= transpose(A1) * λ
                nu > 0 && (Rcon .-= A1 * ΔUi[umap])
                res_c = norm(Rcon, Inf)
            end
            res = _scaled_inf_norm(R, force_scale, umap)
            @printf(data.log,
                    "    it %d  residue: %-10.4e  reference: %-10.4e  limit: %-10.4e\n",
                    it, res, force_ref, force_tol)

            res <= force_tol && res_c < utol && (converged=true; break)

            !isfinite(res) && break

            it>1 && res>lastres && break
        end

        q = 0.0 # increment size factor for autoinc
        
        if syserror
            println(data.info, "$(round(T*100,digits=2))%  " * sysstatus.message)
            println(data.log, sysstatus.message)
            converged = false
        end

        if converged
            data.residue = res

            inc += 1
            data.inc = inc
            
            # Update forces and displacement for the current stage
            U .+= ΔUa
            F .= Fint_tr
            Fext_acc .= Fext_tr
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
                println(data.info, solstatus.message)
                break
            end

            # Adjust increment size for next increment
            if autoinc
                q = _iteration_step_factor(nit, target_its)
                base_ΔT = ΔTbk > 0.0 ? ΔTbk : ΔT
                ΔTbk = 0.0
                ΔTtr = min(q*base_ΔT, ΔTmax, 1.0 - T)

                if T + ΔTtr > Tcheck + ΔTmin
                    ΔTbk = ΔTtr
                    ΔT = Tcheck-T
                    @assert ΔT>=0.0
                else
                    ΔT = ΔTtr
                end
            end

        else
            data.residue = -1.0

            reset_state(active_elems)

            if autoinc
                println(data.log, "      increment failed")
                data.nfails += 1

                q = 0.5
                ΔT = q*ΔT
                ΔTbk = 0.0
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
