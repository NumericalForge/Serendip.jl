# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export ThermoAnalysis, ThermoMechAnalysis, ThermoContext

const ThermoContext = Context


mutable struct ThermoMechAnalysis <: Analysis
    name::String
    model::FEModel
    data::AnalysisData

    function ThermoMechAnalysis(model::FEModel; outdir="./output", outkey="out")
        name = "Thermo-mechanical analysis"
        data = AnalysisData(outdir=outdir, outkey=outkey)
        this = new(name, model, data)

        if model.ctx.stress_state == :auto && model.ctx.ndim == 2
            model.ctx.stress_state = :plane_strain
        end
        this.data.transient = true

        return this
    end
end

const ThermoAnalysis = ThermoMechAnalysis
const ThermalAnalysis = ThermoMechAnalysis

get_natural_keys(::ThermoMechAnalysis) = [:fx, :fy, :fz, :ft]
get_essential_keys(::ThermoMechAnalysis) = [:ux, :uy, :uz, :ut]


@inline function tm_dof_map(elem)
    has_u = isdefined(@__MODULE__, :elem_map_u) && applicable(elem_map_u, elem)
    has_t = isdefined(@__MODULE__, :elem_map_t) && applicable(elem_map_t, elem)

    if has_u && has_t
        return [elem_map_u(elem); elem_map_t(elem)]
    elseif has_u
        return elem_map_u(elem)
    elseif has_t
        return elem_map_t(elem)
    end

    return Int[]
end


# Assemble the coupled thermo-mechanical tangent matrix.
function tm_mount_global_matrices(elems::AbstractVector, ndofs::Int, Δt::Float64)
    α = 1.0

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]
        RHS = zeros(ndofs)

        for elem in elems
            ty = typeof(elem)
            has_stiffness_matrix = hasmethod(elem_stiffness, (ty,))
            has_coupling_matrix = hasmethod(elem_coupling_matrix, (ty,))
            has_conductivity_matrix = hasmethod(elem_conductivity_matrix, (ty,))
            has_mass_matrix = hasmethod(elem_mass_matrix, (ty,))

            if has_stiffness_matrix
                K, rmap, cmap = elem_stiffness(elem)
                nr, nc = size(K)
                for i in 1:nr
                    for j in 1:nc
                        val = K[i, j]
                        abs(val) < eps() && continue
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, val)
                    end
                end
            end

            if has_coupling_matrix
                T0k = elem.ctx.T0 + 273.15
                Cut, rmap, cmap = elem_coupling_matrix(elem)
                nr, nc = size(Cut)
                for i in 1:nr
                    for j in 1:nc
                        val = Cut[i, j]
                        abs(val) < eps() && continue

                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, val)

                        push!(R, cmap[j])
                        push!(C, rmap[i])
                        push!(V, T0k * val)
                    end
                end
            end

            if has_conductivity_matrix
                H, rmap, cmap = elem_conductivity_matrix(elem)
                nr, nc = size(H)
                for i in 1:nr
                    for j in 1:nc
                        val = α * Δt * H[i, j]
                        abs(val) < eps() && continue
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, val)
                    end
                end

                Ut = [dof.vals[:ut] for node in elem.nodes for dof in node.dofs if dof.name == :ut]
                RHS[rmap] .-= Δt .* (H * Ut)
            end

            if has_mass_matrix
                M, rmap, cmap = elem_mass_matrix(elem)
                nr, nc = size(M)
                for i in 1:nr
                    for j in 1:nc
                        val = M[i, j]
                        abs(val) < eps() && continue
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, val)
                    end
                end
            end
        end
    end

    G = sparse(R, C, V, ndofs, ndofs)
    yield()

    return G, RHS
end


function tm_update_state(elems::AbstractVector, ΔUt::Vector{Float64}, Δt::Float64)
    ndofs = length(ΔUt)

    @withthreads begin
        ΔFin = zeros(ndofs)
        statuses = ReturnStatus[]

        for elem in elems
            ΔF, map, status = update_elem!(elem, ΔUt, Δt)
            if failed(status)
                push!(statuses, status)
                break
            end
            ΔFin[map] .+= ΔF
        end
    end

    yield()

    length(statuses) > 0 && return ΔFin, statuses[1]
    any(isnan.(ΔFin)) && return ΔFin, failure("tm_update_state: NaN values in internal forces vector")
    return ΔFin, success()
end


function tm_internal_forces(elems::AbstractVector, ndofs::Int)
    Fin = zeros(ndofs)

    for elem in elems
        map = tm_dof_map(elem)
        isempty(map) && continue

        Fe = zeros(ndofs)
        elem_internal_forces(elem, Fe)
        Fin[map] .+= Fe[map]
    end

    return Fin
end


function complete_ut_T(model::FEModel)
    haskey(model.node_fields, "ut") || return
    Ut = model.node_fields["ut"]
    T0 = model.ctx.T0

    for elem in model.elems
        elem.role == :solid || continue
        elem.shape == elem.shape.base_shape && continue
        npoints = elem.shape.npoints
        nbpoints = elem.shape.base_shape.npoints
        map = [elem.nodes[i].id for i in 1:nbpoints]
        Ute = Ut[map]
        C = elem.shape.nat_coords
        for i in nbpoints+1:npoints
            id = elem.nodes[i].id
            R = C[i, :]
            N = elem.shape.base_shape.func(R)
            Ut[id] = dot(N, Ute)
        end
    end

    model.node_fields["ut"] = Ut .+ T0
end


function stage_solver(ana::ThermoMechAnalysis, stage::Stage, solver_settings::SolverSettings; quiet=quiet)
    tol = solver_settings.tol
    ΔT0 = solver_settings.dT0
    ΔTmin = solver_settings.dTmin
    ΔTmax = solver_settings.dTmax
    rspan = solver_settings.rspan
    maxits = solver_settings.maxits
    autoinc = solver_settings.autoinc

    model = ana.model
    data = ana.data
    ctx = model.ctx
    println(data.log, "Thermo-mechanical FE analysis: Stage $(stage.id)")

    solstatus = success()

    nincs = stage.nincs
    nouts = stage.nouts
    bcs = stage.bcs
    tspan = stage.tspan
    saveouts = stage.nouts > 0
    T0 = ctx.T0
    ftol = tol

    @check tspan > 0.0 "stage_solver: stage.tspan must be positive for thermo-mechanical analysis"

    stress_state = ctx.stress_state
    ctx.ndim == 3 && @assert(stress_state == :auto)

    for elem in stage.activate
        elem.active = true
    end
    active_elems = filter(elem -> elem.active, model.elems)

    dofs, nu = configure_dofs(model, bcs)
    ndofs = length(dofs)
    umap = 1:nu
    pmap = nu+1:ndofs
    data.nu = nu

    println(data.log, "unknown dofs: $nu")
    quiet || nu == ndofs && println(data.alerts, "No essential boundary conditions")

    if stage.id == 1
        commit_state(active_elems)

        for dof in dofs
            dof.vals[dof.name] = 0.0
            dof.vals[dof.natname] = 0.0
            if dof.name == :ut
                dof.vals[:T] = T0
            end
        end

        update_records!(ana, force=true)
        complete_ut_T(model)
    end

    t = data.t
    ΔTbk = 0.0
    ΔTcheck = saveouts ? 1 / nouts : 1.0
    Tcheck = ΔTcheck

    T = 0.0
    ΔT = 1.0 / nincs
    autoinc && (ΔT = min(ΔT, ΔTmax, ΔTcheck, ΔT0))

    inc = 0
    F = zeros(ndofs)
    U = zeros(ndofs)
    R = zeros(ndofs)
    ΔFin = zeros(ndofs)
    ΔUa = zeros(ndofs)
    ΔUi = zeros(ndofs)
    Rc = zeros(ndofs)
    sysstatus = success()

    Uex = zeros(ndofs)
    Fex = zeros(ndofs)
    for bc in bcs
        compute_bc_values(ana, bc, t, Uex, Fex)
    end

    Fin = tm_internal_forces(active_elems, ndofs)
    Fex .-= Fin

    for (i, dof) in enumerate(dofs)
        U[i] = dof.vals[dof.name]
        F[i] = dof.vals[dof.natname]
    end

    local G::SparseMatrixCSC{Float64, Int64}
    local RHS::Vector{Float64}
    data.ΔT = ΔT

    while T < 1.0 - ΔTmin
        Δt = tspan * ΔT
        data.ΔT = ΔT
        data.t = t + Δt

        inc += 1
        data.inc = inc
        println(data.log, "inc $(inc)   T=$(round(T, digits=4))   ΔT=$(round(ΔT, sigdigits=4))")

        UexN = zeros(ndofs)
        FexN = zeros(ndofs)
        for bc in bcs
            compute_bc_values(ana, bc, t + Δt, UexN, FexN)
        end

        ΔUex = UexN - U
        ΔFex = FexN - F

        ΔTcr = min(rspan, 1 - T)
        αcr = min(ΔT / ΔTcr, 1.0)
        T < 1 - rspan && (ΔFex .+= αcr .* Rc)

        ΔUex[umap] .= 0.0
        ΔFex[pmap] .= 0.0

        R .= ΔFex
        ΔUa .= 0.0
        ΔUi .= ΔUex

        res = 0.0
        res1 = 0.0
        converged = false
        syserror = false

        for it in 1:maxits
            yield()

            it > 1 && (ΔUi .= 0.0)
            lastres = res

            G, RHS = tm_mount_global_matrices(active_elems, ndofs, Δt)
            Rtr = R + RHS
            sysstatus = solve_system(G, ΔUi, Rtr, nu)
            failed(sysstatus) && (syserror = true; break)

            reset_state(active_elems)
            ΔUt = ΔUa + ΔUi
            ΔFin, sysstatus = tm_update_state(active_elems, ΔUt, Δt)
            failed(sysstatus) && (syserror = true; break)

            ΔUa .+= ΔUi

            R .= ΔFex .- ΔFin
            R[pmap] .= 0.0
            res = norm(R, Inf)

            @printf(data.log, "    it %d  residue: %-10.4e\n", it, res)

            it == 1 && (res1 = res)
            res < ftol && (converged = true; break)

            isnan(res) && break
            it > 1 && res > lastres && break
        end

        if syserror
            println(data.log, sysstatus.message)
            quiet || sysstatus.message != "" && println(data.alerts, sysstatus.message)
            converged = false
        end

        if converged
            U .+= ΔUa
            F .+= ΔFin
            Uex .= UexN
            Fex .= FexN
            Rc .= (1.0 - αcr) .* Rc .+ R

            commit_state(active_elems)

            for (i, dof) in enumerate(dofs)
                dof.vals[dof.name] += ΔUa[i]
                dof.vals[dof.natname] += ΔFin[i]
                if dof.name == :ut
                    dof.vals[:T] = U[i] + T0
                end
            end

            t += Δt
            T += ΔT
            data.t = t
            data.T = T
            data.residue = res

            checkpoint = T > Tcheck - ΔTmin
            if checkpoint
                data.out += 1
                Tcheck += ΔTcheck
            end

            complete_ut_T(model)
            rstatus = update_records!(ana, checkpoint=checkpoint)
            failed(rstatus) && return rstatus

            if autoinc
                if ΔTbk > 0.0
                    ΔT = min(ΔTbk, Tcheck - T)
                    ΔTbk = 0.0
                else
                    q = 1 + tanh(log10(ftol / (res1 + eps())))
                    q = max(q, 1.1)

                    ΔTtr = min(q * ΔT, ΔTmax, 1.0 - T)
                    if T + ΔTtr > Tcheck - ΔTmin
                        ΔTbk = ΔT
                        ΔT = Tcheck - T
                    else
                        ΔT = ΔTtr
                        ΔTbk = 0.0
                    end
                end
            end
        else
            data.residue = -1.0
            inc -= 1
            data.inc = inc
            reset_state(active_elems)

            if autoinc
                println(data.log, "      increment failed")

                q = 1 + tanh(log10(ftol / (res1 + eps())))
                q = clamp(q, 0.2, 0.9)
                syserror && (q = 0.7)
                ΔT = round(q * ΔT, sigdigits=3)
                if ΔT < ΔTmin
                    solstatus = failure("Solver did not converge.")
                    break
                end
            else
                solstatus = failure("Solver did not converge. Try `autoinc=true`.")
                break
            end
        end
    end

    complete_ut_T(model)
    failed(solstatus) && update_records!(ana, force=true)

    return solstatus
end


function tm_stage_solver!(ana::ThermoMechAnalysis, stage::Stage; args...)
    kwargs = NamedTuple(args)
    quiet = get(kwargs, :quiet, false)
    solver_kwargs = kwargs
    :quiet in keys(solver_kwargs) && (solver_kwargs = Base.structdiff(solver_kwargs, (quiet=quiet,)))
    :scheme in keys(solver_kwargs) && (solver_kwargs = Base.structdiff(solver_kwargs, (scheme=solver_kwargs.scheme,)))

    return stage_solver(ana, stage, SolverSettings(; solver_kwargs...); quiet=quiet)
end
