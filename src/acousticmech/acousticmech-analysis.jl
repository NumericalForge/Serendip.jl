# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export AcousticAnalysis, AcousticMechAnalysis, AcousticContext

const AcousticContext = Context


"""
    AcousticMechAnalysis(model::FEModel; outdir="./output", outkey="out")

Create a transient acoustic analysis for the given finite element model.
"""
mutable struct AcousticMechAnalysis<:Analysis
    name::String
    model::FEModel
    data::AnalysisData

    function AcousticMechAnalysis(model::FEModel; outdir="./output", outkey="out")
        name = "Acoustic analysis"
        data = AnalysisData(outdir=outdir, outkey=outkey)
        this = new(name, model, data)

        if model.ctx.stress_state == :auto && model.ctx.ndim == 2
            model.ctx.stress_state = :plane_strain
        end
        this.data.transient = true

        return this
    end
end

const AcousticAnalysis = AcousticMechAnalysis

get_natural_keys(::AcousticMechAnalysis) = [:fx, :fy, :fz, :fq]
get_essential_keys(::AcousticMechAnalysis) = [:ux, :uy, :uz, :up]


is_acoustic_fs_bc(bc::BoundaryCondition) = bc.kind in (:face, :edge) && get(bc.conds, :fs, false) == true


function validate_acoustic_bcs(stage::Stage, ctx::Context)
    facet_types = Dict{Int,Symbol}()

    for bc in stage.bcs
        bctype =
            haskey(bc.conds, :up) ? :pressure :
            haskey(bc.conds, :tq) ? :flux :
            is_acoustic_fs_bc(bc) ? :fs :
            :other

        bctype == :other && continue
        bctype == :fs && @check ctx.g > 0.0 "validate_acoustic_bcs: `fs=true` requires positive gravity `g`."

        for facet in bc.target
            if bctype == :fs
                has_up_dofs = all(node -> get_dof(node, :up) !== nothing, facet.nodes)
                has_up_dofs || error("validate_acoustic_bcs: `fs=true` can only be applied to facets whose nodes have `:up` dofs (facet $(facet.id)).")
            end

            old = get(facet_types, facet.id, :none)
            if old != :none && old != bctype
                error("validate_acoustic_bcs: conflicting acoustic BCs on facet $(facet.id): $old and $bctype")
            end
            facet_types[facet.id] = bctype
        end
    end

    return nothing
end


function am_mount_M(elems::Array{<:Element,1}, ndofs::Int)
    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_acoustic_mass, (ty,)) || continue

            Me, rmap, cmap = elem_acoustic_mass(elem)

            nr, nc = size(Me)
            for i in 1:nr
                for j in 1:nc
                    val = Me[i, j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    M = sparse(R, C, V, ndofs, ndofs)
    yield()

    return M
end


function acoustic_free_surface_mass(elem::Element, facet::Cell)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    g      = elem.ctx.g
    nnodes = length(facet.nodes)
    C      = get_coords(facet.nodes, ndim)
    M      = zeros(nnodes, nnodes)

    @check g > 0.0 "acoustic_free_surface_mass: `fs=true` requires positive gravity `g`."

    for ip in get_ip_coords(facet.shape)
        R = ip.coord
        w = ip.w
        N = facet.shape.func(R)
        D = facet.shape.deriv(R)
        J = C' * D
        X = C' * N

        coef = 0.0
        if ndim == 2
            elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * X[1])
            coef = norm(J) * w * th / g
        else
            coef = norm(cross(J[:, 1], J[:, 2])) * w / g
        end

        M .+= coef * N * N'
    end

    map = Int[get_dof(node, :up).eq_id for node in facet.nodes]
    return M, map, map
end


function am_mount_Mfs(bcs::Vector{BoundaryCondition}, ndofs::Int)
    isempty(bcs) && return spzeros(ndofs, ndofs)

    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]
        seen = Set{Int}()

        for bc in bcs
            is_acoustic_fs_bc(bc) || continue

            for facet in bc.target
                facet.id in seen && continue
                push!(seen, facet.id)

                Me, rmap, cmap = acoustic_free_surface_mass(facet.owner, facet)
                nr, nc = size(Me)
                for i in 1:nr
                    for j in 1:nc
                        val = Me[i, j]
                        abs(val) < eps() && continue
                        push!(R, rmap[i])
                        push!(C, cmap[j])
                        push!(V, val)
                    end
                end
            end
        end
    end

    M = sparse(R, C, V, ndofs, ndofs)
    yield()
    return M
end


function am_mount_K(elems::Array{<:Element,1}, ndofs::Int)
    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_acoustic_stiffness, (ty,)) || continue

            Ke, rmap, cmap = elem_acoustic_stiffness(elem)

            nr, nc = size(Ke)
            for i in 1:nr
                for j in 1:nc
                    val = Ke[i, j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    K = sparse(R, C, V, ndofs, ndofs)
    yield()

    return K
end


function am_mount_Ks(elems::Array{<:Element,1}, ndofs::Int)
    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_stiffness, (ty,)) || continue
            hasmethod(elem_acoustic_stiffness, (ty,)) && continue
            hasmethod(elem_coupling_matrix, (ty,)) && continue

            Ke, rmap, cmap = elem_stiffness(elem)

            nr, nc = size(Ke)
            for i in 1:nr
                for j in 1:nc
                    val = Ke[i, j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    K = sparse(R, C, V, ndofs, ndofs)
    yield()

    return K
end


function am_mount_Ms(elems::Array{<:Element,1}, ndofs::Int)
    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_mass, (ty,)) || continue
            hasmethod(elem_acoustic_mass, (ty,)) && continue
            hasmethod(elem_interface_inertia_matrix, (ty,)) && continue

            Me, rmap, cmap = elem_mass(elem)

            nr, nc = size(Me)
            for i in 1:nr
                for j in 1:nc
                    val = Me[i, j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    M = sparse(R, C, V, ndofs, ndofs)
    yield()

    return M
end


function am_mount_Kup(elems::Array{<:Element,1}, ndofs::Int)
    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_coupling_matrix, (ty,)) || continue

            Ce, pmap, umap = elem_coupling_matrix(elem)

            nr, nc = size(Ce)
            for i in 1:nr
                for j in 1:nc
                    val = Ce[i, j]
                    abs(val) < eps() && continue
                    push!(R, umap[j])
                    push!(C, pmap[i])
                    push!(V, val)
                end
            end
        end
    end

    K = sparse(R, C, V, ndofs, ndofs)
    yield()

    return K
end


function am_mount_Mpu(elems::Array{<:Element,1}, ndofs::Int)
    @withthreads begin
        R, C, V = Int64[], Int64[], Float64[]

        for elem in elems
            ty = typeof(elem)
            hasmethod(elem_interface_inertia_matrix, (ty,)) || continue

            Me, rmap, cmap = elem_interface_inertia_matrix(elem)

            nr, nc = size(Me)
            for i in 1:nr
                for j in 1:nc
                    val = Me[i, j]
                    abs(val) < eps() && continue
                    push!(R, rmap[i])
                    push!(C, cmap[j])
                    push!(V, val)
                end
            end
        end
    end

    M = sparse(R, C, V, ndofs, ndofs)
    yield()

    return M
end


function stage_solver(ana::AcousticMechAnalysis, stage::Stage, solver_settings::SolverSettings; quiet=quiet)
    tol     = solver_settings.tol
    ΔT0     = solver_settings.dT0
    ΔTmin   = solver_settings.dTmin
    ΔTmax   = solver_settings.dTmax
    maxits  = solver_settings.maxits
    autoinc = solver_settings.autoinc

    model = ana.model
    ctx = model.ctx
    data = ana.data
    println(data.log, "Acoustic FE analysis: Stage $(stage.id)")

    solstatus = success()

    nincs    = stage.nincs
    nouts    = stage.nouts
    bcs      = stage.bcs
    tspan    = stage.tspan
    saveouts = stage.nouts > 0

    @check tspan > 0.0 "stage_solver: stage.tspan must be positive for acoustic analysis"

    stress_state = ctx.stress_state
    ctx.ndim == 3 && @check stress_state == :d3
    validate_acoustic_bcs(stage, ctx)

    for elem in stage.activate
        elem.active = true
    end
    for elem in stage.deactivate
        elem.active = false
    end
    active_elems = filter(elem -> elem.active, model.elems)

    dofs, nu = configure_dofs(model, bcs)
    ndofs = length(dofs)
    umap = 1:nu
    pmap = nu+1:ndofs
    data.nu = nu

    println(data.log, "unknown dofs: $nu")
    quiet || nu == ndofs && println(data.alerts, "No essential boundary conditions")

    components_dict = Dict(
        :up => (:up, :fq, :vp, :ap),
        :ux => (:ux, :fx, :vx, :ax),
        :uy => (:uy, :fy, :vy, :ay),
        :uz => (:uz, :fz, :vz, :az),
        :rx => (:rx, :mx, :vrx, :arx),
        :ry => (:ry, :my, :vry, :ary),
        :rz => (:rz, :mz, :vrz, :arz),
    )

    if stage.id == 1
        commit_state(active_elems)

        for dof in dofs
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[us] = 0.0
            dof.vals[fs] = 0.0
            dof.vals[vs] = 0.0
            dof.vals[as] = 0.0
        end

        A0 = zeros(ndofs)
        V0 = zeros(ndofs)
        Uex0 = zeros(ndofs)
        Fex0 = zeros(ndofs)
        for bc in bcs
            compute_bc_values(ana, bc, data.t, Uex0, Fex0)
        end

        M0 = am_mount_Ms(active_elems, ndofs) + am_mount_M(active_elems, ndofs) + am_mount_Mfs(bcs, ndofs) + am_mount_Mpu(active_elems, ndofs)
        sysstatus = solve_system!(M0, A0, Fex0, nu)
        failed(sysstatus) && return sysstatus

        for (i, dof) in enumerate(dofs)
            us, fs, vs, as = components_dict[dof.name]
            dof.vals[vs] = V0[i]
            dof.vals[as] = A0[i]
            dof.vals[fs] = Fex0[i]
        end

        update_records!(ana, force=true)
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
    ΔUi = zeros(ndofs)
    ΔUk = zeros(ndofs)

    ΔUexi = zeros(ndofs)
    Uexi = zeros(ndofs)
    Uexi_prev = zeros(ndofs)
    Fexi_prev = zeros(ndofs)

    A = zeros(ndofs)
    V = zeros(ndofs)
    for (i, dof) in enumerate(dofs)
        us, fs, vs, as = components_dict[dof.name]
        U[i] = dof.vals[us]
        F[i] = dof.vals[fs]
        V[i] = dof.vals[vs]
        A[i] = dof.vals[as]
    end

    data.ΔT = ΔT

    while T < 1.0 - ΔTmin
        Δt = tspan * ΔT
        data.ΔT = ΔT
        data.t = t + Δt

        inc += 1
        data.inc = inc
        println(data.log, "inc $(inc)   T=$(round(T, digits=4))   ΔT=$(round(ΔT, sigdigits=4))")

        Uexi .= 0.0
        Fexi = zeros(ndofs)
        for bc in bcs
            compute_bc_values(ana, bc, t + Δt, Uexi, Fexi)
        end

        ΔUexi .= Uexi .- Uexi_prev
        ΔF = Fexi .- Fexi_prev

        ΔUexi[umap] .= 0.0
        ΔF[pmap] .= 0.0

        ΔUi .= 0.0
        ΔUk .= ΔUexi

        res = 0.0
        res1 = 0.0
        converged = false
        syserror = false
        local At = copy(A)
        local Vt = copy(V)
        local sysstatus = success()

        for it in 1:maxits
            yield()

            it > 1 && (ΔUk .= 0.0)
            lastres = res

            K = am_mount_Ks(active_elems, ndofs) + am_mount_K(active_elems, ndofs) + am_mount_Kup(active_elems, ndofs)
            M = am_mount_Ms(active_elems, ndofs) + am_mount_M(active_elems, ndofs) + am_mount_Mfs(bcs, ndofs) + am_mount_Mpu(active_elems, ndofs)

            Kdyn = K + 4 / Δt^2 * M
            ΔFdyn = ΔF + M * (A + 4 * V / Δt - 4 * ΔUi / Δt^2)

            sysstatus = solve_system!(Kdyn, ΔUk, ΔFdyn, nu)
            failed(sysstatus) && (syserror = true; break)

            ΔUit = ΔUi + ΔUk
            ΔFin .= 0.0
            ΔFin, sysstatus = update_state(active_elems, ΔUit, Δt)
            failed(sysstatus) && (syserror = true; break)

            Vt = -V + 2 / Δt * ΔUit
            At = -A + 4 / Δt^2 * (ΔUit - V * Δt)

            R .= ΔF .- (ΔFin + M * At)
            res = maximum(abs, R[umap])

            ΔUi .+= ΔUk

            @printf(data.log, "    it %d  residue: %-10.4e\n", it, res)

            it == 1 && (res1 = res)
            res < tol && (converged = true; break)
            isnan(res) && break
            it > 1 && res > lastres && break
        end

        if syserror
            println(data.log, sysstatus.message)
            quiet || sysstatus.message != "" && println(data.alerts, sysstatus.message)
            converged = false
        end

        if converged
            U .+= ΔUi
            F .+= ΔFin
            Fexi_prev .= Fexi
            Uexi_prev .= Uexi

            A .= At
            V .= Vt
            commit_state(active_elems)

            for (i, dof) in enumerate(dofs)
                us, fs, vs, as = components_dict[dof.name]
                dof.vals[us] = U[i]
                dof.vals[fs] = F[i]
                dof.vals[vs] = V[i]
                dof.vals[as] = A[i]
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

            rstatus = update_records!(ana, checkpoint=checkpoint)
            failed(rstatus) && return rstatus

            if autoinc
                if ΔTbk > 0.0
                    ΔT = min(ΔTbk, Tcheck - T)
                    ΔTbk = 0.0
                else
                    q = res1 > 0 ? 1.0 + tanh(log10(tol / res1)) : 2.0
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
                q = res1 > 0 ? 1 + tanh(log10(tol / res1)) : 0.5
                q = clamp(q, 0.2, 0.9)
                syserror && (q = 0.7)
                ΔT = round(q * ΔT, sigdigits=3)
                if ΔT < ΔTmin
                    solstatus = failure("Solver did not converge.")
                    break
                end
            else
                solstatus = failure("Solver did not converge.")
                break
            end
        end
    end

    failed(solstatus) && update_records!(ana, force=true)

    return solstatus
end
