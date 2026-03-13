# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export AcousticModalAnalysis


"""
    AcousticModalAnalysis(model::FEModel; outdir="./output", outkey="out")

Create an acoustic modal analysis for the given finite element model.
"""
mutable struct AcousticModalAnalysis<:Analysis
    name::String
    model::FEModel
    data::AnalysisData

    frequencies::Vector{Float64}
    modes::Matrix{Float64}

    function AcousticModalAnalysis(model::FEModel; outdir="./output", outkey="out")
        name = "Acoustic modal analysis"
        data = AnalysisData(outdir=outdir, outkey=outkey)
        this = new(name, model, data)

        if model.ctx.stress_state == :auto && model.ctx.ndim == 2
            model.ctx.stress_state = :plane_strain
        end
        this.data.transient = false
        this.frequencies = zeros(0)
        this.modes = zeros(0, 0)

        return this
    end
end

get_natural_keys(::AcousticModalAnalysis) = [:fx, :fy, :fz, :fq]
get_essential_keys(::AcousticModalAnalysis) = [:ux, :uy, :uz, :up]


function scale_modal_matrices(K11, M11)
    n = size(K11, 1)
    d = zeros(Float64, n)
    Kd = abs.(diag(K11))
    Md = abs.(diag(M11))

    for i in 1:n
        s = max(Kd[i], Md[i], eps())
        d[i] = 1.0 / sqrt(s)
    end

    D = Diagonal(d)
    return D * K11 * D, D * M11 * D, D
end


function lump_modal_mass(M11)
    lumps = vec(sum(M11, dims=2))
    return Diagonal(Float64.(lumps))
end


function solve_modal_eigenproblem(ana::AcousticModalAnalysis, K11, M11, nmodes::Int, eig_method::Symbol)
    Ks, Ms, D = scale_modal_matrices(K11, M11)

    if eig_method in (:auto, :arpack)
        try
            λ, V = eigs(Ks, Ms; nev=nmodes, which=:SM, ncv=40, tol=1e-4, check=1, maxiter=10000)
            eig_filter = [i for i in eachindex(λ) if isreal(λ[i]) && real(λ[i]) > 0]
            n = length(eig_filter)
            n < nmodes && println(ana.data.alerts, "Arpack found $n valid modes, expected $nmodes.")

            if n >= min(4, nmodes)
                return λ, D * V
            end
            eig_method == :arpack && error("solve_modal_eigenproblem: Arpack found only $n valid modes, expected at least $(min(4, nmodes)).")
        catch
            eig_method == :arpack && rethrow()
            println(ana.data.alerts, "Arpack failed, falling back to dense solve.")
        end
    end

    F = eigen(Matrix(Ks), Matrix(Ms))
    λ = ComplexF64.(F.values)
    V = ComplexF64.(D * F.vectors)
    return λ, V
end


function stage_solver(ana::AcousticModalAnalysis, stage::Stage, solver_settings::SolverSettings; quiet=quiet)
    nmodes = Int(solver_settings.nmodes)
    eig_method = solver_settings.eig_method

    model = ana.model
    ctx = model.ctx
    data = ana.data
    println(data.log, "Acoustic modal FE analysis")

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

    for elem in active_elems
        ty = typeof(elem)
        hasmethod(elem_mass, (ty,)) || continue
        hasmethod(elem_acoustic_mass, (ty,)) && continue
        hasmethod(elem_interface_inertia_matrix, (ty,)) && continue
        hasproperty(elem.etype, :ρ) || continue
        elem.etype.ρ == 0 && error("stage_solver: density should not be zero")
    end

    dofs, nu = configure_dofs(model, stage.bcs)
    ndofs = length(dofs)
    data.nu = nu

    println(data.log, "unknown dofs: $nu")
    quiet || nu == ndofs && println(data.alerts, "No essential boundary conditions")

    for dof in dofs
        dof.vals[dof.name] = 0.0
        dof.vals[dof.natname] = 0.0
    end

    K = am_mount_Ks(active_elems, ndofs) + am_mount_K(active_elems, ndofs) + am_mount_Kup(active_elems, ndofs)
    M = am_mount_Ms(active_elems, ndofs) + am_mount_M(active_elems, ndofs) + am_mount_Mfs(stage.bcs, ndofs) + am_mount_Mpu(active_elems, ndofs)

    K11 = K[1:nu, 1:nu]
    M11 = lump_modal_mass(M[1:nu, 1:nu])
    λ, V = solve_modal_eigenproblem(ana, K11, M11, nmodes, eig_method)

    eig_filter = [i for i in eachindex(λ) if isreal(λ[i]) && real(λ[i]) > 0]
    nfound = min(nmodes, length(eig_filter))
    eig_filter = eig_filter[sortperm(real(λ[eig_filter]))[1:nfound]]

    λ = λ[eig_filter]
    V = V[:, eig_filter]

    ω = Float64.(sqrt.(real.(λ)))
    modes = Float64.(real.(V))

    update_output_data!(model)
    save(model, joinpath(data.outdir, "$(data.outkey)-0.vtu"), quiet=true)

    data.inc = 1
    data.ΔT = 1.0

    for i in 1:nfound
        U = modes[:, i]

        for (k, dof) in enumerate(dofs[1:nu])
            dof.vals[dof.name] = U[k]
        end

        data.T = i / max(nfound, 1)
        data.out = i
        update_output_data!(model)
        save(model, joinpath(data.outdir, "$(data.outkey)-$i.vtu"), quiet=true)
    end

    for dof in dofs
        dof.vals[dof.name] = 0.0
    end

    ana.frequencies = ω
    ana.modes = modes

    return success()
end
