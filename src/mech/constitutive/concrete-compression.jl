# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export ConcreteCompression


"""
    ConcreteCompression(; E, nu=0.2, fc, epsref, Epost, rho=0.0)

Compression-only Rankine plasticity for the diffuse nonlinear response of concrete solids.
The uniaxial response is elastic up to `0.4fc`, follows a quadratic Bezier transition to
`(epsref, fc)`, and then continues with the positive tangent `Epost`. The model has no local
compressive peak or softening branch; those mechanisms are intended to be represented by
cohesive elements.

The implementation supports three-dimensional, plane-strain, and axisymmetric `MechSolid`
elements. Plane stress is not supported.
"""
mutable struct ConcreteCompression <: Constitutive
    E::Float64
    ν::Float64
    fc::Float64
    εp::Float64
    Ep::Float64
    ρ::Float64
    fe::Float64
    εe::Float64
    εq::Float64
    fq::Float64
    ε̅cp::Float64
    Bσ::Float64
    Aσ::Float64
    Hp::Float64

    function ConcreteCompression(;
        E::Real=NaN,
        nu::Real=0.2,
        fc::Real=NaN,
        epsref::Real=NaN,
        Epost::Real=NaN,
        rho::Real=0.0,
    )
        @check E > 0.0 "ConcreteCompression: Young's modulus E must be > 0. Got $(repr(E))."
        @check 0.0 <= nu < 0.5 "ConcreteCompression: Poisson's ratio nu must be in [0, 0.5). Got $(repr(nu))."
        @check fc < 0.0 "ConcreteCompression: Reference stress fc must be < 0. Got $(repr(fc))."
        @check epsref < 0.0 "ConcreteCompression: Reference strain epsref must be < 0. Got $(repr(epsref))."
        @check 0.0 < Epost < E "ConcreteCompression: Epost must satisfy 0 < Epost < E. Got $(repr(Epost))."
        @check rho >= 0.0 "ConcreteCompression: Density rho must be >= 0. Got $(repr(rho))."

        fe = 0.4*fc
        εe = fe/E
        @check epsref < εe "ConcreteCompression: epsref must be smaller than 0.4fc/E ($(repr(εe))). Got $(repr(epsref))."

        Es = (fc-fe)/(epsref-εe)
        @check Epost < Es < E "ConcreteCompression: The Bezier transition requires Epost < Es < E, where Es=$(repr(Es))."

        εq   = (fc-Epost*epsref)/(E-Epost)
        fq   = E*εq
        ε̅cp  = fc/E-epsref
        Bσ   = 2*(fq-fe)
        Aσ   = fe-2*fq+fc
        Hp   = E*Epost/(E-Epost)

        return new(E, nu, fc, epsref, Epost, rho, fe, εe, εq, fq, ε̅cp, Bσ, Aσ, Hp)
    end
end


const ConcreteCompressionProjectors = SMatrix{6,3,Float64,18}


mutable struct ConcreteCompressionState <: ConstState
    ctx::Context
    σ::Vec6
    ε::Vec6
    εp::Vec6
    ε̅c::Float64
    Δλ::Vec3
    active::Int
    projectors::ConcreteCompressionProjectors

    function ConcreteCompressionState(ctx::Context)
        return new(
            ctx,
            zeros(Vec6),
            zeros(Vec6),
            zeros(Vec6),
            0.0,
            zeros(Vec3),
            0,
            zero(ConcreteCompressionProjectors),
        )
    end
end


compat_state_type(::Type{ConcreteCompression}, ::Type{MechSolid}) = ConcreteCompressionState


function compression_strength(mat::ConcreteCompression, ε̅c::Float64)
    if ε̅c <= 0.0
        return mat.fe
    elseif ε̅c < mat.ε̅cp
        t = sqrt(ε̅c/mat.ε̅cp)
        return mat.fe + mat.Bσ*t + mat.Aσ*t^2
    else
        return mat.fc - mat.Hp*(ε̅c-mat.ε̅cp)
    end
end


function compression_plastic_modulus(mat::ConcreteCompression, ε̅c::Float64)
    if ε̅c <= eps(Float64)*mat.ε̅cp
        return Inf
    elseif ε̅c < mat.ε̅cp
        t = sqrt(ε̅c/mat.ε̅cp)
        return -(mat.Bσ + 2*mat.Aσ*t)/(2*mat.ε̅cp*t)
    else
        return mat.Hp
    end
end


function compression_projector(n::AbstractVector{<:Real})
    n1, n2, n3 = n
    return Vec6(n1^2, n2^2, n3^2, SR2*n2*n3, SR2*n1*n3, SR2*n1*n2)
end


function compression_trial_spectrum(σtr::Vec6)
    values_desc, vectors_desc = eigen(σtr)
    values = Vec3(values_desc[3], values_desc[2], values_desc[1])

    projectors = Matrix{Float64}(undef, 6, 3)
    for (a, j) in enumerate((3, 2, 1))
        projectors[:,a] = compression_projector(view(vectors_desc, :, j))
    end

    return values, ConcreteCompressionProjectors(projectors)
end


function compression_transition_root(
    mat::ConcreteCompression,
    ε̅cn::Float64,
    σ̅tr::Float64,
    Cm::Float64,
    tol::Float64,
)
    tmin = sqrt(clamp(ε̅cn/mat.ε̅cp, 0.0, 1.0))
    a = mat.Aσ-Cm*mat.ε̅cp
    b = mat.Bσ
    c = mat.fe-σ̅tr+Cm*ε̅cn

    residual(t) = (a*t+b)*t+c
    residual(tmin) >= -tol || return NaN
    residual(1.0) <= tol || return NaN

    scale = max(abs(a), abs(b), abs(c), 1.0)
    tol_t = 100*eps(Float64)
    if abs(a) <= eps(Float64)*scale
        abs(b) > eps(Float64)*scale || return NaN
        t = -c/b
        return tmin-tol_t <= t <= 1.0+tol_t ? clamp(t, tmin, 1.0) : NaN
    end

    disc = b^2-4*a*c
    disc >= -eps(Float64)*scale^2 || return NaN
    sqrt_disc = sqrt(max(disc, 0.0))
    q = -0.5*(b+copysign(sqrt_disc, b))
    roots = abs(q) > eps(Float64)*scale ? (q/a, c/q) : (-b/(2*a), -b/(2*a))

    best_t = NaN
    best_r = Inf
    for t in roots
        if tmin-tol_t <= t <= 1.0+tol_t
            tc = clamp(t, tmin, 1.0)
            r = abs(residual(tc))
            if r < best_r
                best_t = tc
                best_r = r
            end
        end
    end
    return best_t
end


function compression_return_state(
    mat::ConcreteCompression,
    ε̅cn::Float64,
    σtr_values::Vec3,
    m::Int,
    λe::Float64,
    G::Float64,
    tolσ::Float64,
)
    Cm = λe+2*G/m
    σ̅tr = sum(σtr_values[1:m])/m
    ε̅c = NaN

    if ε̅cn < mat.ε̅cp
        t = compression_transition_root(mat, ε̅cn, σ̅tr, Cm, tolσ)
        if isfinite(t)
            ε̅c = mat.ε̅cp*t^2
        end
    end

    if !isfinite(ε̅c)
        Δγ = (mat.fc-mat.Hp*(ε̅cn-mat.ε̅cp)-σ̅tr)/(mat.Hp+Cm)
        ε̅c = ε̅cn+Δγ
        ε̅c >= mat.ε̅cp-tolσ/mat.E || return NaN, Vec3(NaN, NaN, NaN)
    end

    Δγ = ε̅c-ε̅cn
    c = compression_strength(mat, ε̅c)
    Δλ = zeros(MVector{3,Float64})
    for a in 1:m
        Δλ[a] = (c-σtr_values[a]-λe*Δγ)/(2*G)
    end
    return ε̅c, Vec3(Δλ)
end


function calcD(mat::ConcreteCompression, state::ConcreteCompressionState)
    state.ctx.stress_state == :plane_stress && error("ConcreteCompression: plane stress is not supported.")
    De = calcDe(mat.E, mat.ν)
    state.active == 0 && return De

    H = compression_plastic_modulus(mat, state.ε̅c)
    isfinite(H) || return De

    m = state.active
    P = Matrix(state.projectors[:,1:m])
    Q = Matrix(De)*P
    A = P'*Q .+ H
    D = Matrix(De)-Q*(A\Q')
    return Mat6x6(D)
end


function update_state(
    mat::ConcreteCompression,
    state::ConcreteCompressionState,
    cstate::ConcreteCompressionState,
    Δε::Vector{Float64},
)
    state.ctx.stress_state == :plane_stress && error("ConcreteCompression: plane stress is not supported.")

    De = calcDe(mat.E, mat.ν)
    σtr = Vec6(cstate.σ+De*Δε)
    σtr_values, projectors = compression_trial_spectrum(σtr)
    ctrial = compression_strength(mat, cstate.ε̅c)
    stress_scale = max(abs(mat.fc), maximum(abs, σtr), 1.0)
    tolσ = 1e-9*stress_scale

    state.ε = cstate.ε+Δε
    state.projectors = projectors

    if ctrial-σtr_values[1] <= tolσ
        state.σ = σtr
        state.εp = cstate.εp
        state.ε̅c = cstate.ε̅c
        state.Δλ = zeros(Vec3)
        state.active = 0
        return state.σ-cstate.σ, success()
    end

    E, ν = mat.E, mat.ν
    G = E/(2*(1+ν))
    λe = E*ν/((1+ν)*(1-2ν))
    tolλ = tolσ/E

    for m in 1:3
        ε̅c, Δλ = compression_return_state(mat, cstate.ε̅c, σtr_values, m, λe, G, tolσ)
        isfinite(ε̅c) || continue

        Δγ = ε̅c-cstate.ε̅c
        Δγ >= -tolλ || continue
        all(Δλ[a] >= -tolλ for a in 1:m) || continue

        c = compression_strength(mat, ε̅c)
        σvalues = MVector{3,Float64}(undef)
        for a in 1:m
            σvalues[a] = c
        end
        for a in m+1:3
            σvalues[a] = σtr_values[a]+λe*Δγ
        end
        all(c-σvalues[a] <= tolσ for a in m+1:3) || continue

        Δεp = zeros(MVector{6,Float64})
        for a in 1:m
            Δεp .-= Δλ[a].*projectors[:,a]
        end

        state.σ = Vec6(projectors*Vec3(σvalues))
        state.εp = Vec6(cstate.εp+Δεp)
        state.ε̅c = max(ε̅c, cstate.ε̅c)
        state.Δλ = Δλ
        state.active = m
        return state.σ-cstate.σ, success()
    end

    return state.σ-cstate.σ, failure("ConcreteCompression: no admissible face, edge, or vertex return was found.")
end


function state_values(mat::ConcreteCompression, state::ConcreteCompressionState)
    values = stress_strain_dict(state.σ, state.ε, state.ctx.stress_state)
    values[:εpc] = state.ε̅c
    values[:fc_current] = compression_strength(mat, state.ε̅c)
    values[:active] = Float64(state.active)
    return values
end
