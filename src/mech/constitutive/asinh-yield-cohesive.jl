# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export AsinhYieldCohesive


"""
    AsinhYieldCohesive(; E, nu=0.0, fc, ft, wc=NaN, GF=NaN, psi=1.0,
                       ft_law=:hordijk, alpha=0.5, beta=0.0,
                       theta=1.0, zeta=5.0, tangent=:continuum)

Cohesive constitutive model with tensile softening and an asinh shear envelope
with fixed offset, with yield function:

`f = τ - (A(σmax) + B)*asinh(α*(σmax - σn)/ft)`

Here `σmax` is the current tensile strength and
`A(σmax) = A0*(σmax/ft)^θ`. The initial amplitude `A0` is fitted in the constructor
from a tangency condition with the compressive Mohr circle at `σmax = ft`. The fixed
offset is derived as `B = (beta*fc)/asinh(alpha*fc/ft)` and remains constant for given
`alpha`, `beta`, `fc`, and `ft`.

The model uses the tensile softening law returned by `setup_tensile_strength(ft, GF, wc, ft_law)`.
The local stiffness is obtained from `E`, `nu`, the local characteristic length `h`, and a
scaling factor `zeta`. The stiffness factor is further degraded linearly with normal opening from
`zeta` down to `0.1*zeta`. For robustness in strongly softening problems, the algorithmic tangent
limits the softening modulus to `-ft/(0.5wc)`.

The global equilibrium iterations use the continuum elastoplastic tangent by default. Set
`tangent=:consistent` to use the consistent linearization of the discrete stress update instead.

Plastic flow is non-associated in tension through `psi`, and in compression the potential has no
normal component, so plastic flow acts only in shear.

# Keyword arguments
- `E::Real`: Young's modulus (`E > 0`).
- `nu::Real=0.0`: Poisson ratio (`0 ≤ nu < 0.5`).
- `fc::Real`: Compressive strength (`fc < 0`).
- `ft::Real`: Tensile strength (`ft > 0`).
- `wc::Real=NaN`: Critical crack opening. Provide `wc > 0`, or provide `GF` instead.
- `GF::Real=NaN`: Fracture energy. Used by `setup_tensile_strength` when `wc` is not given.
- `psi::Real=1.0`: Dilatancy coefficient (`psi > 0`).
- `ft_law::Union{Symbol,AbstractSpline}=:hordijk`: Tensile softening law (`:linear`,
  `:bilinear`, `:hordijk`, or a spline).
- `alpha::Real=0.5`: Dimensionless `asinh` parameter (`alpha > 0`).
- `beta::Real=0.0`: Residual shear fraction used to derive the fixed offset `B`
  (`0 ≤ beta ≤ 1`).
- `theta::Real=1.0`: Degradation exponent of `A(σmax)` (`theta ≥ 0`).
- `zeta::Real=5.0`: Maximum elastic stiffness scaling factor (`zeta ≥ 0`).
- `tangent::Symbol=:continuum`: Tangent operator used by the global solver. Supported values are
  `:continuum` and `:consistent`.
"""
mutable struct AsinhYieldCohesive <: Constitutive
    E ::Float64
    ν ::Float64
    fc::Float64
    ft::Float64
    wc::Float64
    ψ ::Float64
    ft_law::Symbol
    ft_fun::Union{AbstractSpline,Nothing}
    α ::Float64
    β ::Float64
    B ::Float64
    θ ::Float64
    A0::Float64
    ζ ::Float64
    tangent::Symbol

    function AsinhYieldCohesive(;
        E::Real = NaN,
        nu::Real = 0.0,
        fc::Real = NaN,
        ft::Real = NaN,
        wc::Real = NaN,
        GF::Real = NaN,
        psi::Real = 1.0,
        ft_law::Union{Symbol,AbstractSpline} = :hordijk,
        alpha::Real = 0.5,
        beta::Real = 0.2,
        theta::Real = 1.0,
        zeta::Real = 5.0,
        tangent::Symbol = :continuum,
    )
        @check E > 0 "AsinhYieldCohesive: Young's modulus E must be > 0. Got $(repr(E))."
        @check 0 <= nu < 0.5 "AsinhYieldCohesive: Poisson ratio nu must be in the range [0, 0.5). Got $(repr(nu))."
        @check fc < 0 "AsinhYieldCohesive: Compressive strength fc must be < 0. Got $(repr(fc))."
        @check ft > 0 "AsinhYieldCohesive: Tensile strength ft must be > 0. Got $(repr(ft))."
        @check psi > 0 "AsinhYieldCohesive: Dilatancy coefficient psi must be positive. Got $(repr(psi))."
        @check zeta >= 0 "AsinhYieldCohesive: Factor zeta must be non-negative. Got $(repr(zeta))."
        @check alpha > 0 "AsinhYieldCohesive: alpha must be positive. Got $(repr(alpha))."
        @check 0.0 <= beta <= 1.0 "AsinhYieldCohesive: beta must be in the range [0, 1]. Got $(repr(beta))."
        @check theta >= 0 "AsinhYieldCohesive: theta must be non-negative. Got $(repr(theta))."
        @check tangent in (:continuum, :consistent) "AsinhYieldCohesive: tangent must be :continuum or :consistent. Got $(repr(tangent))."
        @check ft_law in (:linear, :bilinear, :hordijk) || ft_law isa AbstractSpline "AsinhYieldCohesive: Unknown ft_law model: $ft_law. Supported models are :linear, :bilinear, :hordijk or a custom AbstractSpline."

        B = (beta*fc)/asinh(alpha*fc/ft)

        wc, ft_law, ft_fun, status = setup_tensile_strength(ft, GF, wc, ft_law)
        failed(status) && throw(ArgumentError("AsinhYieldCohesive: " * status.message))

        A0 = fit_initial_asinh_amplitude(float(alpha), float(B), float(fc), float(ft))

        return new(E, nu, fc, ft, wc, psi, ft_law, ft_fun, alpha, beta, B, theta, A0, zeta, tangent)
    end
end


mutable struct AsinhYieldCohesiveState <: ConstState
    ctx::Context
    σ  ::Vec3
    w  ::Vec3
    up ::Float64
    Δλ ::Float64
    h  ::Float64
    kn ::Float64
    ks ::Float64
    function AsinhYieldCohesiveState(ctx::Context)
        this    = new(ctx)
        this.σ  = zeros(Vec3)
        this.w  = zeros(Vec3)
        this.up = 0.0
        this.Δλ = 0.0
        this.h  = 0.0
        this.kn = NaN
        this.ks = NaN
        return this
    end
end


compat_state_type(::Type{AsinhYieldCohesive}, ::Type{MechCohesive}) = AsinhYieldCohesiveState


function calc_total_A_from_t(t::Float64, α::Float64, σmax::Float64, fc::Float64, ft::Float64)
    σT  = fc/2 * (1 - cos(t))
    τT  = -fc/2 * sin(t)
    x   = (σmax - σT)/ft
    den = asinh(α*x)
    den > eps() || return 0.0
    return τT/den
end


function calc_A_from_t(t::Float64, α::Float64, B::Float64, σmax::Float64, fc::Float64, ft::Float64)
    return calc_total_A_from_t(t, α, σmax, fc, ft) - B
end


function tangency_residual(t::Float64, α::Float64, σmax::Float64, fc::Float64, ft::Float64)
    σT = fc/2 * (1 - cos(t))
    x  = (σmax - σT)/ft
    C  = calc_total_A_from_t(t, α, σmax, fc, ft)
    return (C*α)/(ft*sqrt(1 + α^2*x^2)) - cot(t)
end


function fit_initial_asinh_amplitude(α::Float64, B::Float64, fc::Float64, ft::Float64)
    a = 0.1*pi
    b = 0.5*pi
    f(t) = tangency_residual(t, α, ft, fc, ft)
    t, status = findroot_bisection(f, a, b, 1e-8, 1e-8)
    failed(status) && throw(ArgumentError("AsinhYieldCohesive: " * status.message))
    A0 = calc_A_from_t(t, α, B, ft, fc, ft)
    A0 > 0 || error("AsinhYieldCohesive: invalid initial asinh amplitude A0 = $A0.")
    return A0
end


function calc_∂A∂σmax(mat::AsinhYieldCohesive, σmax::Float64)
    σmax <= 0.0 && return 0.0
    A = mat.A0*(σmax/mat.ft)^mat.θ
    return mat.θ*A/σmax
end


function calc_σmax(mat::AsinhYieldCohesive, up::Float64)
    return calc_tensile_strength(mat, up)
end


function deriv_σmax_up(mat::AsinhYieldCohesive, up::Float64)
    return calc_tensile_strength_derivative(mat, up)
end


function yield_shear_envelope(mat::AsinhYieldCohesive, σn::Float64, σmax::Float64)
    x = (σmax - σn)/mat.ft
    A = mat.A0*(σmax/mat.ft)^mat.θ
    return (A + mat.B)*asinh(mat.α*x)
end


function yield_func(mat::AsinhYieldCohesive, σ::Vec3, σmax::Float64)
    σn, τ1, τ2 = σ
    τ = sqrt(τ1^2 + τ2^2)
    return τ - yield_shear_envelope(mat, σn, σmax)
end


function stress_strength_ratio(mat::AsinhYieldCohesive, σ::AbstractVector)
    σn, τ1, τ2 = σ

    τ    = sqrt(τ1^2 + τ2^2)
    σmax = calc_σmax(mat, 0.0)
    τmax = yield_shear_envelope(mat, σn, σmax)
    ησ   = σmax > 0 ? σn/σmax : 0.0
    ητ   = τmax > 0 ? τ/τmax : (τ > 0 ? Inf : 0.0)
    return max(ησ, ητ)
end


function cap_stress(mat::AsinhYieldCohesive, σ::AbstractVector)
    σn, τ1, τ2 = σ
    σmax = calc_σmax(mat, 0.0)
    return Vec3(min(σn, σmax), τ1, τ2)
end


function yield_derivs(mat::AsinhYieldCohesive, σ::Vec3, σmax::Float64)
    σn, τ1, τ2 = σ

    τ       = sqrt(τ1^2 + τ2^2)
    x       = (σmax - σn)/mat.ft
    A       = mat.A0*(σmax/mat.ft)^mat.θ
    C       = A + mat.B
    ∂A∂σmax = calc_∂A∂σmax(mat, σmax)
    ∂asinh  = (C*mat.α)/(mat.ft*sqrt(1 + mat.α^2*x^2))

    ∂f∂σn = ∂asinh
    ∂f∂σmax = -∂A∂σmax*asinh(mat.α*x) - ∂asinh

    if τ > 0.0
        ∂f∂σ = Vec3(∂f∂σn, τ1/τ, τ2/τ)
    else
        ∂f∂σ = Vec3(∂f∂σn, 0.0, 0.0)
    end
    return ∂f∂σ, ∂f∂σmax
end


function potential_derivs(mat::AsinhYieldCohesive, σ::Vec3)
    σn, τ1, τ2 = σ
    if σn < 0.0
        return Vec3(0.0, τ1, τ2)
    else
        return Vec3(mat.ψ^2*σn, τ1, τ2)
    end
end


function calc_kn_ks(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState)
    ζmax = mat.ζ
    ζmin = 0.1*ζmax

    wn = state.w[1]
    w0 = 0.3*mat.wc
    ζ  = clamp(ζmax - (ζmax - ζmin)*wn/w0, ζmin, ζmax)

    kn = mat.E*ζ/state.h
    G  = mat.E/(2*(1 + mat.ν))
    ks = G*ζ/state.h
    return kn, ks
end


function calc_tangent_softening_modulus(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState)
    H = deriv_σmax_up(mat, state.up)
    Hcap = -mat.ft/(0.5*mat.wc)
    return max(H, Hcap)
end


function calcD(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState)
    σmax = calc_σmax(mat, state.up)
    kn, ks = if isnan(state.kn) || isnan(state.ks)
        calc_kn_ks(mat, state)
    else
        state.kn, state.ks
    end
    tiny = 1e-6*mat.ft

    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    if state.Δλ == 0.0
        return De
    elseif σmax <= tiny && state.w[1] >= 0.0
        return De*1e-6
    elseif mat.tangent == :continuum
        n, ∂f∂σmax = yield_derivs(mat, state.σ, σmax)
        m = potential_derivs(mat, state.σ)
        H = calc_tangent_softening_modulus(mat, state)
        De_m = De*m
        nT_De = n'*De
        den = dot(n, De_m) - ∂f∂σmax*H*norm(m)
        return De - (De_m*nT_De)/den
    else
        σn, τ1, τ2 = state.σ
        q = σn < 0.0 ? 0.0 : mat.ψ^2
        sn = 1.0/(1.0 + state.Δλ*kn*q)
        ss = 1.0/(1.0 + state.Δλ*ks)
        S = @SMatrix [ sn   0.0  0.0
                       0.0  ss   0.0
                       0.0  0.0  ss ]

        n, ∂f∂σmax = yield_derivs(mat, state.σ, σmax)
        m = Vec3(q*σn, τ1, τ2)
        norm_m = norm(m)
        unit_m = m/(norm_m + eps())
        a = Vec3(-sn*kn*q*σn, -ss*ks*τ1, -ss*ks*τ2)
        c = state.Δλ*Vec3(unit_m[1]*q*sn, unit_m[2]*ss, unit_m[3]*ss)
        γ = norm_m + state.Δλ*dot(unit_m, Vec3(q*a[1], a[2], a[3]))
        H = calc_tangent_softening_modulus(mat, state)
        b = S*n + ∂f∂σmax*H*c
        rλ = dot(n, a) + ∂f∂σmax*H*γ
        return (S - a*b'/rλ)*De
    end
end


function plastic_update(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, cstate::AsinhYieldCohesiveState, σtr::Vec3)
    kn, ks = state.kn, state.ks
    σntr, τ1tr, τ2tr = σtr

    maxits    = 50
    converged = false
    Δλ        = σntr < 0 ? 0.0 : max(cstate.Δλ, 0.0)
    up        = cstate.up
    σ         = cstate.σ
    tol       = mat.ft*1e-8
    σtol      = mat.ft*1e-6
    σ0        = copy(σ)
    Δλ_pos    = 0.0
    Δλ_neg    = NaN

    for _ in 1:maxits
        den_σn = 1.0 + Δλ*kn*mat.ψ^2
        den_τ  = 1.0 + Δλ*ks

        σn = σntr < 0 ? σntr : σntr/den_σn
        τ1 = τ1tr/den_τ
        τ2 = τ2tr/den_τ
        σ  = Vec3(σn, τ1, τ2)

        m      = potential_derivs(mat, σ)
        norm_m = norm(m)
        unit_m = m / (norm_m + eps())

        up   = cstate.up + Δλ*norm_m
        σmax = calc_σmax(mat, up)
        H    = deriv_σmax_up(mat, up)

        f = yield_func(mat, σ, σmax)
        if f >= 0.0
            Δλ_pos = Δλ
        else
            Δλ_neg = Δλ
        end

        if abs(f) < tol && maximum(abs, σ-σ0) < σtol
            converged = true
            break
        end

        σ0 = copy(σ)

        if σntr < 0
            ∂σn∂Δλ = 0.0
            ∂m∂Δλ  = Vec3(0.0, -τ1tr*ks/den_τ^2, -τ2tr*ks/den_τ^2)
        else
            ∂σn∂Δλ = -σntr*kn*mat.ψ^2/den_σn^2
            ∂m∂Δλ  = Vec3(-σntr*kn*mat.ψ^4/den_σn^2, -τ1tr*ks/den_τ^2, -τ2tr*ks/den_τ^2)
        end

        ∂σ∂Δλ    = Vec3(∂σn∂Δλ, -τ1tr*ks/den_τ^2, -τ2tr*ks/den_τ^2)
        ∂up∂Δλ   = norm_m + Δλ*dot(unit_m, ∂m∂Δλ)
        ∂σmax∂Δλ = H*∂up∂Δλ

        ∂f∂σ, ∂f∂σmax = yield_derivs(mat, σ, σmax)
        ∂f∂Δλ = dot(∂f∂σ, ∂σ∂Δλ) + ∂f∂σmax*∂σmax∂Δλ
        abs(∂f∂Δλ) < eps() && break

        Δλ_new = Δλ - f/∂f∂Δλ

        if isfinite(Δλ_neg)
            if Δλ_new < Δλ_neg || Δλ_new > Δλ_pos
                Δλ_lo = min(Δλ_pos, Δλ_neg)
                Δλ_hi = max(Δλ_pos, Δλ_neg)
                Δλ_new = 0.5*(Δλ_lo + Δλ_hi)
            end
        elseif !isfinite(Δλ_new) || Δλ_new < 0.0
            Δλ_new = 0.5*Δλ
        end

        Δλ = max(Δλ_new, 0.0)
    end

    if converged
        state.σ  = σ
        state.Δλ = Δλ
        state.up = up
        return success()
    end

    return failure("AsinhYieldCohesive: plastic update failed.")
end


function update_state(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, cstate::AsinhYieldCohesiveState, Δw::Vector{Float64})
    kn, ks = calc_kn_ks(mat, state)
    state.kn = kn
    state.ks = ks
    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    σmax = calc_σmax(mat, cstate.up)

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("AsinhYieldCohesive: Invalid value for relative displacement: Δw = $Δw")
    end

    σtr = cstate.σ + De*Δw
    ftr = yield_func(mat, σtr, σmax)

    if σmax == 0.0 && cstate.w[1] + Δw[1] >= 0.0
        state.σ   = Vec3(0.0, 0.0, 0.0)
        state.Δλ  = 1.0
        state.up  = cstate.up + norm(Δw)
    elseif ftr <= 0.0
        state.Δλ = 0.0
        state.σ  = σtr
    else
        status = plastic_update(mat, state, cstate, σtr)
        failed(status) && return state.σ, status
    end

    state.w = cstate.w + Δw
    return state.σ - cstate.σ, success()
end


function state_values(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState)
    σmax = calc_σmax(mat, state.up)
    σn, τ1, τ2 = state.σ
    τ = sqrt(τ1^2 + τ2^2)

    return Dict(
        :w    => state.w[1],
        :σn   => σn,
        :τ    => τ,
        :up   => state.up,
        :σmax => σmax
    )
end


function output_keys(mat::AsinhYieldCohesive)
    return Symbol[:w, :σn, :τ, :up, :σmax]
end
