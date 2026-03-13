 #This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MohrCoulombCohesive


"""
    MohrCoulombCohesive(; E, nu=0.0, ft, GF, wc, mu, psi=mu, ft_law=:hordijk, zeta=5.0)

Constitutive model for cohesive elements with a Mohr–Coulomb strength criterion.
The tensile branch is regularized with the bulk characteristic length `h` to
ensure mesh-objective dissipation. Normal and shear interface stiffnesses are
derived from `E`, `nu`, `zeta`, and `h`.

# Keyword arguments
- `E::Real`
  Young’s modulus of the bulk material (> 0).
- `nu::Real`
  Poisson’s ratio (0 ≤ ν < 0.5).
- `ft::Real`
  Tensile strength (> 0).
- `wc::Real`
  Critical crack opening (> 0 if provided). May be computed from `GF`.
- `GF::Real`
  Mode-I fracture energy (> 0 if provided). May be used to compute `wc`.
- `mu::Real`
  Friction coefficient (> 0).
- `ψ::Real`
  Dilarancy coefficient (> 0)
- `ft_law::Union{Symbol,AbstractSpline} = :hordijk`
  Tensile softening law. Symbols: `:linear`, `:bilinear`, `:hordijk`; or a custom Spline.
- `zeta::Real = 5.0`
  Dimensionless factor controlling elastic relative displacements (≥ 0).

# Returns
A `MohrCoulombCohesive` object.

# Notes
- Provide either `wc` or `GF`. If only `GF` is given, `wc` is computed from `ft_law`.
- Frictional strength is governed by `mu`.
"""
mutable struct MohrCoulombCohesive<:Constitutive
    E  ::Float64
    ν  ::Float64
    ft ::Float64
    wc ::Float64
    μ  ::Float64
    ψ  ::Float64
    ft_law::Symbol
    ft_fun::Union{AbstractSpline,Nothing}
    ζ  ::Float64

    function MohrCoulombCohesive(; 
        E::Real  = NaN,
        nu::Real = 0.0,
        ft::Real = NaN,
        wc::Real = NaN,
        GF::Real = NaN,
        mu::Real = NaN,
        psi::Real = NaN,
        ft_law::Union{Symbol,AbstractSpline} = :hordijk,
        zeta::Real=5.0
    )
        @check E>0 "MohrCoulombCohesive: Young's modulus E must be > 0. Got $(repr(E))."
        @check 0<=nu<0.5 "MohrCoulombCohesive: Poisson ratio nu must be in the range [0, 0.5). Got $(repr(nu))."
        @check ft>0 "MohrCoulombCohesive: Tensile strength ft must be > 0. Got $(repr(ft))."
        @check mu>0 "MohrCoulombCohesive: Friction coefficient mu must be non-negative. Got $(repr(mu))."

        isnan(psi) && (psi = mu)
        @check psi>0 "MohrCoulombCohesive: Dilatancy coefficient psi must be non-negative. Got $(repr(psi))."
        @check zeta>=0 "MohrCoulombCohesive: Factor zeta must be non-negative. Got $(repr(zeta))."

        wc, ft_law, ft_fun, status = setup_tensile_strength(ft, GF, wc, ft_law)
        failed(status) && throw(ArgumentError("MohrCoulombCohesive: " * status.message))

        return new(E, nu, ft, wc, mu, psi, ft_law, ft_fun, zeta)
    end
end


mutable struct MohrCoulombCohesiveState<:ConstState
    ctx::Context
    σ  ::Vec3        # stress
    w  ::Vec3        # relative displacements
    up ::Float64     # effective plastic relative displacement
    Δλ ::Float64     # plastic multiplier
    h  ::Float64     # characteristic length from bulk elements
    function MohrCoulombCohesiveState(ctx::Context)
        this    = new(ctx)
        this.σ  = zeros(Vec3)
        this.w  = zeros(Vec3)
        this.up = 0.0
        this.Δλ = 0.0
        this.h  = 0.0
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{MohrCoulombCohesive}, ::Type{MechCohesive}) = MohrCoulombCohesiveState


function yield_func(mat::MohrCoulombCohesive, σ::Vec3, σmax::Float64)
    σn, τ1, τ2 = σ
    return √(τ1^2 + τ2^2) + (σn - σmax)*mat.μ
end


function stress_strength_ratio(mat::MohrCoulombCohesive, σ::AbstractVector)
    σn, τ1, τ2 = σ

    σmax = calc_σmax(mat, 0.0)
    τmax = (σmax - σn)*mat.μ
    τ    = √(τ1^2 + τ2^2) + eps()
    return max(σn/σmax, τ/τmax)
end


function cap_stress(mat::MohrCoulombCohesive, σ::AbstractVector)
    σn, τ1, τ2 = σ
    σmax = calc_σmax(mat, 0.0)
    τmax = max( (σmax - σn)*mat.μ, 0.0)
    return Vec3(min(σn, σmax), min(τ1, τmax), min(τ2, τmax))
end


function yield_derivs(mat::MohrCoulombCohesive, σ::Vec3)
    σn, τ1, τ2 = σ
    τ = √(τ1^2 + τ2^2) + eps()

    ∂f∂σmax = -mat.μ
    return Vec3( mat.μ, τ1/τ, τ2/τ ), ∂f∂σmax
end


function potential_derivs(mat::MohrCoulombCohesive, σ::Vec3)
    σn, τ1, τ2 = σ

    if σn < 0.0 
        return Vec3( 0.0, τ1, τ2 )
    else
        return Vec3( mat.ψ^2*σn + eps(), τ1, τ2 )
    end
end


function calc_σmax(mat::MohrCoulombCohesive, up::Float64)
    return calc_tensile_strength(mat, up)
end


function deriv_σmax_up(mat::MohrCoulombCohesive, up::Float64)
    # ∂σmax/∂up
    return calc_tensile_strength_derivative(mat, up)
end


function calc_kn_ks(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState)
    ζmax = mat.ζ
    ζmin = 0.1*ζmax # minimum value to prevent excessive stiffness degradation

    wn = state.w[1]
    w0 = 0.3*mat.wc # characteristic relative displacement for stiffness degradation

    ζ  = clamp(ζmax - (ζmax-ζmin)*wn/w0, ζmin, ζmax) # linear degradation of stiffness with opening displacement
    
    kn = mat.E*ζ/state.h
    G  = mat.E/(2*(1 + mat.ν))
    ks = G*ζ/state.h

    return kn, ks
end


function calcD(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState)
    σmax   = calc_σmax(mat, state.up)
    kn, ks = calc_kn_ks(mat, state)
    tiny   = 1e-6*mat.ft

    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    if state.Δλ == 0.0  # Elastic 
        return De
    # elseif σmax <= tiny && state.w[1] >= 0.0
    #     Dep = De*1e-6
    #     return Dep
    else
        n, ∂f∂σmax = yield_derivs(mat, state.σ)
        m = potential_derivs(mat, state.σ)
        H = deriv_σmax_up(mat, state.up)  # ∂σmax/∂up
        Hcap = -mat.ft/(0.5*mat.wc)
        H    = max(H, Hcap) # cap degradation to prevent numerical issues
        
        De_m  = De*m
        nT_De = n'*De
        den   = dot(n, De_m) - ∂f∂σmax*H*norm(m)
        Dep   = De - (De_m*nT_De)/den

        return Dep
    end
end


function plastic_update(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState, cstate::MohrCoulombCohesiveState, σtr::Vec3)
    kn, ks = calc_kn_ks(mat, state)
    σntr, τ1tr, τ2tr = σtr

    μ = mat.μ
    ψ = mat.ψ

    τtr = √(τ1tr^2 + τ2tr^2 + eps())

    maxits    = 50
    converged = false
    Δλ        = 0.0
    up        = cstate.up
    σ         = cstate.σ
    σmax      = calc_σmax(mat, up)
    tol       = mat.ft*1e-8
    σtol      = mat.ft*1e-6
    σ0        = copy(σ)

    for i in 1:maxits
        den_σn = 1.0 + Δλ*kn*ψ^2
        den_τ  = 1.0 + Δλ*ks

        # stresses at current iterate
        σn = (σntr < 0) ? σntr : σntr/den_σn
        τ1 = τ1tr/den_τ
        τ2 = τ2tr/den_τ
        τ  = √(τ1^2 + τ2^2 + eps())
        σ  = Vec3(σn, τ1, τ2)

        # m at current iterate stress
        m      = potential_derivs(mat, σ)
        norm_m = norm(m)
        unit_m = m / (norm_m + eps())

        # softening variable at current iterate
        up   = cstate.up + Δλ*norm_m
        σmax = calc_σmax(mat, up)
        H    = deriv_σmax_up(mat, up)

        # residual
        f = τ + μ*(σn - σmax)
        if abs(f) < tol && maximum(abs, σ-σ0) < σtol
            converged = true
            break
        end

        σ0 = copy(σ)

        # derivatives
        if σntr<0
            ∂σn∂Δλ = 0.0
            ∂m∂Δλ  = Vec3( 0.0, -τ1tr*ks/den_τ^2, -τ2tr*ks/den_τ^2 )
        else
            ∂σn∂Δλ = -σntr*kn*ψ^2/den_σn^2
            ∂m∂Δλ  = Vec3( -σntr*kn*ψ^4/den_σn^2, -τ1tr*ks/den_τ^2, -τ2tr*ks/den_τ^2 )
        end

        ∂τ∂Δλ    = -τtr*ks/den_τ^2
        ∂up∂Δλ   = norm_m + Δλ*dot(unit_m, ∂m∂Δλ)
        ∂σmax∂Δλ = H*∂up∂Δλ

        ∂f∂Δλ = ∂τ∂Δλ + ∂σn∂Δλ*μ - ∂σmax∂Δλ*μ
        Δλ    = max(Δλ - f/∂f∂Δλ, 0.0)
    end

    if converged
        state.σ  = σ
        state.Δλ = Δλ
        state.up = up
        return success()
    else
        failure("MohrCoulombCohesive: plastic update failed.")
    end
end


function update_state(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState, cstate::MohrCoulombCohesiveState, Δw::Vector{Float64})
    kn, ks = calc_kn_ks(mat, state)
    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    σmax = calc_σmax(mat, cstate.up)

    # σ trial and f trial
    σtr = cstate.σ + De*Δw
    ftr = yield_func(mat, σtr, σmax)

    # Elastic and EP integration
    if ftr <= 0.0
        # Pure elastic increment
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # Plastic increment
        status = plastic_update(mat, state, cstate, σtr)
        failed(status) && return state.σ, status
    end

    state.w = cstate.w + Δw
    Δσ      = state.σ - cstate.σ
    return Δσ, success()
end


function state_values(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState)
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


function output_keys(::MohrCoulombCohesive)
    return Symbol[:w, :σn, :τ, :up]
end