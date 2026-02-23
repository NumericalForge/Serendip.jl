 #This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MohrCoulombContact


"""
    MohrCoulombContact(; ft, wc, GF, mu, kn, ks, ft_law=:hordijk)

Constitutive model for interface/contact elements with a Mohr–Coulomb strength
criterion. It combines normal and shear stiffness, tensile strength, friction,
and a post-peak tensile softening law defined either by the critical crack
opening `wc` or by the fracture energy `GF`.

# Keyword arguments
- `ft::Real`
  Tensile strength (≥ 0).
- `mu::Real`
  Friction coefficient (> 0).
- `kn::Real`
  Normal stiffness per unit area (> 0).
- `ks::Real`
  Shear stiffness per unit area (> 0).
- `wc::Real`
  Critical crack opening (> 0 if provided). May be computed from `GF`.
- `GF::Real`
  Mode-I fracture energy (> 0 if provided). May be used to compute `wc`.
- `ft_law::Union{Symbol,AbstractSpline} = :hordijk`
  Tensile softening law. Use a symbol `:linear`, `:bilinear`, or `:hordijk`,
  or pass a Spline.

# Returns
An `MohrCoulombContact` object.

# Notes
- Provide either `wc` or `GF`. If only `GF` is given, `wc` is computed based on `ft_law`.
- `kn` and `ks` control the elastic response before reaching the strength envelope.
"""
mutable struct MohrCoulombContact<:Constitutive
    ft ::Float64
    wc ::Float64
    μ  ::Float64
    ft_law::Symbol
    ft_fun::Union{AbstractSpline,Nothing}
    kn ::Float64
    ks ::Float64

    function MohrCoulombContact(; 
            ft::Real=NaN,
            wc::Real=NaN,
            GF::Real=NaN,
            mu::Real=NaN,
            kn::Real=NaN,
            ks::Real=NaN,
            ft_law::Union{Symbol,AbstractSpline} = :hordijk,
        )

        @check ft>=0 "MohrCoulombContact: Tensile strength ft must be >= 0. Got $(repr(ft))."
        @check mu>0 "MohrCoulombContact: Friction coefficient mu must be non-negative. Got $(repr(mu))."
        @check kn>0 "MohrCoulombContact: Normal stiffness per area kn must be non-negative. Got $(repr(kn))."
        @check ks>0 "MohrCoulombContact: Shear stiffness per area ks must be non-negative. Got $(repr(ks))."

        wc, ft_law, ft_fun, status = setup_tensile_strength(ft, GF, wc, ft_law)
        failed(status) && throw(ArgumentError("MohrCoulombContact: " * status.message))

        this = new(ft, wc, mu, ft_law, ft_fun, kn, ks)
        return this
    end
end


mutable struct MohrCoulombContactState<:ConstState
    ctx::Context
    σ  ::Vec3        # stress
    w  ::Vec3        # relative displacements
    up ::Float64     # effective plastic relative displacement
    Δλ ::Float64     # plastic multiplier
    function MohrCoulombContactState(ctx::Context)
        this    = new(ctx)
        this.σ  = zeros(Vec3)
        this.w  = zeros(Vec3)
        this.up = 0.0
        this.Δλ = 0.0
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{MohrCoulombContact}, ::Type{MechContact}) = MohrCoulombContactState


function yield_func(mat::MohrCoulombContact, σ::Vec3, σmax::Float64)
        σn, τ1, τ2 = σ
    return √(τ1^2 + τ2^2) + (σn - σmax)*mat.μ
end


function yield_derivs(mat::MohrCoulombContact, σ::Vec3)
    σn, τ1, τ2 = σ

    τ = √(τ1^2 + τ2^2 + eps())

    ∂f∂σmax = -mat.μ
    return Vec3( mat.μ, τ1/τ, τ2/τ ), ∂f∂σmax
end


function potential_derivs(mat::MohrCoulombContact, σ::Vec3)
    σn, τ1, τ2 = σ

    if σn < 0.0 
        return Vec3( 0.0, τ1, τ2 )
    else
        μ = mat.μ
        return Vec3( μ^2*σn + eps()^0.5, τ1, τ2 )
    end
end


function calc_σmax(mat::MohrCoulombContact, up::Float64)
    return calc_tensile_strength(mat, up)
end


function deriv_σmax_up(mat::MohrCoulombContact, up::Float64)
    # ∂σmax/∂up
    return calc_tensile_strength_derivative(mat, up)
end


function calcD(mat::MohrCoulombContact, state::MohrCoulombContactState)
    σmax   = calc_σmax(mat, state.up)
    kn, ks = mat.kn, mat.ks
    tiny   = 1e-6*mat.ft

    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax <= tiny && state.w[1] >= 0.0
        Dep = De*1e-3
        return Dep
    else
        n, ∂f∂σmax = yield_derivs(mat, state.σ)
        m = potential_derivs(mat, state.σ)
        H = deriv_σmax_up(mat, state.up)  # ∂σmax/∂up
        
        De_m  = De*m
        nT_De = n'*De
        den   = dot(n, De_m) - ∂f∂σmax*H*norm(m)
        Dep   = De - (De_m*nT_De)/den

        return Dep
    end
end


function plastic_update(mat::MohrCoulombContact, state::MohrCoulombContactState, cstate::MohrCoulombContactState, σtr::Vec3)
    kn, ks = mat.kn, mat.ks
    σntr, τ1tr, τ2tr = σtr

    μ = mat.μ

    τtr = √(τ1tr^2 + τ2tr^2 + eps())

    maxits    = 30
    converged = false
    Δλ        = 0.0
    up        = cstate.up
    σ         = cstate.σ
    σmax      = calc_σmax(mat, up)
    tol       = mat.ft*1e-6

    for i in 1:maxits
        den_σn = 1.0 + Δλ*kn*μ^2
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
        if abs(f) < tol
            converged = true
            break
        end

        # derivatives
        if σntr<0
            ∂σn∂Δλ = 0.0
            ∂m∂Δλ  = Vec3( 0.0, -τ1tr*ks/den_τ^2, -τ2tr*ks/den_τ^2 )
        else
            ∂σn∂Δλ = -σntr*kn*μ^2/den_σn^2
            ∂m∂Δλ  = Vec3( -σntr*kn*μ^4/den_σn^2, -τ1tr*ks/den_τ^2, -τ2tr*ks/den_τ^2 )
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


function update_state(mat::MohrCoulombContact, state::MohrCoulombContactState, cstate::MohrCoulombContactState, Δw::Vector{Float64})
    ks, kn = mat.ks, mat.kn
    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    σmax   = calc_σmax(mat, cstate.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("MohrCoulombCohesive: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr = cstate.σ + De*Δw
    ftr = yield_func(mat, σtr, σmax)

    # Elastic and EP integration
    if σmax == 0.0 && cstate.w[1] + Δw[1] >= 0.0
        # traction-free after full decohesion
        state.σ   = Vec3(0.0, 0.0, 0.0)
        state.Δλ  = 1.0
        state.up  = max(cstate.up, norm(cstate.w + Δw))
    elseif ftr <= 0.0
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


function state_values(mat::MohrCoulombContact, state::MohrCoulombContactState)
    σmax = calc_σmax(mat, state.up)
    σn, τ1, τ2 = state.σ
    wn, s1, s2 = state.w
    τ = sqrt(τ1^2 + τ2^2)
    s = sqrt(s1^2 + s2^2)
    
    return Dict(
        :w    => wn,
        :s    => s,
        :σn   => σn,
        :τ    => τ,
        :up   => state.up,
        :σmax => σmax
      )
end


function output_keys(::MohrCoulombContact)
    return Symbol[:w, :s, :σn, :τ, :up]
end