 #This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export AsinhYieldCohesive


"""
    AsinhYieldCohesive(; E, nu=0.0, ft, fc, zeta=5.0, wc, GF, ft_law=:hordijk, alpha=1.5, gamma=0.1, theta=1.5)

Constitutive model for cohesive elements with a power-lay yield surface ans ft_law in tension.  
The tensile ft_law branch is regularized through a measure of the
bulk element size `h` to ensure mesh-objective fracture energy dissipation.

# Keyword arguments
- `E::Real`:  
  Young’s modulus from the bulk material (must be > 0).
- `nu::Real`:  
  Poisson’s ratio (0 ≤ ν < 0.5).
- `fc::Real`:  
  Compressive strength (< 0).
  - `ft::Real`:  
  Tensile strength (> 0).
- `wc::Real`:  
  Critical crack opening (must be > 0 if given). Can be specified alternatively to `GF`.
- `mu::Real`:  
  Friction coefficient (> 0).
- `GF::Real`:  
  Fracture energy (must be > 0 if given). Can be specified alternatively to `wc`.
- `ft_law::Union{Symbol,AbstractSpline} = :hordijk`:  
  Softening law for post-peak tensile response. Options are:
  `:linear`, `:bilinear`, `:hordijk`, `:soft` or a custom function.
- `alpha::Real = 0.6`:  
  Parameter to control the shape of the yield surface (α > 0.5).
- `gamma::Real = 0.1`:  
  Parameter to control the residual shear strength (γ ≥ 0).
- `theta::Real = 1.5`:  
  Parameter to control the rate of reduction of shear strength (θ ≥ 0).
- `zeta::Real = 5.0`:  
  Factor to control elastic relative displacements in cohesive formulations (≥ 0).

# Returns
A `AsinhYieldCohesive` object.

# Notes
- Either `wc` or `GF` must be provided. If only `GF` is given, `wc` is computed
  internally based on the chosen ft_law.
- The frictional contribution is governed by `mu`.
- Normal and shear stiffnesses (`kn`, `ks`) are computed from the mechanical properties of
  the bulk material and the characteristic length `h` of the adjacent bulk elements.
"""
mutable struct AsinhYieldCohesive<:Constitutive
    E ::Float64
    ν ::Float64
    fc::Float64
    ft::Float64
    wc::Float64
    ψ ::Float64
    ft_law::Symbol
    ft_fun::Union{AbstractSpline,Nothing}
    α::Float64
    γ::Float64
    θ::Float64
    βini::Float64
    ζ ::Float64

    function AsinhYieldCohesive(; 
        E::Real = NaN,
        nu::Real = 0.0,
        fc::Real = NaN,
        ft::Real = NaN,
        wc::Real = NaN,
        GF::Real = NaN,
        psi::Real  = 1.0,
        ft_law::Union{Symbol,AbstractSpline} = :hordijk,
        alpha::Real = 0.6,
        gamma::Real = 0.1,
        theta::Real = 1.5,
        zeta::Real = 5.0,
    )

        @check E>0 "AsinhYieldCohesive: Young's modulus E must be > 0. Got $(repr(E))."
        @check 0<=nu<0.5 "AsinhYieldCohesive: Poisson ratio nu must be in the range [0, 0.5). Got $(repr(nu))."
        @check fc<0 "AsinhYieldCohesive: Compressive strength fc must be < 0. Got $(repr(fc))."
        @check ft>0 "AsinhYieldCohesive: Tensile strength ft must be > 0. Got $(repr(ft))."
        @check psi>0 "AsinhYieldCohesive: Dilatancy coefficient psi must be non-negative. Got $(repr(psi))."
        @check zeta>=0 "AsinhYieldCohesive: Factor zeta must be non-negative. Got $(repr(zeta))."
        @check alpha > 0.5 "AsinhYieldCohesive: alpha must be greater than 0.5. Got $(repr(alpha))."
        @check gamma >= 0.0 "AsinhYieldCohesive: gamma must be non-negative. Got $(repr(gamma))."
        @check theta >= 0.0 "AsinhYieldCohesive: theta must be non-negative. Got $(repr(theta))."
        @check ft_law in (:linear, :bilinear, :hordijk, :soft) || ft_law isa AbstractSpline "AsinhYieldCohesive: Unknown ft_law model: $ft_law. Supported models are :linear, :bilinear, :hordijk, :soft or a custom AbstractSpline."

        wc, ft_law, ft_fun, status = setup_tensile_strength(ft,  GF, wc, ft_law)
        failed(status) && throw(ArgumentError("AsinhYieldCohesive: " * status.message))


        α = alpha

        ta = 0.05*pi
        tb = 0.5*pi

        f(t) = begin
            a = fc/2*(1-cos(t)) # negative value
            b = -fc/2*(sin(t))  # positive value
            χ = (ft-a)/ft
            βini = b/asinh(α*χ)
            fc/2 - a + α*βini*b/(ft*√(α^2*χ^2 + 1))
        end
        
        t, _ = findroot(f, ta, tb, tol=1e-4, method=:default)
        t>pi/2 && throw(SerendipException("Invalid value for βini was found. Check fc and ft values"))
        
        a    = fc/2*(1-cos(t)) # negative value
        b    = -fc/2*(sin(t))  # positive value
        χ    = (ft-a)/ft
        βini = b/asinh(α*χ)

        return new(E, nu, fc, ft, wc, psi, ft_law, ft_fun, alpha, gamma, theta, βini, zeta)
    end
end


mutable struct AsinhYieldCohesiveState<:ConstState
    ctx::Context
    σ  ::Vec3        # stress
    w  ::Vec3        # relative displacements
    up ::Float64     # effective plastic relative displacement
    Δλ ::Float64     # plastic multiplier
    h  ::Float64     # characteristic length from bulk elements
    function AsinhYieldCohesiveState(ctx::Context)
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
compat_state_type(::Type{AsinhYieldCohesive}, ::Type{MechCohesive}) = AsinhYieldCohesiveState

function calc_β(mat::AsinhYieldCohesive, σmax::Float64)
    βini = mat.βini
    βres = mat.γ*βini
    return βres + (βini-βres)*(σmax/mat.ft)^mat.θ
end


function yield_func(mat::AsinhYieldCohesive, σ::Vec3, σmax::Float64)
    σn, τ1, τ2 = σ

    β = calc_β(mat, σmax)
    χ = (σmax - σn)/mat.ft
    τ = sqrt(τ1^2 + τ2^2)

    return τ - β*asinh(mat.α*χ)
end


function stress_strength_ratio(mat::AsinhYieldCohesive, σ::AbstractVector)
    σn, τ1, τ2 = σ

    σmax = calc_σmax(mat, 0.0)
    β    = calc_β(mat, σmax)
    χ    = (σmax - σn)/mat.ft
    τmax = β*asinh(mat.α*χ)
    τ    = √(τ1^2 + τ2^2)
    return max(σn/σmax, τ/τmax)
end


function cap_stress(mat::AsinhYieldCohesive, σ::AbstractVector)
    σn, τ1, τ2 = σ
    σmax = calc_σmax(mat, 0.0)
    return Vec3(min(σn, σmax), τ1, τ2)
end


function calc_∂β∂σmax(mat::AsinhYieldCohesive, σmax::Float64)
    σmax == 0.0 && return 0.0
    βini = mat.βini
    βres = mat.γ*βini
    return (βini - βres)*mat.θ/mat.ft*(σmax/mat.ft)^(mat.θ-1)
end


function yield_derivs(mat::AsinhYieldCohesive, σ::Vec3, σmax::Float64)
    σn, τ1, τ2 = σ
    τ = √(τ1^2 + τ2^2 + eps())
    ft = mat.ft
    α  = mat.α
    β  = calc_β(mat, σmax)
    χ  = (σmax - σn)/ft
    
    ∂f∂σn   = α*β/(ft*√(α^2*χ^2 + 1))
    ∂f∂σ    = Vec3( ∂f∂σn, τ1/τ, τ2/τ )
    ∂β∂σmax = calc_∂β∂σmax(mat, σmax)
    ∂f∂σmax = -∂β∂σmax*asinh(α*χ) - α*β/(ft*√(α^2*χ^2 + 1))

    return ∂f∂σ, ∂f∂σmax
end


function potential_derivs(mat::AsinhYieldCohesive, σ::Vec3)
    σn, τ1, τ2 = σ

    if σn < 0.0 
        return Vec3( 0.0, τ1, τ2 )
    else
        ψ = mat.ψ
        return Vec3( ψ^2*σn + eps()^0.5, τ1, τ2 )
    end
end


function calc_σmax(mat::AsinhYieldCohesive, up::Float64)
    return calc_tensile_strength(mat, up)
end


function deriv_σmax_up(mat::AsinhYieldCohesive, up::Float64)
    # ∂σmax/∂up
    return calc_tensile_strength_derivative(mat, up)
end


function calc_kn_ks(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2*(1+mat.ν))
    ks = G*mat.ζ/state.h
    return kn, ks
end


function calcD(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState)
    σmax   = calc_σmax(mat, state.up)
    kn, ks = calc_kn_ks(mat, state)
    tiny   = 1e-6*mat.ft

    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax <= tiny && state.w[1] >= 0.0
        Dep = De*1e-8
        return Dep
    else
        n, ∂f∂σmax = yield_derivs(mat, state.σ, σmax)
        m = potential_derivs(mat, state.σ)
        H = deriv_σmax_up(mat, state.up)  # ∂σmax/∂up
        
        De_m  = De*m
        nT_De = n'*De
        den   = dot(n, De_m) - ∂f∂σmax*H*norm(m)
        Dep   = De - (De_m*nT_De)/den

        return Dep
    end
end


function plastic_update(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, cstate::AsinhYieldCohesiveState, σtr::Vec3)
    kn, ks = calc_kn_ks(mat, state)
    σntr, τ1tr, τ2tr = σtr

    ψ = mat.ψ

    τtr = √(τ1tr^2 + τ2tr^2 + eps())

    maxits    = 30
    converged = false
    Δλ        = 0.0
    up        = cstate.up
    σ         = cstate.σ
    σmax      = calc_σmax(mat, up)
    tol       = mat.ft*1e-6

    for i in 1:maxits
        den_σn = 1.0 + Δλ*kn*ψ^2
        den_τ  = 1.0 + Δλ*ks

        # stresses at current iterate
        σn = (σntr < 0) ? σntr : σntr/den_σn
        τ1 = τ1tr/den_τ
        τ2 = τ2tr/den_τ
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
        f = yield_func(mat, σ, σmax)
        if abs(f) < tol
            converged = true
            break
        end

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

        ∂β∂σmax = calc_∂β∂σmax(mat, σmax)
        ∂β∂Δλ   = ∂β∂σmax*∂σmax∂Δλ
        χ       = (σmax - σn)/mat.ft
        β       = calc_β(mat, σmax)
        ∂f∂Δλ   = ∂τ∂Δλ - ∂β∂Δλ*asinh(mat.α*χ) - mat.α*β/mat.ft/√(mat.α^2*χ^2+1)*(∂σmax∂Δλ - ∂σn∂Δλ)
        Δλ      = max(Δλ - f/∂f∂Δλ, 0.0)
    end

    if converged
        state.σ  = σ
        state.Δλ = Δλ
        state.up = up
        return success()
    else
        failure("AsinhYieldCohesive: plastic update failed.")
    end
end


function update_state(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, cstate::AsinhYieldCohesiveState, Δw::Vector{Float64})
    kn, ks = calc_kn_ks(mat, state)
    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    σmax   = calc_σmax(mat, cstate.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("AsinhYieldCohesive: Invalid value for relative displacement: Δw = $Δw")
    end

    # σ trial and f trial
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
    return Symbol[:w, :σn, :τ, :up]
end