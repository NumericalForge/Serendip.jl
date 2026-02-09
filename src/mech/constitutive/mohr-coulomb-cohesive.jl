 #This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MohrCoulombCohesive


"""
    MohrCoulombCohesive(; E, nu=0.0, ft, GF, wc, mu, ft_law=:hordijk, zeta=5.0)

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
        ft_law::Union{Symbol,AbstractSpline} = :hordijk,

        zeta::Real=5.0
    )
        @check E>0 "MohrCoulombCohesive: Young's modulus E must be > 0. Got $(repr(E))."
        @check 0<=nu<0.5 "MohrCoulombCohesive: Poisson ratio nu must be in the range [0, 0.5). Got $(repr(nu))."
        @check ft>0 "MohrCoulombCohesive: Tensile strength ft must be > 0. Got $(repr(ft))."
        @check mu>0 "MohrCoulombCohesive: Friction coefficient mu must be non-negative. Got $(repr(mu))."
        @check zeta>=0 "MohrCoulombCohesive: Factor zeta must be non-negative. Got $(repr(zeta))."

        wc, ft_law, ft_fun, status = setup_tensile_strength(ft, GF, wc, ft_law)
        failed(status) && throw(ArgumentError("MohrCoulombCohesive: " * status.message))

        return new(E, nu, ft, wc, mu, ft_law, ft_fun, zeta)
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
    return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
end


function stress_strength_ratio(mat::MohrCoulombCohesive, σ::AbstractVector)
    σmax = calc_σmax(mat, 0.0)
    τmax = (σmax - σ[1])*mat.μ
    τ    = sqrt(σ[2]^2 + σ[3]^2)
    return max(σ[1]/σmax, τ/τmax)
end


function yield_deriv(mat::MohrCoulombCohesive, σ::Vec3)
    τ = sqrt(σ[2]^2 + σ[3]^2)
    return Vec3( mat.μ, σ[2]/τ, σ[3]/τ )
end


function potential_derivs(mat::MohrCoulombCohesive, σ::Vec3)
    if σ[1] >= 0.0 
        # G1:
        r = Vec3( 2*σ[1]*mat.μ^2, 2*σ[2], 2*σ[3] )
    else
        # G2:
        r = Vec3( 0, 2*σ[2], 2*σ[3] )
    end

    return r
end


function calc_σmax(mat::MohrCoulombCohesive, up::Float64)
    return calc_tensile_strength(mat, up)
end


function deriv_σmax_upa(mat::MohrCoulombCohesive, up::Float64)
    # ∂σmax/∂up
    return calc_tensile_strength_derivative(mat, up)
end


function calc_kn_ks(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2*(1 + mat.ν))
    ks = G*mat.ζ/state.h
    return kn, ks
end


function calcD(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState)
    σmax   = calc_σmax(mat, state.up)
    kn, ks = calc_kn_ks(mat, state)

    De = @SMatrix [ kn   0.0  0.0
                    0.0  ks   0.0
                    0.0  0.0  ks ]

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 && state.w[1] >= 0.0
        # Dep  = De*1e-4
        Dep  = De.*1e-3
        return Dep
    else
        v = yield_deriv(mat, state.σ)
        r = potential_derivs(mat, state.σ)
        y = -mat.μ # ∂F/∂σmax
        m = deriv_σmax_upa(mat, state.up)  # ∂σmax/∂up

        #Dep  = De - De*r*v'*De/(v'*De*r - y*m*norm(r))
        den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - y*m*norm(r)

        Dep = @SMatrix [  kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                         -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                         -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]

        return Dep
    end
end


function nonlinear_update(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState, cstate::MohrCoulombCohesiveState, σtr::Vec3)
    # Compute Δλ

    maxits = 50
    Δλ     = 0.0
    up     = 0.0
    tol    = 1e-4
    μ      = mat.μ
    kn, ks = calc_kn_ks(mat, state)
    
    for i in 1:maxits
    
        # quantities at n+1
        if σtr[1]>0
            σ     = Vec3( σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) )
            dσdΔλ = Vec3( -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 )
            drdΔλ = Vec3( -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 )
        else
            σ     = Vec3( σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) )
            dσdΔλ = Vec3( 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 )
            drdΔλ = Vec3( 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 )
        end
                 
        r        = potential_derivs(mat, σ)
        norm_r   = norm(r)
        up       = cstate.up + Δλ*norm_r
        σmax     = calc_σmax(mat, up)
        m        = deriv_σmax_upa(mat, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
    
        f = sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*μ
        if (σ[2]==0 && σ[3]==0) 
            dfdΔλ = (dσdΔλ[1] - dσmaxdΔλ)*μ
        else
            dfdΔλ = 1/sqrt(σ[2]^2 + σ[3]^2) * (σ[2]*dσdΔλ[2] + σ[3]*dσdΔλ[3]) + (dσdΔλ[1] - dσmaxdΔλ)*μ
        end
    
        Δλ = Δλ - f/dfdΔλ
    
        abs(f) < tol && break
    
        if i == maxits || isnan(Δλ)
            # warn("""MohrCoulombCohesive: Could not find Δλ. This may happen when the system
            # becomes hypostatic and thus the global stiffness matrix is nearly singular.
            # Increasing the mesh refinement may result in a nonsingular matrix.
            # """)
            # warn("iterations=$i Δλ=$Δλ")
            return failure("MohrCoulombCohesive: Could nof find Δλ.")
        end
    end

    # Update σ and upa
    if σtr[1] > 0
        σ = Vec3( σtr[1]/(1 + 2*Δλ*kn*(μ^2)), σtr[2]/(1 + 2*Δλ*ks), σtr[3]/(1 + 2*Δλ*ks) )
    else
        σ = Vec3( σtr[1], σtr[2]/(1 + 2*Δλ*ks), σtr[3]/(1 + 2*Δλ*ks) )
    end    
    
    state.Δλ = Δλ
    state.σ  = σ
    r        = potential_derivs(mat, σ)
    state.up = cstate.up + state.Δλ*norm(r)

    return success()
end


function update_state(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState, cstate::MohrCoulombCohesiveState, Δw::Vector{Float64})
    kn, ks = calc_kn_ks(mat, state)
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
    if σmax == 0.0 && cstate.w[1] >= 0.0
        # Return to apex:
        r1 = Vec3( σtr[1]/kn, σtr[2]/ks, σtr[3]/ks )
        r  = r1/norm(r1)
        state.Δλ = norm(r1)
        
        state.up = cstate.up + state.Δλ
        state.σ  = σtr - state.Δλ*De*r
    elseif ftr <= 0.0
        # Pure elastic increment
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # Plastic increment
        status = nonlinear_update(mat, state, cstate, σtr)
        failed(status) && return state.σ, status
    end

    state.w = cstate.w + Δw
    Δσ      = state.σ - cstate.σ
    return Δσ, success()
end


function state_values(mat::MohrCoulombCohesive, state::MohrCoulombCohesiveState)
    σmax = calc_σmax(mat, state.up)
    τ = sqrt(state.σ[2]^2 + state.σ[3]^2)

    return Dict(
        :w    => state.w[1],
        :σn   => state.σ[1],
        :τ    => τ,
        :up   => state.up,
        :σmax => σmax
      )
end


function output_keys(::MohrCoulombCohesive)
    return Symbol[:w, :σn, :τ, :up]
end