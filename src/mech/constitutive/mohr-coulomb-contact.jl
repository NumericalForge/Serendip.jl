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


mutable struct MohrCoulombContactState<:IpState
    ctx::Context
    σ  ::Vector{Float64} # stress
    w  ::Vector{Float64} # relative displacements
    up::Float64           # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    function MohrCoulombContactState(ctx::Context)
        this    = new(ctx)
        ndim    = ctx.ndim
        this.σ  = zeros(ndim)
        this.w  = zeros(ndim)
        this.up = 0.0
        this.Δλ = 0.0
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{MohrCoulombContact}, ::Type{MechContact}, ctx::Context) = MohrCoulombContactState


function yield_func(mat::MohrCoulombContact, state::MohrCoulombContactState, σ::Vector{Float64})
    ndim = state.ctx.ndim
    σmax = calc_σmax(mat, state.up)
    if ndim == 3
        return sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*mat.μ
    else
        return abs(σ[2]) + (σ[1]-σmax)*mat.μ
    end
end


function yield_deriv(mat::MohrCoulombContact, state::MohrCoulombContactState)
    ndim = state.ctx.ndim
    if ndim == 3
        return [ mat.μ, state.σ[2]/sqrt(state.σ[2]^2 + state.σ[3]^2), state.σ[3]/sqrt(state.σ[2]^2 + state.σ[3]^2)]
    else
        return [ mat.μ, sign(state.σ[2]) ]
    end
end


function potential_derivs(mat::MohrCoulombContact, state::MohrCoulombContactState, σ::Vector{Float64})
    ndim = state.ctx.ndim
    if ndim == 3
        if σ[1] >= 0.0 
            # G1:
            r = [ 2.0*σ[1]*mat.μ^2, 2.0*σ[2], 2.0*σ[3]]
        else
            # G2:
            r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
        end
    else
        if σ[1] >= 0.0 
            # G1:
            r = [ 2*σ[1]*mat.μ^2, 2*σ[2]]
        else
            # G2:
            r = [ 0.0, 2*σ[2] ]
        end
    end
    return r
end


function calc_σmax(mat::MohrCoulombContact, up::Float64)
    return calc_tensile_strength(mat, up)
end


function deriv_σmax_upa(mat::MohrCoulombContact, up::Float64)
    # ∂σmax/∂up
    return calc_tensile_strength_derivative(mat, up)
end


function calc_De(mat::MohrCoulombContact, state::MohrCoulombContactState)
    ndim = state.ctx.ndim
    ks, kn = mat.ks, mat.kn

    if ndim == 3
        De = [  kn  0.0  0.0
               0.0   ks  0.0
               0.0  0.0   ks ]
    else
        De = [  kn   0.0
                0.0  ks  ]
    end

    return De
end


function calc_Δλ(mat::MohrCoulombContact, state::MohrCoulombContactState, σtr::Vector{Float64})
    ndim   = state.ctx.ndim
    maxits = 100
    Δλ     = 0.0
    f      = 0.0
    up     = 0.0
    tol    = 1e-4
    μ      = mat.μ
    ks, kn = mat.ks, mat.kn

    for i in 1:maxits

        # quantities at n+1
        if ndim == 3
            if σtr[1]>0
                 σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
                 drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
                 drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2,  -4*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            end
        else
            if σtr[1]>0
                 σ     = [ σtr[1]/(1+2*Δλ*kn*μ^2),  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ -2*kn*μ^2*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
                 drdΔλ = [ -4*kn*μ^4*σtr[1]/(1+2*Δλ*kn*μ^2)^2,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
                 drdΔλ = [ 0,  -4*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
             end
        end
                 
        r        = potential_derivs(mat, state, σ)
        norm_r   = norm(r)
        up       = state.up + Δλ*norm_r
        σmax     = calc_σmax(mat, up)
        m        = deriv_σmax_upa(mat, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))

        if ndim == 3
            f = sqrt(σ[2]^2 + σ[3]^2) + (σ[1]-σmax)*μ
            if (σ[2]==0 && σ[3]==0) 
                dfdΔλ = (dσdΔλ[1] - dσmaxdΔλ)*μ
            else
                dfdΔλ = 1/sqrt(σ[2]^2 + σ[3]^2) * (σ[2]*dσdΔλ[2] + σ[3]*dσdΔλ[3]) + (dσdΔλ[1] - dσmaxdΔλ)*μ
            end
        else
            f = abs(σ[2]) + (σ[1]-σmax)*mat.μ
            dfdΔλ = sign(σ[2])*dσdΔλ[2] + (dσdΔλ[1] - dσmaxdΔλ)*μ
        end

        Δλ = Δλ - f/dfdΔλ

        abs(f) < tol && break

        if i == maxits || isnan(Δλ)
            # warn("""MohrCoulombContact: Could not find Δλ. This may happen when the system
            # becomes hypostatic and thus the global stiffness matrix is nearly singular.
            # Increasing the mesh refinement may result in a nonsingular matrix.
            # """)
            # warn("iterations=$i Δλ=$Δλ")
            return 0.0, failure("MohrCoulombContact: Could nof find Δλ.")
        end
    end
    return Δλ, success()
end


function calc_σ_upa(mat::MohrCoulombContact, state::MohrCoulombContactState, σtr::Vector{Float64})
    ndim = state.ctx.ndim
    μ = mat.μ
    ks, kn = mat.ks, mat.kn

    if ndim == 3
        if σtr[1] > 0
            σ = [σtr[1]/(1 + 2*state.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*state.Δλ*ks), σtr[3]/(1 + 2*state.Δλ*ks)]
        else
            σ = [σtr[1], σtr[2]/(1 + 2*state.Δλ*ks), σtr[3]/(1 + 2*state.Δλ*ks)]
        end    
    else
        if σtr[1] > 0
            σ = [σtr[1]/(1 + 2*state.Δλ*kn*(μ^2)), σtr[2]/(1 + 2*state.Δλ*ks)]
        else
            σ = [σtr[1], σtr[2]/(1 + 2*state.Δλ*ks)]
        end    
    end
    state.σ = σ
    r = potential_derivs(mat, state, state.σ)
    state.up += state.Δλ*norm(r)
    return state.σ, state.up
end


function calcD(mat::MohrCoulombContact, state::MohrCoulombContactState)
    ndim = state.ctx.ndim
    ks, kn = mat.ks, mat.kn
    De = calc_De(mat, state)
    σmax = calc_σmax(mat, state.up)

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 
        # Dep  = De*1e-10 
        # Dep  = De*1e-5
        # Dep  = De*1e-4
        Dep  = De*1e-3
        return Dep
    else
        v = yield_deriv(mat, state)
        r = potential_derivs(mat, state, state.σ)
        y = -mat.μ # ∂F/∂σmax
        m = deriv_σmax_upa(mat, state.up)  # ∂σmax/∂up

        #Dep  = De - De*r*v'*De/(v'*De*r - y*m*norm(r))

        if ndim == 3
            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - y*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                     -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]
        else
            den = kn*r[1]*v[1] + ks*r[2]*v[2] - y*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  ]
        end

        return Dep
    end
end


function update_state(mat::MohrCoulombContact, state::MohrCoulombContactState, Δw::Vector{Float64})
    ndim = state.ctx.ndim
    σini = copy(state.σ)

    De = calc_De(mat, state)
    σmax = calc_σmax(mat, state.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("MohrCoulombContact: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr = state.σ + De*Δw
    Ftr = yield_func(mat, state, σtr)

    # Elastic and EP integration
    if σmax == 0.0 && state.w[1] >= 0.0
        # Return to apex:
        if ndim==3
            r1 = [ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ]
            r = r1/norm(r1)
            state.Δλ = norm(r1)
        else
            r1 = [ σtr[1]/kn, σtr[2]/ks ]
            r = r1/norm(r1)
            state.Δλ = norm(r1)  
        end

        state.up += state.Δλ
        state.σ = σtr - state.Δλ*De*r     

    elseif Ftr <= 0.0
        # Pure elastic increment
        state.Δλ = 0.0
        state.σ  = copy(σtr) 

    else
        # Plastic increment
        state.Δλ, status = calc_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.up = calc_σ_upa(mat, state, σtr)
                      
        # Return to surface:
        F  = yield_func(mat, state, state.σ)   
        F > 1e-3 && alert("MohrCoulombContact: Yield function value ($F) outside tolerance")

    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function state_values(mat::MohrCoulombContact, state::MohrCoulombContactState)
    ndim = state.ctx.ndim
    σmax = calc_σmax(mat, state.up)
    τ = norm(state.σ[2:ndim])
    s = norm(state.w[2:ndim])

    return Dict(
       :w => state.w[1],
       :s  => s,
       :σn => state.σ[1],
       :τ  => τ,
       :up => state.up,
       :σmax => σmax
       )
end


function output_keys(::MohrCoulombContact)
    return Symbol[:w, :s, :σn, :τ, :up]
end