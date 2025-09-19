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
    ft_law::Symbol
    ft_fun::Union{AbstractSpline,Nothing}
    α::Float64
    γ::Float64
    θ::Float64
    β0::Float64
    ζ ::Float64

    function AsinhYieldCohesive(; 
        E::Real = NaN,
        nu::Real = 0.0,
        fc::Real = NaN,
        ft::Real = NaN,
        wc::Real = NaN,
        GF::Real = NaN,
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
            β0 = b/asinh(α*χ)
            fc/2 - a + α*β0*b/(ft*√(α^2*χ^2 + 1))
        end
        
        t, _ = findroot(f, ta, tb, tol=1e-4, method=:default)
        t>pi/2 && throw(SerendipException("Invalid value for β0 was found. Check fc and ft values"))
        
        a  = fc/2*(1-cos(t)) # negative value
        b  = -fc/2*(sin(t))  # positive value
        χ  = (ft-a)/ft
        β0 = b/asinh(α*χ)

        return new(E, nu, fc, ft, wc, ft_law, ft_fun, alpha, gamma, theta, β0, zeta)
    end
end


mutable struct AsinhYieldCohesiveState<:IpState
    ctx::Context
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up ::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    w_rate::Float64       # w rate wrt T
    function AsinhYieldCohesiveState(ctx::Context)
        this = new(ctx)
        ndim = ctx.ndim
        this.σ   = zeros(ndim)
        this.w   = zeros(ndim)
        this.up  = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{AsinhYieldCohesive}, ::Type{MechInterface}, ctx::Context) = AsinhYieldCohesiveState


function paramsdict(mat::AsinhYieldCohesive)
    mat = OrderedDict( string(field) => getfield(mat, field) for field in fieldnames(typeof(mat)) )

    if mat.ft_law in (:hordijk, :soft)
        mat["GF"] = 0.1943*mat.ft*mat.wc
    elseif mat.ft_law == :bilinear
        mat["GF"] = mat.ft*mat.wc/5
    elseif mat.ft_law == :linear
        mat["GF"] = mat.ft*mat.wc/2
    else
        mat["GF"] = NaN
    end
    return mat
end


function yield_func(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, σ::Array{Float64,1}, σmax::Float64)
    β = calc_β(mat, σmax)
    χ = (σmax-σ[1])/mat.ft

    if state.ctx.ndim == 3
        τnorm = sqrt(σ[2]^2 + σ[3]^2)
    else
        τnorm = abs(σ[2])
    end

    return τnorm - β*asinh(mat.α*χ)
end


function yield_derivs(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, σ::Array{Float64,1}, σmax::Float64)
    ft   = mat.ft
    α    = mat.α
    β    = calc_β(mat, σmax)
    βres = mat.γ*mat.β0
    χ    = (σmax-σ[1])/ft
    
    dfdσn  = α*β/(ft*√(α^2*χ^2 + 1))
    
    if state.ctx.ndim == 3
        τnorm = sqrt(σ[2]^2 + σ[3]^2)
        dfdσ = [ dfdσn, σ[2]/τnorm, σ[3]/τnorm]
    else
        dfdσ = [ dfdσn, sign(σ[2]) ]
    end

    if σmax>0
        θ = mat.θ
        dβdσmax = (β-βres)*θ/ft*(σmax/ft)^(θ-1)
    else 
        dβdσmax = 0.0
    end
    dfdσmax = -dβdσmax*asinh(α*χ) - α*β/(ft*√(α^2*χ^2 + 1))

    return dfdσ, dfdσmax
end


function potential_derivs(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, σ::Array{Float64,1})
    ndim = state.ctx.ndim
    if ndim == 3
        if σ[1] > 0.0 
            # G1:
            r = Float64[ 2.0*σ[1], 2.0*σ[2], 2.0*σ[3]]
        else
            # G2:
            r = Float64[ 0.0, 2.0*σ[2], 2.0*σ[3] ]
        end
        if r[1]==r[2]==r[3]==0.0
            r = Float64[ 1.0, 0.0, 0.0 ] # important
        end
    else
        if σ[1] > 0.0 
            # G1:
            r = Float64[ 2*σ[1], 2*σ[2]]
        else
            # G2:
            r = Float64[ 0.0, 2*σ[2] ]
        end
        if r[1]==r[2]==0.0
            r = Float64[ 1.0, 0.0 ] # important
        end
    end
    return r
end


function calc_β(mat::AsinhYieldCohesive, σmax::Float64)
    βres  = mat.γ*mat.β0
    return βres + (mat.β0-βres)*(σmax/mat.ft)^mat.θ
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
    G  = mat.E/(2.0*(1.0+mat.ν))
    ks = G*mat.ζ/state.h
    return kn, ks
end


function calcD(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState)

    ndim = state.ctx.ndim
    kn, ks = calc_kn_ks(mat, state)
    σmax = calc_σmax(mat, state.up)

    De = diagm([kn, ks, ks][1:ndim])

    if state.Δλ == 0.0  # Elastic 
        return De
    elseif σmax == 0.0 && state.w[1] >= 0.0
        Dep = De*1e-4
        # Dep = De*1e-3
        return Dep
    else
        r = potential_derivs(mat, state, state.σ) # ∂g/∂σ
        dfdσ, dfdσmax = yield_derivs(mat, state, state.σ, σmax)
        dσmaxdup = deriv_σmax_up(mat, state.up)  # ∂σmax/∂up

        if ndim == 3
            den = kn*r[1]*dfdσ[1] + ks*r[2]*dfdσ[2] + ks*r[3]*dfdσ[3] - dfdσmax*dσmaxdup*norm(r)

            Dep = [   kn - kn^2*r[1]*dfdσ[1]/den    -kn*ks*r[1]*dfdσ[2]/den      -kn*ks*r[1]*dfdσ[3]/den
                     -kn*ks*r[2]*dfdσ[1]/den         ks - ks^2*r[2]*dfdσ[2]/den  -ks^2*r[2]*dfdσ[3]/den
                     -kn*ks*r[3]*dfdσ[1]/den        -ks^2*r[3]*dfdσ[2]/den        ks - ks^2*r[3]*dfdσ[3]/den ]
        else
            den = kn*r[1]*dfdσ[1] + ks*r[2]*dfdσ[2] - dfdσmax*dσmaxdup*norm(r)

            Dep = [   kn - kn^2*r[1]*dfdσ[1]/den    -kn*ks*r[1]*dfdσ[2]/den      
                     -kn*ks*r[2]*dfdσ[1]/den         ks - ks^2*r[2]*dfdσ[2]/den  ]
        end
        return Dep
    end
end


function calc_σ_up_Δλ(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, σtr::Array{Float64,1})
    ndim = state.ctx.ndim
    Δλ   = 0.0
    up   = 0.0
    σ    = zeros(ndim)
    σ0   = zeros(ndim)

    tol    = 1e-6
    maxits = 50
    for i in 1:maxits
        kn, ks = calc_kn_ks(mat, state)

        # quantities at n+1
        if ndim == 3
            if σtr[1]>0
                 σ     = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ -2*kn*σtr[1]/(1+2*Δλ*kn)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2,  -2*ks*σtr[3]/(1+2*Δλ*ks)^2 ]
            end
        else
            if σtr[1]>0
                 σ     = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ -2*kn*σtr[1]/(1+2*Δλ*kn)^2,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
            else
                 σ     = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
                 dσdΔλ = [ 0,  -2*ks*σtr[2]/(1+2*Δλ*ks)^2 ]
             end
        end

        drdΔλ = 2*dσdΔλ
                 
        r      = potential_derivs(mat, state, σ)
        norm_r = norm(r)
        up     = state.up + Δλ*norm_r
        σmax   = calc_σmax(mat, up)
        f      = yield_func(mat, state, σ, σmax)
        dfdσ, dfdσmax = yield_derivs(mat, state, σ, σmax)

        dσmaxdup = deriv_σmax_up(mat, up)
        dσmaxdΔλ = dσmaxdup*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        
        if Δλ<=0 || isnan(Δλ) || i==maxits
            # return 0.0, state.σ, 0.0, failure("AsinhYieldCohesive: failed to find Δλ")
            # switch to bissection method
            # return calc_σ_up_Δλ_bis(mat, state, σtr)
            σ, up, Δλ, status = calc_σ_up_Δλ_bis(mat, state, σtr)
            failed(status) && return state.σ, 0.0, 0.0, failure("AsinhYieldCohesive: failed to find Δλ")
        end

        if maximum(abs, σ-σ0) <= tol
            break
        end
        σ0 .= σ
    end

    return σ, up, Δλ, success()
end


function calc_σ_up(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, σtr::Array{Float64,1}, Δλ::Float64)
    ndim = state.ctx.ndim
    kn, ks  = calc_kn_ks(mat, state)

    if ndim == 3
        if σtr[1]>0
            σ = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        else
            σ = [ σtr[1],  σtr[2]/(1+2*Δλ*ks),  σtr[3]/(1+2*Δλ*ks) ]
        end
    else
        if σtr[1]>0
            σ = [ σtr[1]/(1+2*Δλ*kn),  σtr[2]/(1+2*Δλ*ks) ]
        else
            σ = [ σtr[1],  σtr[2]/(1+2*Δλ*ks) ]
        end
    end

    r  = potential_derivs(mat, state, σ)
    up = state.up + Δλ*norm(r)
    return σ, up
end

function calc_σ_up_Δλ_bis(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, σtr::Array{Float64,1})
    ndim    = state.ctx.ndim
    kn, ks  = calc_kn_ks(mat, state)
    De      = diagm([kn, ks, ks][1:ndim])
    r       = potential_derivs(mat, state, state.σ)

    ff(Δλ)  = begin
        # quantities at n+1
        σ, up = calc_σ_up(mat, state, σtr, Δλ)
        σmax = deriv_σmax_upa(mat, up)
        yield_func(mat, state, σ, σmax)
    end

    # find root interval from Δλ estimative
    Δλ0 = norm(σtr-state.σ)/norm(De*r)
    a, b, status = findrootinterval(ff, 0.0, Δλ0)
    failed(status) && return state.σ, 0.0, 0.0, status

    Δλ, status = findroot(ff, a, b, ftol=1e-5, method=:bisection)
    failed(status) && return state.σ, 0.0, 0.0, status

    σ, up = calc_σ_up(mat, state, σtr, Δλ)
    return σ, up, Δλ, success()  

end


function update_state(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState, Δw::Array{Float64,1})

    ndim = state.ctx.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(mat, state)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, state.up)


    if isnan(Δw[1]) || isnan(Δw[2])
        alert("AsinhYieldCohesive: Invalid value for joint displacement: Δw = $Δw")
    end
    

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr, σmax)

    # Elastic and EP integration
    if Ftr <= 0.0
        state.Δλ  = 0.0
        state.σ  .= σtr
    elseif state.up>=mat.wc && σtr[1]>0
        if ndim==3
            Δup = norm([ σtr[1]/kn, σtr[2]/ks, σtr[3]/ks ])
        else
            Δup = norm([ σtr[1]/kn, σtr[2]/ks ])
        end
        state.up += Δup
        state.σ  .= 0.0
        state.Δλ  = 1.0
    else
        # Plastic increment
        σ, up, Δλ, status = calc_σ_up_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.up, state.Δλ = σ, up, Δλ
    end
    state.w += Δw
    Δσ = state.σ - σini
    return Δσ, success()
end


function state_values(mat::AsinhYieldCohesive, state::AsinhYieldCohesiveState)
    ndim = state.ctx.ndim
    σmax = calc_σmax(mat, state.up)
    τ = norm(state.σ[2:ndim])
    if ndim == 3
        return Dict(
            :w => state.w[1],
            :σn => state.σ[1],
            :τ  => τ,
            :s2 => state.σ[2],
            :s3 => state.σ[3],
            :up => state.up,
            :σmax => σmax
          )
    else
        return Dict(
            :w => state.w[1],
            :σn => state.σ[1],
            :τ  => τ,
            :s2 => state.σ[2],
            :up => state.up,
            :σmax => σmax
        )
    end
end


function output_keys(mat::AsinhYieldCohesive)
    return Symbol[:w, :σn, :τ, :up]
end