 #This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export PowerYieldCohesive


"""
    PowerYieldCohesive(; E, nu=0.0, ft, fc, zeta=5.0, wc, GF, ft_law=:hordijk, alpha=1.5, gamma=0.1, theta=1.5)

Constitutive model for cohesive elements with a power-lay yield surface ans softening in tension.  
The tensile softening branch is regularized through a measure of the
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
- `softening::Union{Symbol,PathFunction} = :hordijk`:  
  Softening law for post-peak tensile response. Options are:
  `:linear`, `:bilinear`, `:hordijk`, `:soft` or a custom function.
- `alpha::Real = 1.5`:  
  Parameter to control the shape of the yield surface (α > 0.5).
- `gamma::Real = 0.1`:  
  Parameter to control the residual shear strength (γ ≥ 0).
- `theta::Real = 1.5`:  
  Parameter to control the rate of reduction of shear strength (θ ≥ 0).
- `zeta::Real = 5.0`:  
  Factor to control elastic relative displacements in cohesive formulations (≥ 0).

# Returns
A `PowerYieldCohesive` object.

# Notes
- Either `wc` or `GF` must be provided. If only `GF` is given, `wc` is computed
  internally based on the chosen softening law.
- The frictional contribution is governed by `mu`.
- Normal and shear stiffnesses (`kn`, `ks`) are computed from the mechanical properties of
  the bulk material and the characteristic length `h` of the adjacent bulk elements.
"""
mutable struct PowerYieldCohesive<:Constitutive
    E ::Float64
    ν ::Float64
    fc::Float64
    ft::Float64
    wc::Float64
    softening::Symbol
    ft_fun::Union{PathFunction,Nothing}
    α::Float64
    γ::Float64
    θ::Float64
    βini::Float64
    ζ ::Float64

    function PowerYieldCohesive(; 
        E::Real = NaN,
        nu::Real = 0.0,
        fc::Real = NaN,
        ft::Real = NaN,
        wc::Real = NaN,
        GF::Real = NaN,
        softening::Union{Symbol,PathFunction} = :hordijk,
        alpha::Real = 1.5,
        gamma::Real = 0.1,
        theta::Real = 1.5,
        zeta::Real = 5.0,
    )

        @check E>0 "PowerYieldCohesive: Young's modulus E must be > 0. Got $(repr(E))."
        @check 0<=nu<0.5 "PowerYieldCohesive: Poisson ratio nu must be in the range [0, 0.5). Got $(repr(nu))."
        @check fc<0 "PowerYieldCohesive: Compressive strength fc must be < 0. Got $(repr(fc))."
        @check ft>0 "PowerYieldCohesive: Tensile strength ft must be > 0. Got $(repr(ft))."
        @check zeta>=0 "PowerYieldCohesive: Factor zeta must be non-negative. Got $(repr(zeta))."
        @check alpha > 0.5 "PowerYieldCohesive: alpha must be greater than 0.5. Got $(repr(alpha))."
        @check gamma >= 0.0 "PowerYieldCohesive: gamma must be non-negative. Got $(repr(gamma))."
        @check theta >= 0.0 "PowerYieldCohesive: theta must be non-negative. Got $(repr(theta))."
        @check softening in (:linear, :bilinear, :hordijk, :soft) || softening isa PathFunction "PowerYieldCohesive: Unknown softening model: $softening. Supported models are :linear, :bilinear, :hordijk, :soft or a custom PathFunction."

        ft_fun = nothing
        if softening isa PathFunction
            ft_fun = softening
            softening = :custom
            ft = softening(0.0)
            if softening.points[end][2] == 0.0
                wc = softening.points[end][1]
            else
                wc = Inf
            end
        else
            if isnan(wc)
                isnan(GF) && error("PowerYieldCohesive: wc or GF must be defined when using a predefined softening model")
                @check GF>0 "PowerYieldCohesive: GF must be positive. Got GF=$GF"

                if softening == :linear
                    wc = round(2*GF/ft, sigdigits=5)
                elseif softening == :bilinear
                    wc = round(5*GF/ft, sigdigits=5)
                elseif softening==:hordijk
                    wc = round(GF/(0.1947019536*ft), sigdigits=5)
                elseif softening==:soft
                    wc = round(GF/(0.1947019536*ft), sigdigits=5)
                end
            end
        end
        @check wc > 0.0 "PowerYieldCohesive: wc must be greater than zero"

        a     = (2*alpha*ft + alpha*fc - fc - √(alpha^2*fc^2 - 4*alpha^2*fc*ft + 4*alpha^2*ft^2 - 2*alpha*fc^2 + fc^2)) / (4*alpha-2)
        b     = √(alpha*(2*a-fc)*(ft-a))
        βini  = (b^2/ft^2)^alpha/(ft-a)

        return new(
            E, nu, fc, ft, wc, softening, ft_fun,
            alpha, gamma, theta, βini, zeta
        )
    end
end

@deprecate PowerYieldCrack PowerYieldCohesive


mutable struct PowerYieldCohesiveState<:IpState
    ctx::Context
    σ  ::Array{Float64,1} # stress
    w  ::Array{Float64,1} # relative displacements
    up ::Float64          # effective plastic relative displacement
    Δλ ::Float64          # plastic multiplier
    h  ::Float64          # characteristic length from bulk elements
    function PowerYieldCohesiveState(ctx::Context)
        this = new(ctx)
        this.σ   = zeros(ctx.ndim)
        this.w   = zeros(ctx.ndim)
        this.up  = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0
        return this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{PowerYieldCohesive}, ::Type{MechInterface}, ctx::Context) = PowerYieldCohesiveState


function paramsdict(mat::PowerYieldCohesive)
    mat = OrderedDict( string(field)=> getfield(mat, field) for field in fieldnames(typeof(mat)) )

    mat.softening in (:hordijk, :soft) && ( mat["GF"] = 0.1943*mat.ft*mat.wc )
    return mat
end


function beta(mat::PowerYieldCohesive, σmax::Float64)
    βini = mat.βini
    βres = mat.γ*βini
    return βres + (βini-βres)*(σmax/mat.ft)^mat.θ
end


function yield_func(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, σ::Array{Float64,1}, σmax::Float64)
    α  = mat.α
    β = beta(mat, σmax)
    ft = mat.ft
    if state.ctx.ndim == 3
        return β*(σ[1] - σmax) + ((σ[2]^2 + σ[3]^2)/ft^2)^α
    else
        return β*(σ[1] - σmax) + (σ[2]^2/ft^2)^α
    end
end


function yield_derivs(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, σ::Array{Float64,1}, σmax::Float64)
    α = mat.α
    β = beta(mat, σmax)
    ft = mat.ft

    if state.ctx.ndim == 3
        tmp = 2*α/ft^2*((σ[2]^2+σ[3]^2)/ft^2)^(α-1)
        σ[2]==σ[3]==0.0 && (tmp=0)
        return [ β , σ[2]*tmp, σ[3]*tmp ]
    else
        tmp = 2*α/ft^2*(σ[2]^2/ft^2)^(α-1)
        σ[2]==0.0 && (tmp=0)
        return [ β , σ[2]*tmp ]
    end
end


function potential_derivs(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, σ::Array{Float64,1})
    ndim = state.ctx.ndim
    if ndim == 3
        if σ[1] > 0.0 
            # G1:
            r = [ 2.0*σ[1], 2.0*σ[2], 2.0*σ[3]]
        else
            # G2:
            r = [ 0.0, 2.0*σ[2], 2.0*σ[3] ]
        end
        if r[1]==r[2]==r[3]==0.0
            r = [ 1.0, 0.0, 0.0]
        end
    else
        if σ[1] > 0.0 
            # G1:
            r = [ 2*σ[1], 2*σ[2]]
        else
            # G2:
            r = [ 0.0, 2*σ[2] ]
        end
        if r[1]==r[2]==0.0
            r = [ 1.0, 0.0]
        end
    end
    return r
end


function calc_σmax(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, up::Float64)
    if mat.softening == :linear
        if up < mat.wc
            a = mat.ft 
            b = mat.ft /mat.wc
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif mat.softening == :bilinear
        σs = 0.25*mat.ft
        ws = mat.wc*0.15
        if up < ws
            a  = mat.ft  
            b  = (mat.ft  - σs)/ws
        elseif up < mat.wc
            a  = mat.wc*σs/(mat.wc-ws)
            b  = σs/(mat.wc-ws)
        else
            a = 0.0
            b = 0.0
        end
        σmax = a - b*up
    elseif mat.softening == :hordijk
        if up < mat.wc
            z = (1 + 27*(up/mat.wc)^3)*exp(-6.93*up/mat.wc) - 28*(up/mat.wc)*exp(-6.93)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    elseif mat.softening == :soft
        dσmaxdup = 0.55
        a = 1.30837
        if up == 0.0
            z = 1.0
        elseif 0.0 < up < mat.wc
            x = up/mat.wc
            z = 1.0 - a^(1.0 - 1.0/x^dσmaxdup)
        else
            z = 0.0
        end
        σmax = z*mat.ft
    else
        σmax = mat.ft_fun(up)
    end

    return σmax
end


function deriv_σmax_upa(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, up::Float64)
    if mat.softening == :linear
        if up < mat.wc
            b = mat.ft /mat.wc
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softening == :bilinear
        ws = mat.wc*0.15
        σs = 0.25*mat.ft 
        if up < ws
            b  = (mat.ft  - σs)/ws
        elseif up < mat.wc
            b  = σs/(mat.wc-ws)
        else
            b = 0.0
        end
        dσmax = -b
    elseif mat.softening == :hordijk
        if up < mat.wc
            dz = ((81*up^2*exp(-6.93*up/mat.wc)/mat.wc^3) - (6.93*(1 + 27*up^3/mat.wc^3)*exp(-6.93*up/mat.wc)/mat.wc) - 0.02738402432/mat.wc)
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    elseif mat.softening == :soft
        m = 0.55
        a = 1.30837

        if up == 0.0
            dz = 0.0
        elseif up < mat.wc
            x = up/mat.wc
            dz =  -m*log(a)*a^(1-x^-m)*x^(-m-1)/mat.wc
        else
            dz = 0.0
        end
        dσmax = dz*mat.ft 
    else
        dσmax = derive(mat.ft_fun, up)
    end

    return dσmax
end


function calc_kn_ks(mat::PowerYieldCohesive, state::PowerYieldCohesiveState)
    kn = mat.E*mat.ζ/state.h
    G  = mat.E/(2*(1+mat.ν))
    ks = G*mat.ζ/state.h

    return kn, ks
end


function consistentD(mat::PowerYieldCohesive, state::PowerYieldCohesiveState)
    # numerical approximation
    # seems not to work under compressive loads

    ndim = state.ctx.ndim
    σmax = calc_σmax(mat, state, state.up)

    if state.Δλ == 0.0
        kn, ks = calc_kn_ks(mat, state)
        De = diagm([kn, ks, ks][1:ndim])
        return De
    # elseif σmax == 0.0 && state.w[1] >= 0.0
    #     kn, ks = calc_kn_ks(mat, state)
    #     De = diagm([kn, ks, ks][1:ndim])
    #     Dep = De*1e-4
    #     # Dep = De*1e-3
    #     return Dep
    end

    Dep = zeros(ndim, ndim)
    V = zeros(ndim)
    h = √eps()
    h = eps()^(1/3)

    # iteration for all w components
    for j in 1:ndim
        statej = copy(state)
        V[j] = 1.0
        Δw = h*V
        Δσ, succeeded = update_state(mat, statej, Δw)
        Dep[:,j] .= Δσ./h
        V[j] = 0.0
    end

    return Dep
end


function calcD(mat::PowerYieldCohesive, state::PowerYieldCohesiveState)
    # return consistentD(mat, state)

    ndim = state.ctx.ndim
    kn, ks = calc_kn_ks(mat, state)
    θ = mat.θ
    σmax = calc_σmax(mat, state, state.up)

    De = diagm([kn, ks, ks][1:ndim])

    if state.Δλ == 0.0  # Elastic 
        # @show "Elastic"
        return De
    elseif σmax == 0.0 && state.w[1] >= 0.0
        # @show "Plast"
        Dep = De*1e-4
        # Dep = De*1e-3
        return Dep
    else
        # @show "Elastic Pla"

        ft = mat.ft
        βini = mat.βini
        βres = mat.γ*βini
        β = beta(mat, σmax)
        dfdσmax = (βini-βres)/ft*(state.σ[1]-σmax)*θ*(σmax/ft)^(θ-1) - β

        r = potential_derivs(mat, state, state.σ)
        v = yield_derivs(mat, state, state.σ, σmax)
        m = deriv_σmax_upa(mat, state, state.up)  # ∂σmax/∂up

        if ndim == 3
            den = kn*r[1]*v[1] + ks*r[2]*v[2] + ks*r[3]*v[3] - dfdσmax*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      -kn*ks*r[1]*v[3]/den
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  -ks^2*r[2]*v[3]/den
                     -kn*ks*r[3]*v[1]/den        -ks^2*r[3]*v[2]/den        ks - ks^2*r[3]*v[3]/den ]
        else
            den = kn*r[1]*v[1] + ks*r[2]*v[2] - dfdσmax*m*norm(r)

            Dep = [   kn - kn^2*r[1]*v[1]/den    -kn*ks*r[1]*v[2]/den      
                     -kn*ks*r[2]*v[1]/den         ks - ks^2*r[2]*v[2]/den  ]
        end

        return Dep
    end
end


function calc_σ_up_Δλ(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, σtr::Array{Float64,1})
    ndim = state.ctx.ndim
    Δλ   = 0.0
    up   = 0.0
    σ    = zeros(ndim)
    σ0   = zeros(ndim)
    βini = mat.βini
    βres = mat.γ*βini
    θ    = mat.θ
    ft   = mat.ft
    
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
        up    = state.up + Δλ*norm_r
        σmax   = calc_σmax(mat, state, up)
        β      = beta(mat, σmax)

        f    = yield_func(mat, state, σ, σmax)
        dfdσ = yield_derivs(mat, state, σ, σmax)

        dfdσmax = (βini-βres)/ft*(σ[1]-σmax)*θ*(σmax/ft)^(θ-1) - β
        m = deriv_σmax_upa(mat, state, up)
        dσmaxdΔλ = m*(norm_r + Δλ*dot(r/norm_r, drdΔλ))
        dfdΔλ = dot(dfdσ, dσdΔλ) + dfdσmax*dσmaxdΔλ
        Δλ = Δλ - f/dfdΔλ
        
        if Δλ<=0 || isnan(Δλ) || i==maxits
            # return 0.0, state.σ, 0.0, failure("PowerYieldCohesive: failed to find Δλ")
            # switch to bissection method
            return calc_σ_up_Δλ_bis(mat, state, σtr)
        end

        if maximum(abs, σ-σ0) <= tol
            break
        end
        σ0 .= σ
    end

    return σ, up, Δλ, success()
end


function calc_σ_up(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, σtr::Array{Float64,1}, Δλ::Float64)
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

function calc_σ_up_Δλ_bis(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, σtr::Array{Float64,1})
    ndim    = state.ctx.ndim
    kn, ks  = calc_kn_ks(mat, state)
    De      = diagm([kn, ks, ks][1:ndim])
    r       = potential_derivs(mat, state, state.σ)

    # Δλ estimative
    Δλ0 = norm(σtr-state.σ)/norm(De*r)
    
    # find initial interval
    # a  = 0.0
    # b  = Δλ0

    ff(Δλ)  = begin
        # quantities at n+1
        σ, up = calc_σ_up(mat, state, σtr, Δλ)
        σmax = calc_σmax(mat, state, up)
        yield_func(mat, state, σ, σmax)
    end

    a, b, status = findrootinterval(ff, 0.0, Δλ0)
    failed(status) && return state.σ, 0.0, 0.0, status
    # @show a, b

    Δλ, status = findroot(ff, a, b, ftol=1e-5, method=:bisection)
    failed(status) && return state.σ, 0.0, 0.0, status

    σ, up = calc_σ_up(mat, state, σtr, Δλ)
    return σ, up, Δλ, success()  

end


function yield_func_from_Δλ(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, σtr::Array{Float64,1}, Δλ::Float64)
    ndim = state.ctx.ndim
    kn, ks = calc_kn_ks(mat, state)

    # quantities at n+1
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
    nr = norm(r)
    up = state.up + Δλ*nr
    
    σmax = calc_σmax(mat, state, up)
    f    = yield_func(mat, state, σ, σmax)

    return f
end



function update_state(mat::PowerYieldCohesive, state::PowerYieldCohesiveState, Δw::Array{Float64,1})

    ndim = state.ctx.ndim
    σini = copy(state.σ)

    kn, ks = calc_kn_ks(mat, state)
    De = diagm([kn, ks, ks][1:ndim])
    σmax = calc_σmax(mat, state, state.up)  

    if isnan(Δw[1]) || isnan(Δw[2])
        alert("PowerYieldCohesive: Invalid value for joint displacement: Δw = $Δw")
    end

    # σ trial and F trial
    σtr  = state.σ + De*Δw

    Ftr  = yield_func(mat, state, σtr, σmax)

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


function state_values(mat::PowerYieldCohesive, state::PowerYieldCohesiveState)
    ndim = state.ctx.ndim
    σmax = calc_σmax(mat, state, state.up)
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


function output_keys(mat::PowerYieldCohesive)
    return Symbol[:w, :σn, :up]
end