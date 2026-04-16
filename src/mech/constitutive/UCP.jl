# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export UCP

"""
    UCP(; E, nu, fc, epsc, eta=2.2, ft, GF, wc, p0=NaN, alpha=0.666, beta=1.15,
         ft_law=:hordijk, fc_law=:popovics, H=0.0)

Unified Concrete Plasticity model.

This constitutive model defines a three-invariant plasticity surface for concrete,
with a closed cap in compression and fracture-energy regularization in tension.
It combines elastic isotropy, nonlinear hardening/softening in compression, and
tension softening controlled by the fracture energy.

# Keyword arguments
- `E::Real`:  
  Young’s modulus (must be > 0).
- `nu::Real`:  
  Poisson’s ratio (0 ≤ ν < 0.5).
- `fc::Real`:  
  Uniaxial compressive strength (< 0).
- `epsc::Real`:  
  Strain at the compressive peak (< 0).
- `eta::Real = 2.2`:  
  Shape parameter for the compression hardening/softening curve (eta > 1).
- `ft::Real`:  
  Uniaxial tensile strength (> 0).
- `GF::Real`:  
  Tensile fracture energy (> 0). Can be given alternatively to `wc`.
- `wc::Real`:  
  Critical crack opening displacement (≥ 0). Can be given alternatively to `GF`.
- `p0::Real = NaN`:  
  Elastic limit in isotropic compression. If not given, computed internally
  from `fc` and `beta`.
- `alpha::Real = 0.666`:  
  Shape parameter of the meridional section (0.2 < α ≤ 1.0).
- `beta::Real = 1.15`:  
  Factor relating biaxial to uniaxial compressive strength (1 ≤ β ≤ 1.5).
- `H::Real = 0.0`:  
  Plastic modulus for isotropic compression (≥ 0).

# Returns
A `UCP` material object that can be attached to mechanical bulk elements
for 2D (plane strain) or 3D analyses. For plane stress analyses, the
integration is performed decoupled with Δε33=0 and subsequently updated according
to the resulting stress increment.

# Notes
- The yield function uses a nonlinear function that matches `fc` and `ft`.
- The tensile law is regularized through `GF` and `wc` to ensure energy
  dissipation is independent of element size.
- The compressive response follows a nonlinear curve defined by `fc`, `epsc`, and `eta`.
- The surface excentricity is computed internally to match the biaxial strength.
- The surface section follows the Willam-Warnke ellipsoidal shape.
- The cap position is adjusted by `beta` and `p0`.
"""
mutable struct UCP<:Constitutive
    E::Float64
    ν::Float64
    fc::Float64
    εc::Float64
    η::Float64
    ft::Float64
    wc::Float64
    fb::Float64
    ξc0::Float64
    ξt0::Float64
    ft_law::Symbol
    ft_fun::Union{Nothing,AbstractSpline}
    fc_law::Symbol
    fc_fun::Union{Nothing,AbstractSpline}
    α::Float64
    e::Float64
    H::Float64

    function UCP(;
        E::Real    = NaN,
        nu::Real   = NaN,
        alpha::Real= 0.666,
        fc::Real   = NaN,
        epsc::Real = NaN,
        eta::Real  = 2.2,
        ft::Real   = NaN,
        GF::Real   = NaN,
        wc::Real   = NaN,
        beta::Real = 1.15,
        p0::Real   = NaN,
        ft_law     = :hordijk,
        fc_law     = :popovics,
        H::Real    = 0.0,
    )
        @check E>0 "UCP: Young's modulus E must be > 0. Got $E."
        @check 0<=nu<0.5 "UCP: Poisson's ratio nu must be in the range [0, 0.5). Got $nu."
        @check 0.2<alpha<=1.0 "UCP: Curvature coefficient alpha must be in the range (0.2, 1.0]. Got $alpha."
        @check 1<=beta<=1.5 "UCP: Factor beta must be in the range [1.0, 1.5]. Got $beta."
        @check eta>1 "UCP: Shape parameter eta must be > 1. Got $eta."

        @check ft>0 "UCP: Tensile strength ft must be > 0. Got $ft."
        @check H>=0 "UCP: Plastic modulus H must be >= 0. Got $H."

        wc, ft_law, ft_fun, status = setup_tensile_strength(ft, GF, wc, ft_law)
        failed(status) && throw(ArgumentError("UCP: " * status.message))


        fc_law, fc_fun, status = setup_compressive_strength(E, fc, epsc, fc_law)
        failed(status) && throw(ArgumentError("UCP: " * status.message))

        if fc_law isa AbstractSpline
            # fc_fun = fc_law
            # fc_law = :custom
            fc = fc_law(0.0) # TODO: get the maximun value ?
        end

        @check abs(epsc)>abs(fc)/E "UCP: epsc should be greater than fc/E."

        α = alpha
        β = beta

        # value of exentricity to match fb in a biaxial trajectory, assuming the state when ξt=0
        e  = β/(2*β)^α
        fb = β*fc
        
        if isnan(wc)
            @check GF>0 "UCP: Fracture energy GF must be > 0. Got $(repr(GF))."
            wc = round(GF/(0.1947*ft), sigdigits=5)  # inverse of Hordijk approximation
            notify("UCP: Using Hordijk's approximation wc=$(repr(wc)).")
        else
            @check wc>=0 "UCP: Critical crack opening wc must be >= 0. Got $(repr(wc))."
        end
        wc > 1e-5 || notify("UCP: Warning: very low value of wc=$(repr(wc)).")

        if isnan(p0)
            ξb  = 2*fb/√3
            ξc0 = 1.5*ξb
        else
            @check p0<0 "UCP: Elastic limit in isotropic compression p0 must be < 0. Got $(repr(p0))."
            ξc0 = √3*p0
        end

        fc0 = 0.4*fc
        ft0 = ft
        Ω   = (-ft0/(fc0*e))^(1/α)
        ξt0 = 1/√3*(fc0*Ω - ft0)/(Ω-1)

        return new(E, nu, fc, epsc, eta, ft, wc, fb, ξc0, ξt0, ft_law, ft_fun, fc_law, fc_fun, α, e, H)
    end
end


mutable struct UCPState<:ConstState
    ctx::Context
    σ  ::Vec6
    ε  ::Vec6
    εtp::Float64
    εcp::Float64
    εvp::Float64
    Δλ ::Float64
    h  ::Float64
    function UCPState(ctx::Context)
        this     = new(ctx)
        this.σ   = zeros(Vec6)
        this.ε   = zeros(Vec6)
        this.εtp = 0.0 # plastic strain in tension
        this.εcp = 0.0 # plastic strain in compression
        this.εvp = 0.0 # plastic volumetric strain in compression
        this.Δλ  = 0.0 # increment of plastic multiplier
        this.h   = 0.0
        this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{UCP}, ::Type{MechSolid}) = UCPState


function calc_θ(mat::UCP, ρ::Float64, j3::Float64)

    ρtol = 1e-8*abs(mat.fc)
    ctol = 1e-6
    
    if ρ < ρtol
        # hydrostatic axis / apex
        θ = 0.0
    else
        c = clamp(3*√6*j3/ρ^3, -1.0, 1.0)
        if 1 - abs(c) < ctol 
            # meridians
            θ = c > 0 ? 0.0 : π/3
        else
            θ = acos(c)/3
        end
    end

    return θ
end


function calc_rθ(mat::UCP, θ::Float64)
    e = mat.e

    rnum   = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden   = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    r      = rnum/rden

    return r
end


function calc_rξ(mat::UCP, ξc::Float64, ξt::Float64, ξ::Float64)
    fca = abs(mat.fc)
    return spow((ξt-ξ)/fca, mat.α)
end


function calc_rc(mat::UCP, ξc::Float64, ξ::Float64)
    ξb = 2*mat.fb/√3
    ξ>=ξb && return 1.0
    ξ<ξc  && return 0.0
    return √(1 - ((ξb-ξ)/(ξb-ξc))^2)
end


function is_apex_state(mat::UCP, ξ::Float64, ξt::Float64)
    ξtol = √eps(Float64)*max(abs(mat.fc), abs(ξt), 1.0)
    return ξ >= ξt - ξtol
end


function calc_fc(mat::UCP, εcp::Float64)
    fc0 = 0.4*mat.fc
    fcr = 0.1*mat.fc

    return calc_compressive_strength(mat, fc0, fcr, εcp)
end

function calc_fc_derivative(mat::UCP, εcp::Float64)
    fc0 = 0.4*mat.fc
    fcr = 0.1*mat.fc

    return calc_compressive_strength_derivative(mat, fc0, fcr, εcp)
end


function calc_ft(mat::UCP, w::Float64)
    return calc_tensile_strength(mat, w)
end


function calc_ft_derivative(mat::UCP, w::Float64)
    ∂ft∂w = calc_tensile_strength_derivative(mat, w)
    Hcap  = isfinite(mat.wc) ? -mat.ft/(0.5*mat.wc) : -Inf
    return max(∂ft∂w, Hcap)
end


function calc_ξc_ξt_m(mat::UCP, h::Float64, εtp::Float64, εcp::Float64, εvp::Float64)
    α  = mat.α
    w  = εtp*h

    ft = calc_ft(mat, w)
    fc = calc_fc(mat, εcp)
    
    # p = p0 + H*εvp  -> ξc = √3*p0 + √3*H*εvp
    ξc = mat.ξc0 + √3*mat.H*εvp # hardening in isotropic compression
    ξt = mat.ξt0*ft/mat.ft # this may give wrong results for the peak tensile stress

    @assert ξc<0
    @assert ξc<fc/√3

    fca = abs(mat.fc)
    m  = -√(2/3)*fc*((ξt - fc/√3)/fca)^-α  # fc is current fc value
    @assert m>0

    return ξc, ξt, m
end


function yield_func(mat::UCP, h::Float64, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    # f(σ) = ρ - rθ⋅rc⋅rξ⋅κ

    i1, j2 = tr(σ), J2(σ)

    ξ = i1/√3
    ρ = √(2*j2)

    ξc, ξt, m = calc_ξc_ξt_m(mat, h, εtp, εcp, εvp)
    θ = calc_θ(mat, ρ, J3(σ))
    rθ = calc_rθ(mat, θ)
    rc = calc_rc(mat, ξc, ξ)
    rξ = calc_rξ(mat, ξc, ξt, ξ)

    return ρ - rθ*rc*rξ*m
end


function yield_derivs(mat::UCP, h::Float64, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    i1, j2, j3 = tr(σ), J2(σ), J3(σ)
    ξc, ξt, m = calc_ξc_ξt_m(mat, h, εtp, εcp, εvp)
    
    α  = mat.α
    ξb = 2*mat.fb/√3
    ξ  = i1/√3
    ρ  = √(2*j2)
    θ  = calc_θ(mat, ρ, j3)
    rθ = calc_rθ(mat, θ)
    rc = calc_rc(mat, ξc, ξ)
    rξ = calc_rξ(mat, ξc, ξt, ξ)

    fca = abs(mat.fc)

    dfdrc = -rθ*rξ*m
    dfdrξ = -rθ*rc*m

    # ∂f/∂εtp
    fc     = calc_fc(mat, εcp)
    drξdξt = ξt - ξ > 0.0 ? α/fca * ((ξt - ξ)/fca)^(α-1) : 0.0
    ∂f∂m   = -rθ*rc*rξ
    ∂m∂ξt  = √(2/3)*α*fc/fca*((ξt - fc/√3)/fca)^(-α-1)
    dfdξt  = dfdrξ*drξdξt + ∂f∂m*∂m∂ξt
    ∂ξt∂ft = mat.ξt0/mat.ft
    w      = εtp*h
    ∂ft∂w  = calc_tensile_strength_derivative(mat, w)
    ∂w∂εtp = h
    ∂f∂εtp = dfdξt*∂ξt∂ft*∂ft∂w*∂w∂εtp

    # ∂f/∂εcp
    ∂f∂m    = -rθ*rc*rξ
    fc      = calc_fc(mat, εcp)
    ∂mdfc   = -√(2/3) * ((ξt - fc/√3)/fca)^-α  -  α*√2/3*fc/fca * ((ξt - fc/√3)/fca)^(-α-1)
    dfcdεcp = calc_fc_derivative(mat, εcp)
    ∂f∂εcp  = ∂f∂m*∂mdfc*dfcdεcp

    # ∂f/∂εvp
    if mat.H!=0.0
        ∂rc∂ξc  = ξc < ξ < ξb ? -(ξb-ξ)^2/(ξb-ξc)^3/√(1 - ((ξb-ξ)/(ξb-ξc))^2) : 0.0
        ∂f∂rc   = -rθ*rξ*m
        ∂f∂ξc   = ∂f∂rc*∂rc∂ξc
        ∂ξc∂εvp = √3*mat.H
        ∂f∂εvp  = ∂f∂ξc*∂ξc∂εvp
    else
        ∂f∂εvp = 0.0
    end
    
    # check apex condition
    ξ = tr(σ)/√3
    is_apex_state(mat, ξ, ξt) && return √3/3*I2, ∂f∂εtp, ∂f∂εcp, ∂f∂εvp
    

    ξb = 2*mat.fb/√3
    rc = calc_rc(mat, ξc, ξ)
    rξ = calc_rξ(mat, ξc, ξt, ξ)

    # f derivative w.r.t. σ:
    ∂f∂ρ  = 1.0
    dfdrc = -rθ*rξ*m
    dfdrξ = -rθ*rc*m
    drcdξ = ξc<ξ<ξb ? (ξb-ξ)/(ξb-ξc)^2/√(1-((ξb-ξ)/(ξb-ξc))^2) : 0.0

    drξdξ = ξ < ξt ? -α/fca * abs((ξt-ξ)/fca)^(α-1) : 0.0
    
    ∂f∂ξ  = dfdrc*drcdξ + dfdrξ*drξdξ
    dξdσ = √3/3*I2

    θ  = calc_θ(mat, ρ, j3)
    rθ = calc_rθ(mat, θ)

    use_lode_derivative = !( θ == 0.0 || θ == π/3 ) # avoid singularity at meridians (θ=0, π/3) and apex (ρ=0)

    if use_lode_derivative
        s = dev(σ)
        e = mat.e
        rnum = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
        rden = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2

        ∂rθ∂numdθ = (2*sin(2*θ)*(2*e-1)*(e^2-1))/√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e) - 2*(1 - e^2)*sin(θ)
        ∂rθ∂endθ = 4*sin(2*θ)*(e^2-1)
        ∂rθ∂θ    = (∂rθ∂numdθ*rden - rnum*∂rθ∂endθ)/rden^2
        ∂f∂rθ    = -rc*rξ*m
        ∂f∂θ     = ∂f∂rθ*∂rθ∂θ
        ∂ρ∂σ     = s/ρ
        ∂s∂σ     = Psd
        adj_s    = adj(s)
        ∂θ∂s     = -√6*(adj_s/ρ^3 - 3*s*j3/ρ^5)/√abs(1 - 54*j3^2/ρ^6)
        ∂θ∂σ     = ∂s∂σ*∂θ∂s

        ∂f∂σ = ∂f∂ρ*∂ρ∂σ + ∂f∂ξ*dξdσ + ∂f∂θ*∂θ∂σ
    else
        s    = dev(σ)
        ∂ρ∂σ = ρ > 0.0 ? s/ρ : zero(s)
        ∂f∂σ = ∂f∂ρ*∂ρ∂σ + ∂f∂ξ*dξdσ
    end

    return ∂f∂σ, ∂f∂εtp, ∂f∂εcp, ∂f∂εvp
end


function potential_derivs(mat::UCP, h::Float64, σ::AbstractArray, εtp::Float64)
    # g(σ) = ρ^2 - 4 tan^2(ψ) (ξ_t - ξ_f'_c) (ξ_t - ξ) = 0
    # g(σ) = ρ^2 - 4 χ^2 (ξ_t - ξ_f'_c) (ξ_t - ξ) = 0
    ξfc = mat.fc/√3

    w  = εtp*h
    ft = calc_ft(mat, w)
    ξt = mat.ξt0*ft/mat.ft
    ψ  = 0.2
    χ  = tan(ψ)
    ρ  = √(2*J2(σ))

    dgdξ = 4*χ^2*(ξt - ξfc)
    dξdσ = √3/3*I2
    
    ρtol = 1e-8*abs(mat.fc)
    if ρ > ρtol
        s    = dev(σ)
        dgdρ = 2*ρ
        ∂ρ∂σ = s/ρ
        ∂g∂σ = dgdρ*∂ρ∂σ + dgdξ*dξdσ
    else
        # hydrostatic axis / apex
        ∂g∂σ = dgdξ*dξdσ
    end

    return ∂g∂σ
end


function ucp_plastic_flow_invariant_rates(∂g∂σ::Vec6)
    # Recover the principal values analytically from invariants, then apply
    # the original positive/negative spectral split.
    Λ1, Λ2, Λ3 = eigvals(∂g∂σ)

    # rate_εtp = max(Λ1, Λ2, Λ3, 0.0)
    rate_εtp = sqrt(max(Λ1, 0.0)^2 + max(Λ2, 0.0)^2 + max(Λ3, 0.0)^2) # what about p-norm with p around 5
    rate_εcp = sqrt(min(Λ1, 0.0)^2 + min(Λ2, 0.0)^2 + min(Λ3, 0.0)^2)
    rate_εvp = abs(min(Λ1, 0.0) + min(Λ2, 0.0) + min(Λ3, 0.0))

    return rate_εtp, rate_εcp, rate_εvp
end


function calcD(mat::UCP, state::UCPState)
    De  = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    h = state.h

    state.Δλ==0.0 && return De

    ∂f∂σ, ∂f∂εtp, ∂f∂εcp = yield_derivs(mat, h, state.σ, state.εtp, state.εcp, state.εvp)
    ∂g∂σ = potential_derivs(mat, h, state.σ, state.εtp)
    rate_εtp, rate_εcp, _ = ucp_plastic_flow_invariant_rates(∂g∂σ)

    De_dgdσ = De*∂g∂σ
    denom = ∂f∂σ'*De_dgdσ - ∂f∂εcp*rate_εcp - ∂f∂εtp*rate_εtp
    Dep = De - De_dgdσ*∂f∂σ'*De / denom

    return Dep
end


function plastic_update(mat::UCP, state::UCPState, cstate::UCPState, σtr::Vec6)
    maxits = 50
    tol    = mat.ft*1e-5
    h      = state.h
    ∂g∂σ   = potential_derivs(mat, h, cstate.σ, cstate.εtp) # ketp frozen
    De     = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    Δλ     = eps()

    σ  = σtr - Δλ*(De*∂g∂σ)

    εcp = cstate.εcp
    εtp = cstate.εtp
    εvp = cstate.εvp

    f = yield_func(mat, h, σ, εtp, εcp, εvp)
    ω = 1.0 # initial damping

    # NR iterations
    for i in 1:maxits
        
        ∂f∂σ, ∂f∂εtp, ∂f∂εcp, ∂f∂εvp = yield_derivs(mat, h, σ, εtp, εcp, εvp)
        rate_εtp, rate_εcp, rate_εvp = ucp_plastic_flow_invariant_rates(∂g∂σ)

        ∂f∂Δλ = -∂f∂σ'*De*∂g∂σ + ∂f∂εcp*rate_εcp + ∂f∂εtp*rate_εtp + ∂f∂εvp*rate_εvp

        function eval_f(Δλtest::Float64)
            σt   = σtr - Δλtest * (De * ∂g∂σ)
            εtpt = cstate.εtp + Δλtest * rate_εtp
            εcpt = cstate.εcp + Δλtest * rate_εcp
            εvpt = cstate.εvp + Δλtest * rate_εvp
            return yield_func(mat, h, σt, εtpt, εcpt, εvpt)
        end

        # Newton step direction
        Δλmin = max(Δλ - ω * f / ∂f∂Δλ, 0.0)
        isfinite(Δλmin) || break
        fmin = eval_f(Δλmin)

        # Backtracking line search: required since ∂f∂Δλ is not the exact derivative
        for ω in 0.9:-0.1:0.3
            Δλtr = Δλ - ω * f / ∂f∂Δλ
        
            Δλtr > 0.0 || continue

            ftr = eval_f(Δλtr)
            isfinite(ftr) || continue

            if abs(ftr) < abs(fmin)
                Δλmin = Δλtr
                fmin  = ftr
            end
        end

        Δλ = Δλmin
        f  = fmin
 
        isfinite(Δλ) || break

        σ  = σtr - Δλ*(De*∂g∂σ)

        εtp = cstate.εtp + Δλ*rate_εtp
        εcp = cstate.εcp + Δλ*rate_εcp
        εvp = cstate.εvp + Δλ*rate_εvp

        f = yield_func(mat, h, σ, εtp, εcp, εvp)
        
        if abs(f) < tol
            Δλ < 0.0 && break

            w  = εtp * state.h
            ft = calc_ft(mat, w)
            fc = calc_fc(mat, εcp)
            abs(fc*mat.e/ft) > 1.1 || break
            @assert εcp >= 0.0
            @assert εtp >= 0.0
            @assert εvp >= 0.0

            state.σ   = σ
            state.εtp = εtp
            state.εcp = εcp
            state.εvp = εvp
            state.Δλ  = Δλ
            
            return success()
        end

    end

    return failure("UCP: plastic update failed")
end


function update_state(mat::UCP, state::UCPState, cstate::UCPState, Δε::AbstractArray)

    De   = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    h    = state.h
    σtr  = cstate.σ + De*Δε
    ftr  = yield_func(mat, h, σtr, cstate.εtp, cstate.εcp, cstate.εvp)

    tol = 0.001

    # Elastic step
    if ftr < tol
        # elastic
        state.Δλ = 0.0
        state.σ  = σtr
        Δσ       = state.σ - cstate.σ
        state.ε  = cstate.ε + Δε
        return Δσ, success()
    end

    # Plastic update
    status = plastic_update(mat, state, cstate, σtr)
    failed(status) && return state.σ, status

    Δσ = state.σ - cstate.σ
    
    # Update Δεzz for plane stress (since update_state in uncoupled with Δεzz for plane stress)
    if state.ctx.stress_state == :plane_stress
        ∂g∂σ = potential_derivs(mat, h, state.σ, state.εtp)
        
        Δε33e = -(mat.ν / mat.E) * (Δσ[1] + Δσ[2])
        Δεp = state.Δλ * ∂g∂σ

        Δε = Vec6(Δε[1], Δε[2], Δε33e + Δεp[3], 0.0, 0.0, Δε[6])
        # σtr  = cstate.σ + De*Δε
        # status = plastic_update(mat, state, cstate, σtr)
        # failed(status) && return state.σ, status

        Δσ = state.σ - cstate.σ
    end

    state.ε = cstate.ε + Δε

    return Δσ, success()
end


function state_values(mat::UCP, state::UCPState)
    σ, ε  = state.σ, state.ε
    h = state.h
    ρ = √(2*J2(σ))
    ξ = tr(σ)/√3
    θ  = calc_θ(mat, ρ, J3(σ))

    w  = state.εtp*state.h
    ft = calc_ft(mat, w)
    fc = calc_fc(mat, state.εcp)

    ξc, ξt, m = calc_ξc_ξt_m(mat, h, state.εtp, state.εcp, state.εvp)
    # rc = calc_rc(mat, ξc, ξ)
    # rξ = calc_rξ(mat, ξc, ξt, ξ)

    vals_d = stress_strain_dict(σ, ε, state.ctx.stress_state)

    vals_d[:εcp] = state.εcp
    vals_d[:εtp] = state.εtp
    vals_d[:ξ]   = ξ
    vals_d[:ρ]   = ρ
    vals_d[:θ]   = θ
    vals_d[:fc]  = fc
    vals_d[:ft]  = ft
    vals_d[:ξc]  = ξc
    vals_d[:ξt]  = ξt
    vals_d[:m]   = m
    # vals_d[:r]   = r
    # vals_d[:rξ]  = rξ
    # vals_d[:rc]  = rc
    # vals_d[:ξt]  = 2*mat.fb/√3
    # vals_d[:fcb] = abs(mat.fc)

    return vals_d
end
