# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export UCP

"""
    UCP(; E, nu, fc, epsc, eta=4, ft, GF, wc, p0, alpha=0.666, beta=1.15, H=0.0)

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
- The tensile law is regularized through `GF` and `wc` to ensure energy
  dissipation is independent of element size.
- The compressive response follows a nonlinear curve defined by `fc`, `epsc`, and `eta`.
- The cap position is adjusted by `beta` and `p0`.
- The surface excentricity is computed internally to match the biaxial strength.
- The surface section follows the Willam-Warnke ellipsoidal shape.
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
    ξa0::Float64
    ξc0::Float64
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
        fc_law     = :default,
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
        
        # value of exentricity to match fb in a biaxial trajectory, assuming the state when ξc=0
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
            ξa0 = 1.5*ξb
        else
            @check p0<0 "UCP: Elastic limit in isotropic compression p0 must be < 0. Got $(repr(p0))."
            ξa0 = √3*p0
        end

        fc0 = 0.4*fc
        ft0 = ft
        Ω   = (-ft0/(fc0*e))^(1/α)
        ξc0 = 1/√3*(fc0*Ω - ft0)/(Ω-1)

        return new(E, nu, fc, epsc, eta, ft, wc, fb, ξa0, ξc0, ft_law, ft_fun, fc_law, fc_fun, α, e, H)
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
# compat_state_type(::Type{UCP}, ::Type{MechBulk}) = ctx.stress_state!=:plane_stress ? UCPState : error("UCP: This model is not compatible with planestress")
compat_state_type(::Type{UCP}, ::Type{MechBulk}) = UCPState


function calc_θ(::UCP, σ::Vec6)
    j2 = J2(σ)
    if j2==0.0
        θ = 0.0
    else
        norm_s = √(2*j2)
        det_s  = J3(σ)
        θ      = 1/3*acos( clamp(3*√6*det_s/norm_s^3, -1.0, 1.0) )
    end
    return θ
end


function calc_rθ(mat::UCP, σ::Vec6)
    e = mat.e
    θ = calc_θ(mat, σ)

    rnum   = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden   = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    r      = rnum/rden

    return r
end


function calc_rξ(mat::UCP, ξa::Float64, ξc::Float64, ξ::Float64)
    abs_fc = abs(mat.fc)
    return spow((ξc-ξ)/abs_fc, mat.α)
end


function calc_rχ(mat::UCP, ξa::Float64, ξ::Float64)
    ξb = 2*mat.fb/√3
    ξ>=ξb && return 1.0
    ξ<ξa  && return 0.0
    return √(1 - ((ξb-ξ)/(ξb-ξa))^2)
end


function calc_fc(mat::UCP, εcp::Float64)
    fc0 = 0.4*mat.fc
    fcr = 0.1*mat.fc

    return calc_compressive_strength(mat, fc0, fcr, εcp)
end


function calc_ft(mat::UCP, w::Float64)
    # return calc_tensile_strength(mat, w)
    return max(calc_tensile_strength(mat, w), 0.01*mat.ft)
end


function calc_ξa_ξc_κ(mat::UCP, h::Float64, εtp::Float64, εcp::Float64, εvp::Float64)
    α  = mat.α
    w  = εtp*h

    ft = calc_ft(mat, w)
    fc = calc_fc(mat, εcp)
    
    # p = p0 + H*εvp  -> ξa = √3*p0 + √3*H*εvp
    ξa = mat.ξa0 + √3*mat.H*εvp # hardening in isotropic compression
    ξc = mat.ξc0*ft/mat.ft

    @assert ξa<0
    @assert ξa<fc/√3

    abs_fc = abs(mat.fc)
    κ  = -√(2/3)*fc*((ξc - fc/√3)/abs_fc)^-α  # fc is current fc value
    @assert κ>0

    return ξa, ξc, κ
end


function yield_func(mat::UCP, h::Float64, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    # f(σ) = ρ - rθ⋅rc⋅rξ⋅κ

    i1, j2 = tr(σ), J2(σ)

    ξ = i1/√3
    ρ = √(2*j2)

    ξa, ξc, κ = calc_ξa_ξc_κ(mat, h, εtp, εcp, εvp)
    rθ = calc_rθ(mat, σ)
    rχ = calc_rχ(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξa, ξc, ξ)

    return ρ - rθ*rχ*rξ*κ
end


function yield_derivs(mat::UCP, h::Float64, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)

    # ∂f/∂εtp, ∂f/∂εcp
    f_εcp  = εcp -> yield_func(mat, h, σ, εtp, εcp, εvp)
    ∂f∂εcp = derive(f_εcp, εcp)

    f_εtp  = εtp -> yield_func(mat, h, σ, εtp, εcp, εvp)
    ∂f∂εtp = derive(f_εtp, εtp)

    ξa, ξc, κ = calc_ξa_ξc_κ(mat, h, εtp, εcp, εvp)
    
    # check apex condition
    ξ = tr(σ)/√3
    ξ >= ξc && return √3/3*I2, ∂f∂εtp, ∂f∂εcp
    
    # deviatoric derivatives
    j2     = J2(σ)
    ρ      = √(2*j2)
    s      = dev(σ)
    det_s  = J3(σ)
    adj_s  = det_s*inv(s)
    norm_s = ρ
    
    # θ and derivatives
    e        = mat.e
    θ        = 1/3*acos( clamp(3*√6*det_s/norm_s^3, -1.0, 1.0) )
    rnum     = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden     = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    rθ       = rnum/rden
    drθnumdθ = (2*sin(2*θ)*(2*e-1)*(e^2-1))/√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e) - 2*(1 - e^2)*sin(θ) 
    drθdendθ = 4*sin(2*θ)*(e^2-1)
    drθdθ    = (drθnumdθ*rden - rnum*drθdendθ)/rden^2

    ϵ = 1e-10
    dθds = -√6*(adj_s/ρ^3 - 3*s*det_s/ρ^5)/√abs(1 - 54*det_s^2/ρ^6 + ϵ ) # denominator approaches zero at meridians (θ=0, π/3)

    ξb = 2*mat.fb/√3
    rχ = calc_rχ(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξa, ξc, ξ)

    # f derivative w.r.t. σ:
    dfdρ  = 1.0
    dfdrχ = -rθ*rξ*κ
    dfdrξ = -rθ*rχ*κ
    drχdξ = ξa<ξ<ξb ? (ξb-ξ)/(ξb-ξa)^2/√(1-((ξb-ξ)/(ξb-ξa))^2) : 0.0

    α = mat.α
    abs_fc = abs(mat.fc)
    # drξdξ = -α/abs_fc * abs((ξc-ξ)/abs_fc)^(α-1)
    # drξdξ =  -α/abs_fc * abs((ξc-ξ)/abs_fc)^(α-1)*sign((ξc-ξ)/abs_fc)
    drξdξ = ξ < ξc ? -α/abs_fc * abs((ξc-ξ)/abs_fc)^(α-1) : 0.0

    
    dfdξ  = dfdrχ*drχdξ + dfdrξ*drξdξ
    dfdrθ = -rχ*rξ*κ
    dfdθ  = dfdrθ*drθdθ

    dρdσ = s/norm(s)
    dξdσ = √3/3*I2
    dsdσ = Psd
    dθdσ = dsdσ*dθds

    ∂f∂σ = dfdρ*dρdσ + dfdξ*dξdσ + dfdθ*dθdσ
    
    return ∂f∂σ, ∂f∂εtp, ∂f∂εcp
end


function potential_derivs(mat::UCP, h::Float64, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    # g(σ) = ρ - rc⋅rξ⋅κ
    abs_fc = abs(mat.fc)

    i1 = tr(σ)
    ξ  = i1/√3
    
    ξa, ξc, κ = calc_ξa_ξc_κ(mat, h, εtp, εcp, εvp)
    ξ >= ξc && return √3/3*I2 # apex

    α = mat.α
    s = dev(σ)
    ρ = norm(s) + eps()

    ξb = 2*mat.fb/√3
    rχ = calc_rχ(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξa, ξc, ξ)

    dgdrχ = -rξ*κ
    dgdrξ = -rχ*κ
    drχdξ = ξa<ξ<ξb ? (ξb-ξ)/(ξb-ξa)^2/√(1-((ξb-ξ)/(ξb-ξa))^2) : 0.0
    drξdξ = ξ < ξc ? -α/abs_fc * abs((ξc-ξ)/abs_fc)^(α-1) : 0.0
    dgdξ  = dgdrχ*drχdξ + dgdrξ*drξdξ

    dξdσ = √3/3*I2
    dgdρ = 1.0

    ∂g∂σ = s/ρ + dgdξ*dξdσ

    # Near apex check    
    if ξc < 0.2*mat.ft && ξ >= 0.0
        ∂g∂σ = √3/3*I2
    end

    return ∂g∂σ
end


function calcD(mat::UCP, state::UCPState)
    De  = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    h = state.h

    state.Δλ==0.0 && return De

    ∂f∂σ, ∂f∂εtp, ∂f∂εcp = yield_derivs(mat, h, state.σ, state.εtp, state.εcp, state.εvp)
    ∂g∂σ = potential_derivs(mat, h, state.σ, state.εtp, state.εcp, state.εvp)

    Λ = eigvals(∂g∂σ)
    Λ1, Λ2, Λ3 = Λ

    max_Λp  = max(Λ1, Λ2, Λ3, 0.0)
    norm_Λn = (min(Λ1, 0.0)^2 + min(Λ2, 0.0)^2 + min(Λ3, 0.0)^2)^0.5
    
    De_dgdσ = De*∂g∂σ
    denom = ∂f∂σ'*De_dgdσ - ∂f∂εcp*norm_Λn - ∂f∂εtp*max_Λp
    Dep = De - De_dgdσ*∂f∂σ'*De / denom

    return Dep
end


function plastic_update(mat::UCP, state::UCPState, cstate::UCPState, σtr::Vec6)
    maxits = 50
    tol    = mat.ft*1e-4
    h      = state.h
    ∂g∂σ   = potential_derivs(mat, h, cstate.σ, cstate.εtp, cstate.εcp, cstate.εvp)
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
        
        ∂f∂σ, ∂f∂εtp, ∂f∂εcp = yield_derivs(mat, h, σ, εtp, εcp, εvp)
        ∂g∂σ = potential_derivs(mat, h, σ, εtp, εcp, εvp)
        Λ  = eigvals(∂g∂σ)
        Λ1, Λ2, Λ3 = Λ
        
        max_Λp  = max(Λ1, Λ2, Λ3, 0.0)
        norm_Λn = (min(Λ1, 0.0)^2 + min(Λ2, 0.0)^2 + min(Λ3, 0.0)^2)^0.5
        sum_Λn  = abs(min(Λ1, 0.0) + min(Λ2, 0.0) + min(Λ3, 0.0) )

        ∂f∂Δλ   = -∂f∂σ'*De*∂g∂σ + ∂f∂εcp*norm_Λn + ∂f∂εtp*max_Λp

        function eval_f(Δλtest::Float64)
            σt   = σtr - Δλtest * (De * ∂g∂σ)
            εtpt = cstate.εtp + Δλtest * max_Λp
            εcpt = cstate.εcp + Δλtest * norm_Λn
            εvpt = cstate.εvp + Δλtest * sum_Λn
            return yield_func(mat, h, σt, εtpt, εcpt, εvpt)
        end

        # Newton step direction
        Δλmin = max(Δλ - ω * f / ∂f∂Δλ, 0.0)
        fmin  = eval_f(Δλmin)

        # Backtracking line search
        for ω in 0.9:-0.1:0.3
            Δλtr = Δλ - ω * f / ∂f∂Δλ
        
            Δλtr > 0.0 || continue

            ftr   = eval_f(Δλtr)
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

        εtp = cstate.εtp + Δλ*max_Λp
        εcp = cstate.εcp + Δλ*norm_Λn
        εvp = cstate.εvp + Δλ*sum_Λn
        
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



function plastic_update_num(mat::UCP, state::UCPState, cstate::UCPState, σtr::Vec6)
    maxits = 50
    tol    = mat.ft * 1e-4
    h      = state.h

    ϵ = 1e-6

    De = calcDe(mat.E, mat.ν, state.ctx.stress_state)

    # Initial guess
    Δλ = eps()

    # Use the same initial direction you had (from cstate)
    ∂g∂σ = potential_derivs(mat, h, cstate.σ, cstate.εtp, cstate.εcp, cstate.εvp)

    σ   = σtr - Δλ * (De * ∂g∂σ)
    εcp = cstate.εcp
    εtp = cstate.εtp
    εvp = cstate.εvp

    f = yield_func(mat, h, σ, εtp, εcp, εvp)
    ω = 1.0

    # NR iterations
    for i in 1:maxits
        # Update flow direction at current iterate
        ∂g∂σ = potential_derivs(mat, h, σ, εtp, εcp, εvp)

        # Spectral split scalars (same as your original code)
        Λ  = eigvals(∂g∂σ)
        Λ1, Λ2, Λ3 = Λ

        max_Λp  = max(Λ1, Λ2, Λ3, 0.0)
        norm_Λn = (min(Λ1, 0.0)^2 + min(Λ2, 0.0)^2 + min(Λ3, 0.0)^2)^0.5
        sum_Λn  = abs(min(Λ1, 0.0) + min(Λ2, 0.0) + min(Λ3, 0.0))
        # norm_Λn = (sneg(Λ1,ϵ)^2 + sneg(Λ2,ϵ)^2 + sneg(Λ3,ϵ)^2)^0.5
        # sum_Λn  = abs(sneg(Λ1,ϵ) + sneg(Λ2,ϵ) + sneg(Λ3,ϵ))

        # f(Δλ) evaluator along the *current* direction (frozen ∂g∂σ and split scalars)
        function eval_f(Δλtest::Float64)
            σt   = σtr - Δλtest * (De * ∂g∂σ)
            εtpt = cstate.εtp + Δλtest * max_Λp
            εcpt = cstate.εcp + Δλtest * norm_Λn
            εvpt = cstate.εvp + Δλtest * sum_Λn
            return yield_func(mat, h, σt, εtpt, εcpt, εvpt)
        end

        # Numerical derivative ∂f/∂Δλ (central difference when possible)
        # Step scaled to magnitude of Δλ (robust default)
        fd_relstep = 1e-8
        fd_minstep = 1e-14
        δ = max(fd_minstep, fd_relstep * max(1.0, abs(Δλ)))

        # If you're enforcing Δλ ≥ 0, switch to forward diff near 0
        # (your current code allows negative during iterations, so central is usually ok)
        if (Δλ - δ) < 0.0
            fp = eval_f(Δλ + δ)
            ∂f∂Δλ = (fp - f) / δ
        else
            fp = eval_f(Δλ + δ)
            fm = eval_f(Δλ - δ)
            ∂f∂Δλ = (fp - fm) / (2.0 * δ)
        end

        # Safety
        if !isfinite(∂f∂Δλ) || abs(∂f∂Δλ) < eps(Float64)
            break
        end

        # Newton step direction
        Δλmin = max(Δλ - ω * f / ∂f∂Δλ, 0.0)
        fmin  = eval_f(Δλmin)

        # Backtracking line search
        for ω in 0.9:-0.1:0.3
            Δλtr = Δλ - ω * f / ∂f∂Δλ
        
            isfinite(Δλtr) || continue
            Δλtr > 0.0 || continue

            ftr   = eval_f(Δλtr)
            isfinite(ftr) || continue

            if abs(ftr) < abs(fmin)
                Δλmin = Δλtr
                fmin  = ftr
            end
        end

        Δλ = Δλmin
        f  = fmin
 
        isfinite(Δλ) || break

        # Update state variables using the same frozen scalars
        σ   = σtr - Δλ * (De * ∂g∂σ)
        εtp = cstate.εtp + Δλ * max_Λp
        εcp = cstate.εcp + Δλ * norm_Λn
        εvp = cstate.εvp + Δλ * sum_Λn

        f = yield_func(mat, h, σ, εtp, εcp, εvp)

        if abs(f) < tol
            Δλ < 0.0 && break

            w  = εtp * state.h
            ft = calc_ft(mat, w)
            fc = calc_fc(mat, εcp)
            abs(fc * mat.e / ft) > 1.1 || break

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


function calculate_apex_potential_slope(mat::UCP, ξa::Float64, ξc::Float64, κ::Float64)

    ϵ = 1e-6 * mat.ft # small offset
    ξ = ξc - ϵ
    
    # Meridian function (rξ) evaluation
    α     = mat.α
    drξdξ = -α / (ξc - ξa) * ((ξc - ξ) / (ξc - ξa))^(α - 1)

    # Potential derivative ∂g/∂ξ
    # g(σ) = ρ - rξ * κ since rχ = 1 in this domain
    dgdrξ = - κ  # rχ = 1 in this domain
    
    # mg = ∂g/∂ξ
    mg  = dgdrξ * drξdξ

    return mg
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
        Δσ      = state.σ - cstate.σ
        state.ε = cstate.ε + Δε
        return Δσ, success()
    end
    
    # Return to apex
    if tr(σtr) > 0.0
        # Calculate trial invariants
        ξ_tr  = tr(σtr) / √3
        s_tr  = dev(σtr)
        ρ_tr  = norm(s_tr)
        el_rt = (1 + mat.ν) / (1 - 2*mat.ν)
        
        ξa, ξc, κ = calc_ξa_ξc_κ(mat, h, cstate.εtp, cstate.εcp, cstate.εvp)
        mg = calculate_apex_potential_slope(mat, ξa, ξc, κ)

        if (ξ_tr - ξc) >= el_rt * ρ_tr * mg 
            state.σ = ξc * √3/3 * I2  # Pure hydrostatic tension
            Δσ      = state.σ - cstate.σ
            state.ε = cstate.ε + Δε

            # Calculate plastic strain geometrically
            Δε_p      = Δε - inv(De) \ Δσ

            state.εtp = cstate.εtp + tr(Δε_p) # Update tensile scalar
            G = mat.E / (2*(1+mat.ν))
            state.Δλ  = ρ_tr/(2*G) # Equivalent multiplier
            
            return Δσ, success()
        end
    end

    # Plastic update
    status = plastic_update(mat, state, cstate, σtr)
    failed(status) && return state.σ, status

    Δσ = state.σ - cstate.σ
    
    # Update Δεzz for plane stress (since update_state in uncoupled with Δεzz for plane stress)
    if state.ctx.stress_state == :plane_stress
        ∂g∂σ = potential_derivs(mat, h, state.σ, state.εtp, state.εcp, state.εvp)
        
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
    θ  = calc_θ(mat, σ)
    # r  = calc_rθ(mat, σ)

    w  = state.εtp*state.h
    ft = calc_ft(mat, w)
    fc = calc_fc(mat, state.εcp)

    ξa, ξc, κ = calc_ξa_ξc_κ(mat, h, state.εtp, state.εcp, state.εvp)
    # rχ = calc_rχ(mat, ξa, ξ)
    # rξ = calc_rξ(mat, ξa, ξc, ξ)

    vals_d = stress_strain_dict(σ, ε, state.ctx.stress_state)

    vals_d[:εcp] = state.εcp
    vals_d[:εtp] = state.εtp
    vals_d[:ξ]   = ξ
    vals_d[:ρ]   = ρ
    vals_d[:θ]   = θ
    vals_d[:fc]  = fc
    vals_d[:ft]  = ft
    vals_d[:ξa]  = ξa
    vals_d[:ξc]  = ξc
    vals_d[:κ]   = κ
    # vals_d[:r]   = r
    # vals_d[:rξ]  = rξ
    # vals_d[:rχ]  = rχ
    # vals_d[:ξb]  = 2*mat.fb/√3
    # vals_d[:fcb] = abs(mat.fc)

    return vals_d
end
