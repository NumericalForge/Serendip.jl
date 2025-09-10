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
  Curvature coefficient of the meridional section (0.2 < α ≤ 1.0).
- `beta::Real = 1.15`:  
  Factor relating biaxial to uniaxial compressive strength (1 ≤ β ≤ 1.5).
- `H::Real = 0.0`:  
  Plastic modulus for isotropic compression (≥ 0).

# Returns
A `UCP` material object that can be attached to mechanical bulk elements
for 2D (plane strain) or 3D analyses. Not compatible with plane stress.

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
    p0::Float64
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

        fc_law, fc_fun, status = setup_compressive_strength(fc, epsc, fc_law)
        failed(status) && throw(ArgumentError("UCP: " * status.message))

        fc_fun = nothing
        if fc_law isa AbstractSpline
            fc_fun = fc_law
            fc_law = :custom
            fc     = fc_law(0.0)
        end

        α = alpha
        β = beta
        
        # value of exentricity to match fb in a biaxial trajectory, assuming the state when ξb=0
        e  = β/(2*β)^α
        fb = β*fc
        
        if isnan(wc)
            @check GF>0 "UCP: Fracture energy GF must be > 0. Got $(repr(GF))."
            wc = round(GF/(0.1947*ft), sigdigits=5)  # inverse of Hordijk approximation
            notify("UCP: Using Hordijk's approximation wc=$(repr(wc)).")
        else
            @check wc>=0 "UCP: Critical crack opening wc must be >= 0. Got $(repr(wc))."
        end

        if isnan(p0)
            ξc = 2*fb/√3
            ξa = 1.5*ξc
            p0 = ξa/√3
        else
            @check p0<0 "UCP: Elastic limit in isotropic compression p0 must be < 0. Got $(repr(p0))."
        end
        @assert p0<0

        return new(E, nu, fc, epsc, eta, ft, wc, fb, p0, ft_law, ft_fun, fc_law, fc_fun, α, e, H)
    end
end


mutable struct UCPState<:IpState
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
compat_state_type(::Type{UCP}, ::Type{MechBulk}, ctx::Context) = ctx.stress_state!=:plane_stress ? UCPState : error("UCP: This model is not compatible with planestress")


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


function calc_rξ(mat::UCP, ξb::Float64, ξ::Float64)
    α  = mat.α
    fc_peak = abs(mat.fc)

    return spow((ξb-ξ)/fc_peak, α)
end

function calc_rc(mat::UCP, ξa::Float64, ξ::Float64)
    ξc = 2*mat.fb/√3
    ξ>=ξc && return 1.0
    ξ<ξa  && return 0.0
    return √(1 - ((ξc-ξ)/(ξc-ξa))^2)
end


function calc_fc(mat::UCP, εcp::Float64)
    fc0 = 0.35*mat.fc
    fcr = 0.1*mat.fc
    return calc_compressive_strength(mat, fc0, fcr, εcp)
end


function calc_ft(mat::UCP, w::Float64)
    return calc_tensile_strength(mat, w)
end


function calc_p(mat::UCP, εvp::Float64)
    return mat.p0 + mat.H*εvp
end


function calc_ξa_ξb_κ(mat::UCP, state::UCPState, εtp::Float64, εcp::Float64, εvp::Float64)
    e  = mat.e
    α  = mat.α
    w  = εtp*state.h

    ft      = calc_ft(mat, w)
    fc      = calc_fc(mat, εcp)
    fc_peak = abs(mat.fc)

    # ξa = √3*mat.p_fun(√3*εcp)    # ξ = √3p ; plastic volumetric strain εvp = √3*εcp in isotropic compression
    p = calc_p(mat, εvp)
    ξa = √3*p    # ξ = √3p
    # ξa = √3*mat.p_fun(εvp)    # ξ = √3p
    @assert ξa<0
    @assert ξa<fc/√3

    Ω  = (-ft/(fc*e))^(1/α)
    ξb = 1/√3*(fc*Ω - ft)/(Ω-1)
    ξb<0 && @show ξb

    if ξb<0
        @show α
        @show Ω
        @show εcp
        @show εtp
        @show fc
        @show ft
        @show ξa
        @show ξb
    end

    κ  = -√(2/3)*fc*((ξb-fc/√3)/fc_peak)^-α
    @assert κ>0

    return ξa, ξb, κ
end


function yield_func(mat::UCP, state::UCPState, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    # f(σ) = ρ - rθ⋅rc⋅rξ⋅κ

    i1, j2 = tr(σ), J2(σ)

    ξ = i1/√3
    ρ = √(2*j2)

    ξa, ξb, κ = calc_ξa_ξb_κ(mat, state, εtp, εcp, εvp)
    rθ = calc_rθ(mat, σ)
    rc = calc_rc(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξb, ξ)

    return ρ - rθ*rc*rξ*κ
end


function yield_derivs(mat::UCP, state::UCPState, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    e = mat.e
    α = mat.α
    fc_peak = abs(mat.fc)

    i1, j2 = tr(σ), J2(σ)

    ρ = √(2*j2)
    ξ = i1/√3

    # deviatoric derivatives
    s      = dev(σ)
    det_s  = J3(σ)
    adj_s  = det_s*inv(s)
    norm_s = ρ

    # θ and derivatives
    θ        = 1/3*acos( clamp(3*√6*det_s/norm_s^3, -1.0, 1.0) )
    rnum     = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden     = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    rθ       = rnum/rden
    drθnumdθ = (2*sin(2*θ)*(2*e-1)*(e^2-1))/√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e) - 2*(1 - e^2)*sin(θ)
    drθdendθ = 4*sin(2*θ)*(e^2-1)
    drθdθ    = (drθnumdθ*rden - rnum*drθdendθ)/rden^2

    if 1-abs(cos(3*θ)) > 1e-6 # condition to avoid division by zero
        dθds = -√6*(adj_s/ρ^3 - 3*s*det_s/ρ^5)/√abs(1 - 54*det_s^2/ρ^6)
    else
        dθds = 0.0*I2
    end

    ξa, ξb, κ = calc_ξa_ξb_κ(mat, state, εtp, εcp, εvp)

    ξc = 2*mat.fb/√3
    rc = calc_rc(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξb, ξ)

    # f derivative w.r.t. σ:
    dfdρ  = 1.0
    dfdrc = -rθ*rξ*κ
    dfdrξ = -rθ*rc*κ
    drcdξ = ξa<ξ<ξc ? (ξc-ξ)/(ξc-ξa)^2/√(1-((ξc-ξ)/(ξc-ξa))^2) : 0.0
    drξdξ = ξb-ξ!=0.0 ? -α/fc_peak * abs((ξb-ξ)/fc_peak)^(α-1) : 0.0
    dfdξ  = dfdrc*drcdξ + dfdrξ*drξdξ
    dfdrθ = -rc*rξ*κ
    dfdθ  = dfdrθ*drθdθ

    dρdσ = s/norm(s)
    dξdσ = √3/3*I2
    dsdσ = Psd
    dθdσ = dsdσ*dθds

    if ρ==0 # apex
        dfdσ = √3/3*I2
    else
        dfdσ = dfdρ*dρdσ + dfdξ*dξdσ + dfdθ*dθdσ
    end

    f_εcp  = εcp -> yield_func(mat, state, σ, εtp, εcp, εvp)
    dfdεcp = derive(f_εcp, εcp)

    f_εtp  = εtp -> yield_func(mat, state, σ, εtp, εcp, εvp)
    dfdεtp = derive(f_εtp, εtp)

    return dfdσ, dfdεtp, dfdεcp
end


function potential_derivs(mat::UCP, state::UCPState, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    # f(σ) = ρ - e⋅rc⋅rξ⋅κ

    e  = mat.e
    α  = mat.α
    fc_peak = abs(mat.fc)

    i1 = tr(σ)
    ξ  = i1/√3
    s  = dev(σ)
    ρ  = norm(s)
    ρ == 0 && return  √3/3*I2

    ξa, ξb, κ = calc_ξa_ξb_κ(mat, state, εtp, εcp, εvp)

    ξc = 2*mat.fb/√3
    rc = calc_rc(mat, ξa, ξ)
    rξ = calc_rξ(mat, ξb, ξ)

    dgdrc = -e*rξ*κ
    dgdrξ = -e*rc*κ
    drcdξ = ξa<ξ<ξc ? (ξc-ξ)/(ξc-ξa)^2/√(1-((ξc-ξ)/(ξc-ξa))^2) : 0.0
    drξdξ = ξb-ξ!=0.0 ? -α/fc_peak * abs((ξb-ξ)/fc_peak)^(α-1) : 0.0
    dgdξ  = dgdrc*drcdξ + dgdrξ*drξdξ

    dξdσ = √3/3*I2
    dgdρ = 1.0

    # if ξ>ξb && (ξ-ξb)*dgdξ > ρ*dgdρ  # apex
    #     dgdσ = ξb*I2 - σ
    # else
    #     dgdσ = s/ρ + dgdξ*dξdσ
    # end
    dgdσ = s/ρ + dgdξ*dξdσ
    return dgdσ

end


function calcD(mat::UCP, state::UCPState)
    De  = calcDe(mat.E, mat.ν, state.ctx.stress_state)

    state.Δλ==0.0 && return De

    dfdσ, dfdεtp, dfdεcp = yield_derivs(mat, state, state.σ, state.εtp, state.εcp, state.εvp)
    dgdσ = potential_derivs(mat, state, state.σ, state.εtp, state.εcp, state.εvp)

    Λ = eigvals(dgdσ, sort=false)
    # Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεcp*norm(min.(0.0, Λ)) - dfdεtp*norm(max.(0.0, Λ)))
    Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεcp*norm(min.(0.0, Λ)) - dfdεtp*maximum(max.(0.0, Λ)))
    
    return Dep
end


function calc_σ_εp_Δλ(mat::UCP, state::UCPState, σtr::Vec6)
    maxits = 60
    tol    = 0.1
    tol    = 1.0
    dgdσ   = potential_derivs(mat, state, state.σ, state.εtp, state.εcp, state.εvp)
    De     = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    Δλ     = eps()

    σ  = σtr - Δλ*(De*dgdσ)

    εcp = state.εcp
    εtp = state.εtp
    εvp = state.εvp

    f   = yield_func(mat, state, state.σ, εtp, εcp, εvp)
    eta   = 1.0 # initial damping

    # iterative process
    for i in 1:maxits
        dfdσ, _ = yield_derivs(mat, state, σ, εtp, εcp, εvp)
        dgdσ    = potential_derivs(mat, state, σ, εtp, εcp, εvp)
        dfdΔλ   = -dfdσ'*De*dgdσ

        Δλ = Δλ - eta*f/dfdΔλ
        if Δλ<0
            # Δλ = abs(Δλ)
            # @show Δλ
        end

        if isnan(Δλ)
            return state.σ, 0.0, 0.0, 0.0, 0.0, failure("UCP: Δλ is NaN")
        end

        σ  = σtr - Δλ*(De*dgdσ)

        Λ   = eigvals(dgdσ, sort=false)
        εtp = state.εtp + Δλ*maximum(max.(0.0, Λ))
        εcp = state.εcp + Δλ*norm(min.(0.0, Λ))
        εvp = state.εvp + Δλ*sum(abs, min.(0.0, Λ))
        f   = yield_func(mat, state, σ, εtp, εcp, εvp)

        if abs(f) < tol
            Δλ < 0.0 && return σ, 0.0, 0.0, 0.0, 0.0, failure("UCP: negative Δλ")

            return σ, εtp, εcp, εvp, Δλ, success()
        end

        # dumping
        i>10 && (eta = 0.6)
        i>15 && (eta = 0.3)
    end

    return state.σ, 0.0, 0.0, 0.0, 0.0, failure("UCP: maximum iterations reached")
end


function update_state(mat::UCP, state::UCPState, Δε::AbstractArray)
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    σtr  = state.σ + De*Δε
    ftr  = yield_func(mat, state, σtr, state.εtp, state.εcp, state.εvp)

    Δλ  = 0.0
    tol = 1.0
    tol = 0.1

    if ftr < tol
        # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # plastic
        state.σ, state.εtp, state.εcp, state.εvp, state.Δλ, status = calc_σ_εp_Δλ(mat, state, σtr)
        @assert state.εcp >= 0.0
        @assert state.εtp >= 0.0
        @assert state.εvp >= 0.0

        Δσ = state.σ - σini

        failed(status) && return state.σ, status
    end

    state.ε += Δε
    Δσ       = state.σ - σini
    return Δσ, success()
end


function state_values(mat::UCP, state::UCPState)
    σ, ε  = state.σ, state.ε
    ρ = √(2*J2(σ))
    ξ = tr(σ)/√3
    # θ  = calc_θ(mat, σ)
    # r  = calc_rθ(mat, σ)

    w  = state.εtp*state.h
    ft = calc_ft(mat, w)
    fc = calc_fc(mat, state.εcp)
    # ft = mat.ft_fun(w)
    # fc = mat.fc_fun(state.εcp)

    # ξa, ξb, κ = calc_ξa_ξb_κ(mat, state, state.εtp, state.εcp, state.εvp)
    # rc = calc_rc(mat, ξa, ξ)
    # rξ = calc_rξ(mat, ξb, ξ)

    vals_d = stress_strain_dict(σ, ε, state.ctx.stress_state)

    vals_d[:εcp] = state.εcp
    vals_d[:εtp] = state.εtp
    vals_d[:ξ]   = ξ
    vals_d[:ρ]   = ρ
    # vals_d[:θ]   = θ
    vals_d[:fc]  = fc
    vals_d[:ft]  = ft
    # vals_d[:ξa]  = ξa
    # vals_d[:ξb]  = ξb
    # vals_d[:κ]   = κ
    # vals_d[:r]   = r
    # vals_d[:rξ]  = rξ
    # vals_d[:rc]  = rc
    # vals_d[:ξc]  = 2*mat.fb/√3
    # vals_d[:fcb] = abs(mat.fc)

    return vals_d
end
