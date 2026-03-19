# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export ECP

"""
    ECP(; E, nu, fc, epsc, ft, GF=NaN, wc=NaN, beta=1.15, chi=0.2, p0=NaN,
         ft_law=:hordijk, fc_law=:popovics, H=0.0)

Evolving Concrete Plasticity model.

`ECP` is an invariant-based plasticity model for concrete with:
- a tensile-compressive yield surface written in terms of `ξ`, `ρ`, and `θ`
- a Willam-Warnke-type deviatoric section controlled by the eccentricity `e`
- a rounded compressive cap controlled by the compressive limit `ξc`
- tensile softening driven by crack opening `w ≈ h·εtp`
- compressive hardening/softening driven by the plastic compressive strain `εcp`
- isotropic compression hardening driven by the plastic volumetric strain `εvp`

The model is calibrated from the uniaxial compressive strength `fc`, tensile strength
`ft`, and biaxial strength ratio `beta`, with `fb = beta*fc` and `e = sqrt(beta/2)`.

# Keyword arguments
- `E::Real`:
  Young's modulus. Must be positive.
- `nu::Real`:
  Poisson's ratio. Must satisfy `0 <= nu < 0.5`.
- `fc::Real`:
  Uniaxial compressive strength. Must be negative.
- `epsc::Real`:
  Strain at the compressive peak. Must be negative.
- `ft::Real`:
  Uniaxial tensile strength. Must be positive.
- `GF::Real = NaN`:
  Tensile fracture energy. Used to infer `wc` when `wc` is not provided.
- `wc::Real = NaN`:
  Critical crack opening. If omitted, it is estimated from `GF`.
- `beta::Real = 1.15`:
  Ratio between biaxial and uniaxial compressive strengths, `fb/fc`.
- `chi::Real = 0.2`:
  Dilatance ratio `χ = tan(ψ)` used in the plastic potential.
- `p0::Real = NaN`:
  Initial elastic limit in isotropic compression. If omitted, `ξc0` is set
  automatically from the biaxial strength level.
- `ft_law = :hordijk`:
  Tensile softening law. A spline curve can be provided.
- `fc_law = :popovics`:
  Compressive evolution law. A spline curve can be provided.
- `H::Real = 0.0`:
  Hardening modulus for isotropic compression.

# Notes
- The model is designed for matching uniaxial and biaxial compressive strengths as well as the tensile strength.
- The tensile limit `ξt` evolves with the current strengths `ft` and `fc`.
- The tensile response is regularized by `wc`/`GF` and the element characteristic length `h`.
- The initial compressive limit is given by `ξc0=2·ξb`, with `ξb = 2fb/√3`.
- The compressive cap is active for `ξc < ξ < 1.2·ξb`.
- The model is intended for `MechSolid` elements.
"""
mutable struct ECP<:Constitutive
    E::Float64
    ν::Float64
    fc::Float64
    εc::Float64
    ft::Float64
    wc::Float64
    fb::Float64
    χ::Float64
    H::Float64
    ξc0::Float64
    ξt0::Float64
    ft_law::Symbol
    ft_fun::Union{Nothing,AbstractSpline}
    fc_law::Symbol
    fc_fun::Union{Nothing,AbstractSpline}
    e::Float64

    function ECP(;
        E::Real    = NaN,
        nu::Real   = NaN,
        fc::Real   = NaN,
        epsc::Real = NaN,
        ft::Real   = NaN,
        GF::Real   = NaN,
        wc::Real   = NaN,
        beta::Real = 1.15,
        chi::Real  = 0.2,
        p0::Real   = NaN,
        ft_law     = :hordijk,
        fc_law     = :popovics,
        H::Real    = 0.0,
    )
        @check E>0 "ECP: Young's modulus E must be > 0. Got $E."
        @check 0<=nu<0.5 "ECP: Poisson's ratio nu must be in the range [0, 0.5). Got $nu."
        @check 0.1<chi<=1.0 "ECP: Dilatance ratio χ=tan(ψ) [0.1, 1.0]. Got $chi."
        @check 1<=beta<=1.5 "ECP: Factor beta must be in the range [1.0, 1.5]. Got $beta."

        @check ft>0 "ECP: Tensile strength ft must be > 0. Got $ft."
        @check H>=0 "ECP: Plastic modulus H must be >= 0. Got $H."

        wc, ft_law, ft_fun, status = setup_tensile_strength(ft, GF, wc, ft_law)
        failed(status) && throw(ArgumentError("ECP: " * status.message))

        fc_law, fc_fun, status = setup_compressive_strength(E, fc, epsc, fc_law)
        failed(status) && throw(ArgumentError("ECP: " * status.message))

        @check abs(epsc)>abs(fc)/E "ECP: epsc should be greater than fc/E."

        # Excentricity matching the biaxial compressive strength.
        β  = beta
        e  = √(β/2)
        fb = β*fc
        
        if isnan(wc)
            @check GF>0 "ECP: Fracture energy GF must be > 0. Got $(repr(GF))."
            wc = round(GF/(0.1947*ft), sigdigits=5)  # inverse of Hordijk approximation
            notify("ECP: Using Hordijk's approximation wc=$(repr(wc)).")
        else
            @check wc>=0 "ECP: Critical crack opening wc must be >= 0. Got $(repr(wc))."
        end
        wc > 1e-5 || notify("ECP: Warning: very low value of wc=$(repr(wc)).")

        if isnan(p0)
            ξb  = (2/√3*fb)*1.2 # cap position (ξb) 20% beyond the biaxial strength (2/√3*fb)
            ξc0 = 1.5*ξb
        else
            @check p0<0 "ECP: Elastic limit in isotropic compression p0 must be < 0. Got $(repr(p0))."
            ξc0 = √3*p0
        end

        den = e^2*fc^2 - ft^2
        @check den > 0 "ECP: invalid strength ratio. Expected e^2*fc^2 > ft^2."
        ξt0 = fc*ft*(e^2*fc - ft)/(√3*den)

        return new(E, nu, fc, epsc, ft, wc, fb, chi, H, ξc0, ξt0, ft_law, ft_fun, fc_law, fc_fun, e)
    end
end


mutable struct ECPState<:ConstState
    ctx::Context
    σ  ::Vec6
    ε  ::Vec6
    εtp::Float64
    εcp::Float64
    εvp::Float64
    Δλ ::Float64
    h  ::Float64
    function ECPState(ctx::Context)
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
# compat_state_type(::Type{ECP}, ::Type{MechSolid}) = ctx.stress_state!=:plane_stress ? ECPState : error("ECP: This model is not compatible with planestress")
compat_state_type(::Type{ECP}, ::Type{MechSolid}) = ECPState


function calc_θ(mat::ECP, ρ::Float64, j3::Float64)

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


function calc_rθ(mat::ECP, θ::Float64)
    e = mat.e

    rnum   = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden   = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    r      = rnum/rden

    return r
end


function calc_rξ(mat::ECP, ξc::Float64, ξt::Float64, ξ::Float64)
    abs_fc = abs(mat.fc)
    return spow((ξt-ξ)/abs_fc, mat.α)
end


function calc_rc(mat::ECP, ξc::Float64, ξ::Float64)
    ξb = 2*mat.fb/√3
    ξ>=ξb && return 1.0
    ξ<ξc  && return 0.0
    return √(1 - ((ξb-ξ)/(ξb-ξc))^2)
end


function is_apex_state(mat::ECP, ξ::Float64, ξt::Float64)
    ξtol = √eps(Float64)*max(abs(mat.fc), abs(ξt), 1.0)
    return ξ >= ξt - ξtol
end


function calc_fc(mat::ECP, εcp::Float64)
    fc0 = 0.4*mat.fc
    fcr = 0.1*mat.fc

    return calc_compressive_strength(mat, fc0, fcr, εcp)
end

function calc_fc_derivative(mat::ECP, εcp::Float64)
    fc0 = 0.4*mat.fc
    fcr = 0.1*mat.fc

    return calc_compressive_strength_derivative(mat, fc0, fcr, εcp)
end


function calc_ft(mat::ECP, w::Float64)
    return calc_tensile_strength(mat, w)
end


function calc_ft_derivative(mat::ECP, w::Float64)
    ∂ft∂w = calc_tensile_strength_derivative(mat, w)
    Hcap  = isfinite(mat.wc) ? -mat.ft/(0.5*mat.wc) : -Inf
    return max(∂ft∂w, Hcap)
end


function calc_ξc_ξt_m(mat::ECP, h::Float64, εtp::Float64, εcp::Float64, εvp::Float64)
    w  = εtp*h
    e  = mat.e
    ft = calc_ft(mat, w)
    fc = calc_fc(mat, εcp)

    den = e^2*fc^2 - ft^2
    ξt  = fc*ft*(e^2*fc - ft)/(√3*den)
    
    # p = p0 + H*εvp  -> ξc = √3*p0 + √3*H*εvp
    ξc = mat.ξc0 + √3*mat.H*εvp # cap hardening

    @assert ξc<0
    @assert ξc<fc/√3

    m = (2/3)*fc^2/(ξt - fc/√3)
    @assert m>0

    return ξc, ξt, m
end


function yield_func(mat::ECP, h::Float64, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    # f(σ) = ρ² - m⋅rθ²⋅rc²⋅(ξt-ξ)

    i1, j2 = tr(σ), J2(σ)

    ξ = i1/√3
    ρ = √(2*j2)

    ξc, ξt, m = calc_ξc_ξt_m(mat, h, εtp, εcp, εvp)
    θ = calc_θ(mat, ρ, J3(σ))
    rθ = calc_rθ(mat, θ)
    rc = calc_rc(mat, ξc, ξ)

    return ρ^2 - m*rθ^2*rc^2*(ξt - ξ)
end


function yield_derivs(mat::ECP, h::Float64, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    i1, j2, j3 = tr(σ), J2(σ), J3(σ)
    ξc, ξt, m = calc_ξc_ξt_m(mat, h, εtp, εcp, εvp)
    
    ξb = 2*mat.fb/√3
    ξ  = i1/√3
    ρ  = √(2*j2)
    θ  = calc_θ(mat, ρ, j3)
    rθ = calc_rθ(mat, θ)
    rc = calc_rc(mat, ξc, ξ)

    # ∂f/∂εtp
    fc = calc_fc(mat, εcp)
    A  = ξt - fc/√3
    ∂f∂m  = -rθ^2*rc^2*(ξt - ξ)
    ∂m∂ξt = -(2/3)*fc^2/A^2
    ∂f∂ξt = -m*rθ^2*rc^2 + ∂f∂m*∂m∂ξt

    dξtdft = mat.ξt0/mat.ft
    w      = εtp*h
    ∂ft∂w  = calc_ft_derivative(mat, w)
    ∂w∂εtp = h
    ∂f∂εtp = ∂f∂ξt*dξtdft*∂ft∂w*∂w∂εtp

    # ∂f/∂εcp
    ∂m∂fc   = (2/3)*fc*(2*ξt - fc/√3)/A^2
    ∂fc∂εcp = calc_fc_derivative(mat, εcp)
    ∂f∂εcp  = ∂f∂m*∂m∂fc*∂fc∂εcp

    # ∂f/∂εvp
    if mat.H!=0.0
        ∂rc∂ξc  = ξc < ξ < ξb ? -(ξb-ξ)^2/(ξb-ξc)^3/√(1 - ((ξb-ξ)/(ξb-ξc))^2) : 0.0
        ∂f∂rc   = -2*m*rθ^2*rc*(ξt - ξ)
        ∂f∂ξc   = ∂f∂rc*∂rc∂ξc
        ∂ξc∂εvp = √3*mat.H
        ∂f∂εvp  = ∂f∂ξc*∂ξc∂εvp
    else
        ∂f∂εvp = 0.0
    end

    # f derivative w.r.t. σ:
    ∂f∂ρ  = 2*ρ
    ∂f∂ξ  = m*rθ^2*rc^2
    ∂rc∂ξ = ξc < ξ < ξb ? (ξb-ξ)/(ξb-ξc)^2/√(1 - ((ξb-ξ)/(ξb-ξc))^2) : 0.0
    ∂f∂ξ += -2*m*rθ^2*rc*(ξt - ξ)*∂rc∂ξ
    dξdσ = √3/3*I2

    use_lode_derivative = !( θ == 0.0 || θ == π/3 ) # avoid singularity at meridians (θ=0, π/3) and apex (ρ=0)

    if use_lode_derivative
        s = dev(σ)
        e = mat.e
        rnum = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
        rden = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2

        ∂rθ∂numdθ = (2*sin(2*θ)*(2*e-1)*(e^2-1))/√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e) - 2*(1 - e^2)*sin(θ)
        ∂rθ∂dendθ = 4*sin(2*θ)*(e^2-1)
        ∂rθ∂θ    = (∂rθ∂numdθ*rden - rnum*∂rθ∂dendθ)/rden^2
        ∂f∂θ     = -2*m*rθ*rc^2*(ξt - ξ)*∂rθ∂θ
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


function potential_derivs(mat::ECP, h::Float64, σ::AbstractArray, εtp::Float64)
    # g(σ) = ρ^2 - 4 tan^2(ψ) (ξ_t - ξ_f'_c) (ξ_t - ξ) = 0
    # g(σ) = ρ^2 - 4 χ^2 (ξ_t - ξ_f'_c) (ξ_t - ξ) = 0
    ξfc = mat.fc/√3

    w  = εtp*h
    ft = calc_ft(mat, w)
    ξt = mat.ξt0*ft/mat.ft
    ρ  = √(2*J2(σ))

    dgdξ = 4*mat.χ^2*(ξt - ξfc)
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


function ecp_plastic_flow_invariant_rates(∂g∂σ::Vec6)
    # Recover the principal values analytically from invariants, then apply
    # the original positive/negative spectral split.
    Λ1, Λ2, Λ3 = eigvals(∂g∂σ)

    # rate_εtp = max(Λ1, Λ2, Λ3, 0.0)
    rate_εtp = sqrt(max(Λ1, 0.0)^2 + max(Λ2, 0.0)^2 + max(Λ3, 0.0)^2) # what about p-norm with p around 5
    rate_εcp = sqrt(min(Λ1, 0.0)^2 + min(Λ2, 0.0)^2 + min(Λ3, 0.0)^2)
    rate_εvp = abs(min(Λ1, 0.0) + min(Λ2, 0.0) + min(Λ3, 0.0))

    return rate_εtp, rate_εcp, rate_εvp
end


function calcD(mat::ECP, state::ECPState)
    De  = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    h = state.h

    state.Δλ==0.0 && return De

    ∂f∂σ, ∂f∂εtp, ∂f∂εcp, ∂f∂εvp = yield_derivs(mat, h, state.σ, state.εtp, state.εcp, state.εvp)
    ∂g∂σ = potential_derivs(mat, h, state.σ, state.εtp)
    rate_εtp, rate_εcp, rate_εvp = ecp_plastic_flow_invariant_rates(∂g∂σ)

    De_dgdσ = De*∂g∂σ
    denom = ∂f∂σ'*De_dgdσ - ∂f∂εcp*rate_εcp - ∂f∂εtp*rate_εtp - ∂f∂εvp*rate_εvp
    Dep = De - De_dgdσ*∂f∂σ'*De / denom

    return Dep
end


function plastic_update(mat::ECP, state::ECPState, cstate::ECPState, σtr::Vec6)
    maxits = 50
    tol    = mat.ft^2*1e-5
    h      = state.h
    ∂g∂σ   = potential_derivs(mat, h, cstate.σ, cstate.εtp)
    De     = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    Δλ     = eps()

    σ  = σtr - Δλ*(De*∂g∂σ)

    εcp = cstate.εcp
    εtp = cstate.εtp
    εvp = cstate.εvp

    f = yield_func(mat, h, σ, εtp, εcp, εvp)
    rate_εtp, rate_εcp, rate_εvp = ecp_plastic_flow_invariant_rates(∂g∂σ)

    # NR iterations considering ∂g∂σ frozen
    for i in 1:maxits
        ∂f∂σ, ∂f∂εtp, ∂f∂εcp, ∂f∂εvp = yield_derivs(mat, h, σ, εtp, εcp, εvp)

        ∂f∂Δλ = -∂f∂σ'*De*∂g∂σ + ∂f∂εcp*rate_εcp + ∂f∂εtp*rate_εtp + ∂f∂εvp*rate_εvp

        # Newton step direction
        Δλ = max(Δλ - f / ∂f∂Δλ, 0.0)
        isfinite(Δλ) || break

        σ   = σtr - Δλ*(De*∂g∂σ)
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

            state.σ   = σ
            state.εtp = εtp
            state.εcp = εcp
            state.εvp = εvp
            state.Δλ  = Δλ
            
            return success()
        end

    end

    return failure("ECP: plastic update failed")
end


function update_state(mat::ECP, state::ECPState, cstate::ECPState, Δε::AbstractArray)

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


function state_values(mat::ECP, state::ECPState)
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
