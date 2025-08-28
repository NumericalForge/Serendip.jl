# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export CSCP

CSCP_params = [
    FunInfo(:CSCP, "Closed surface concrete plasticity"),
    KwArgInfo(:E, "Young's modulus", cond=:(E>0)),
    KwArgInfo(:nu, "Poisson's ratio", cond=:(0<=nu<0.5)),
    KwArgInfo(:alpha, "Curvature coefficient", 0.666, cond=:(0.2<alpha<=1.0)),
    KwArgInfo(:beta, "Factor to get fb from fc_max", 1.15, cond=:(1<=beta<=1.5)),
    KwArgInfo(:fc, "Compressive strength", 0.0, cond=:(fc<0)),
    KwArgInfo(:epsc, "Strain corresponding to the peak value in the uniaxial compression curve", 0.0, cond=:(epsc<0)),
    KwArgInfo(:n, "Shape curve parameter for the compression curve", 4.0, cond=:(n>1)),
    KwArgInfo(:ft, "Tensile strength", 0.0, cond=:(ft>0)),
    KwArgInfo(:wc, "Critical crack opening in uniaxial extension", 0.0, cond=:(wc>=0)),
    KwArgInfo(:GF, "Tensile fracture energy", 0.0, cond=:(GF>0)),
    KwArgInfo(:p0, "Elastic limit for p under isotropic compression", nothing),
    KwArgInfo(:Hp, "Plastic modulus for the isotropic compression", 0.0, cond=:(Hp>=0)),
]

mutable struct CSCP<:Constitutive
    E::Float64
    ν::Float64
    e::Float64
    α::Float64
    ft::Float64
    fc::Float64
    εc::Float64
    n::Float64
    GF::Float64
    wc::Float64
    fb::Float64
    p0::Float64
    Hp::Float64

    function CSCP(; kwargs...)
        args = checkargs(kwargs, CSCP_params)

        α = args.alpha
        β = args.beta

        # value of exentricity to match fb in a biaxial trajectory, assuming the state when ξb=0
        e = β/(2*β)^α
        fb = β*args.fc

        if args.p0 === nothing
            ξc = 2*fb/√3
            ξa = 1.5*ξc
            p0 = ξa/√3
        else
            p0 = args.p0
            @assert p0<0
        end

        this = new(args.E, args.nu, e, α, args.ft, args.fc, args.epsc, args.n, args.GF, args.wc, fb, p0, args.Hp)
        return this
    end
end


mutable struct CSCPState<:IpState
    ctx::Context
    σ  ::Vec6
    ε  ::Vec6
    εtp::Float64
    εcp::Float64
    εvp::Float64
    Δλ ::Float64
    h  ::Float64
    function CSCPState(ctx::Context)
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
compat_state_type(::Type{CSCP}, ::Type{MechBulk}, ctx::Context) = ctx.stress_state!=:plane_stress ? CSCPState : error("CSCP: This model is not compatible with planestress")


function calc_θ(::CSCP, σ::Vec6)
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


function calc_rθ(mat::CSCP, σ::Vec6)
    e = mat.e
    θ = calc_θ(mat, σ)

    rnum   = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden   = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    r      = rnum/rden

    return r
end


function calc_rξ(mat::CSCP, ξb::Float64, ξ::Float64)
    α  = mat.α
    fc_abs = abs(mat.fc)

    return spow((ξb-ξ)/fc_abs, α)
end

function calc_rc(mat::CSCP, ξa::Float64, ξ::Float64)
    ξc = 2*mat.fb/√3
    ξ>=ξc && return 1.0
    ξ<ξa  && return 0.0
    return √(1 - ((ξc-ξ)/(ξc-ξa))^2)
end

function calc_fc(mat::CSCP, εcp::Float64)
    fc     = mat.fc
    εcp_pk = abs(mat.εc) - abs(fc)/mat.E  # εcp_pk is the compression plastic strain at the peak of the uniaxial compression curve
    χ      = εcp/εcp_pk
    n      = mat.n
    fc0    = 0.4*fc
    fcr    = 0.5*fc

    if εcp < 0.0
        fc_cur = fc0
    elseif εcp < εcp_pk
        # before peak
        fc_cur = fc0 + (fc - fc0) * n * χ/(n - 1 + χ^n)
    else
        # after peak
        fc_cur = fcr + (fc - fcr) * n * χ/(n - 1 + χ^n)
    end

    return fc_cur
end


function calc_ft(mat::CSCP, w::Float64)
    wc = mat.wc
    if w < 0
        z = 1.0
    elseif w < wc
        z = (1 + 27*(w/wc)^3)*exp(-6.93*w/wc) - 28*(w/wc)*exp(-6.93)
    else
        z = 0.0
    end
    return z*mat.ft
end

function calc_p(mat::CSCP, εvp::Float64)
    return mat.p0 + mat.Hp*εvp
end


function calc_ξa_ξb_κ(mat::CSCP, state::CSCPState, εtp::Float64, εcp::Float64, εvp::Float64)

    α  = mat.α
    e  = mat.e
    w  = εtp*state.h

    ft_cur = calc_ft(mat, w)
    fc_cur = calc_fc(mat, εcp)
    fc_abs = abs(mat.fc)

    # ξa = √3*mat.p_fun(√3*εcp)    # ξ = √3p ; plastic volumetric strain εvp = √3*εcp in isotropic compression
    p = calc_p(mat, εvp)
    ξa = √3*p    # ξ = √3p
    # ξa = √3*mat.p_fun(εvp)    # ξ = √3p
    @assert ξa<0
    @assert ξa<fc_cur/√3

    Ω  = (-ft_cur/(fc_cur*e))^(1/α)
    ξb = 1/√3*(fc_cur*Ω - ft_cur)/(Ω-1)
    ξb<0 && @show ξb

    if ξb<0
        @show α
        @show Ω
        @show εcp
        @show εtp
        @show fc_cur
        @show ft_cur
        @show ξa
        @show ξb
    end

    κ  = -√(2/3)*fc_cur*((ξb-fc_cur/√3)/fc_abs)^-α
    @assert κ>0

    return return ξa, ξb, κ
end


function yield_func(mat::CSCP, state::CSCPState, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
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


function yield_derivs(mat::CSCP, state::CSCPState, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    e = mat.e
    α = mat.α
    fc_abs = abs(mat.fc)

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
    dfdρ  = 1
    dfdrc = -rθ*rξ*κ
    dfdrξ = -rθ*rc*κ
    drcdξ = ξa<ξ<ξc ? (ξc-ξ)/(ξc-ξa)^2/√(1-((ξc-ξ)/(ξc-ξa))^2) : 0.0
    drξdξ = ξb-ξ!=0.0 ? -α/fc_abs * abs((ξb-ξ)/fc_abs)^(α-1) : 0.0
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


function potential_derivs(mat::CSCP, state::CSCPState, σ::AbstractArray, εtp::Float64, εcp::Float64, εvp::Float64)
    # f(σ) = ρ - e⋅rc⋅rξ⋅κ

    e  = mat.e
    α  = mat.α
    fc_abs = abs(mat.fc)

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
    drξdξ = ξb-ξ!=0.0 ? -α/fc_abs * abs((ξb-ξ)/fc_abs)^(α-1) : 0.0
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


function calcD(mat::CSCP, state::CSCPState)
    De  = calcDe(mat.E, mat.ν, state.ctx.stress_state)

    state.Δλ==0.0 && return De

    dfdσ, dfdεtp, dfdεcp = yield_derivs(mat, state, state.σ, state.εtp, state.εcp, state.εvp)
    dgdσ = potential_derivs(mat, state, state.σ, state.εtp, state.εcp, state.εvp)

    Λ = eigvals(dgdσ, sort=false)
    # Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεcp*norm(min.(0.0, Λ)) - dfdεtp*norm(max.(0.0, Λ)))
    Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεcp*norm(min.(0.0, Λ)) - dfdεtp*maximum(max.(0.0, Λ)))
    return Dep
end


function calc_σ_εp_Δλ(mat::CSCP, state::CSCPState, σtr::Vec6)
    maxits = 60
    tol    = 0.1
    tol    = 1.0
    dgdσ   = potential_derivs(mat, state, state.σ, state.εtp, state.εcp, state.εvp)
    De     = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    # Δλ     = norm(σtr-state.σ)/norm(De*dgdσ)/10
    Δλ     = eps()

    σ  = σtr - Δλ*(De*dgdσ)

    εcp = state.εcp
    εtp = state.εtp
    εvp = state.εvp

    f   = yield_func(mat, state, state.σ, εtp, εcp, εvp)
    η   = 1.0 # initial damping

    # iterative process
    for i in 1:maxits
        dfdσ, _ = yield_derivs(mat, state, σ, εtp, εcp, εvp)
        dgdσ    = potential_derivs(mat, state, σ, εtp, εcp, εvp)
        dfdΔλ   = -dfdσ'*De*dgdσ

        Δλ = Δλ - η*f/dfdΔλ
        if Δλ<0
            # Δλ = abs(Δλ)
            # @show Δλ
        end

        if isnan(Δλ)
            return state.σ, 0.0, 0.0, 0.0, 0.0, failure("CSCP: Δλ is NaN")
        end

        σ  = σtr - Δλ*(De*dgdσ)

        Λ   = eigvals(dgdσ, sort=false)
        εtp = state.εtp + Δλ*maximum(max.(0.0, Λ))
        εcp = state.εcp + Δλ*norm(min.(0.0, Λ))
        εvp = state.εvp + Δλ*sum(abs, min.(0.0, Λ))
        # εtp = state.εtp + Δλ*norm(max.(0.0, Λ))
        f   = yield_func(mat, state, σ, εtp, εcp, εvp)

        if abs(f) < tol
            Δλ < 0.0 && return σ, 0.0, 0.0, 0.0, 0.0, failure("CSCP: negative Δλ")

            # @show Δλ
            # @show f
            return σ, εtp, εcp, εvp, Δλ, success()
        end

        # dumping
        i>10 && (η = 0.6)
        i>15 && (η = 0.3)
        # i>20 && (η = 0.15)
    end

    return state.σ, 0.0, 0.0, 0.0, 0.0, failure("CSCP: maximum iterations reached")
end


function update_state(mat::CSCP, state::CSCPState, Δε::AbstractArray)
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    Δσtr = De*Δε
    σtr  = state.σ + Δσtr
    ftr  = yield_func(mat, state, σtr, state.εtp, state.εcp, state.εvp)

    Δλ  = 0.0
    tol = 1.0
    tol = 0.1


    if ftr < tol
        # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        state.σ, state.εtp, state.εcp, state.εvp, state.Δλ, status = calc_σ_εp_Δλ(mat, state, σtr)
        @assert state.εcp >= 0.0
        @assert state.εtp >= 0.0
        @assert state.εvp >= 0.0

        # @show objectid(state)
        Δσ = state.σ - σini
        # @show Δε
        # @show state.σ
        # @show Δσ
        # @show state.εtp, state.εcp, state.εvp
        # @show failed(status)

        failed(status) && return state.σ, status
    end

    state.ε += Δε
    Δσ       = state.σ - σini
    return Δσ, success()
end


function state_values(mat::CSCP, state::CSCPState)
    σ, ε  = state.σ, state.ε
    ρ  = √(2*J2(σ))
    ξ  = tr(σ)/√3
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
