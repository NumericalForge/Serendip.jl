# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export WillamWarnke


"""
    WillamWarnke(; E, nu, fc, epsc, η=2.2, ft, GF=NaN, wc=NaN,
                  ft_law=:hordijk, fc_law=:popovics, beta=1.15)

Linear‐elastic concrete with a Willam–Warnke yield surface and
nonlinear hardening/softening in compression and tension softening
regularized by fracture energy.

# Keyword arguments
- `E::Real`: Young’s modulus (> 0).
- `nu::Real`: Poisson’s ratio (0 ≤ ν < 0.5).
- `fc::Real`: Uniaxial compressive strength (< 0).
- `epsc::Real`: Strain at the compressive peak (< 0).
- `η::Real=2.2`: Shape parameter for the compressive curve (> 1).
- `ft::Real`: Uniaxial tensile strength (> 0).
- `GF::Real=NaN`: Tensile fracture energy (> 0). Use `GF` or `wc`.
- `wc::Real=NaN`: Critical crack opening (≥ 0). Use `wc` or `GF`.
- `ft_law::Symbol=:hordijk`: Tension softening law. `:hordijk` uses the Hordijk curve.
- `fc_law::Symbol=:popovics`: Compression law. `:popovics` uses a Popovics-type curve.
- `beta::Real=1.15`: Biaxial/uniaxial compressive strength factor (1 ≤ β ≤ 1.5).

Exactly one of `GF` or `wc` must be provided.

# Returns
A `WillamWarnke` material usable with 3D or plane-strain bulk elements.
Not compatible with plane stress.

# Notes
- Tension softening is regularized by `GF` or `wc` to keep energy dissipation mesh-objective.
- Compression response follows the selected `fc_law` shaped by `(fc, epsc, η)`.
- `beta` sets `fb = β·fc` internally.

# Example
```julia
mat = WillamWarnke(E=30e9, nu=0.2, fc=-30e6, epsc=-0.002, ft=3e6,
                   GF=120.0, ft_law=:hordijk, fc_law=:popovics, beta=1.15)
```
"""
mutable struct WillamWarnke<:Constitutive
    E::Float64
    ν::Float64
    fc::Float64
    εc::Float64
    η::Float64
    ft::Float64
    wc::Float64
    ft_law::Symbol
    ft_fun::Union{Nothing,AbstractSpline}
    fc_law::Symbol
    fc_fun::Union{Nothing,AbstractSpline}
    fb::Float64
    e::Float64

    function WillamWarnke(;
        E::Real    = NaN,
        nu::Real   = NaN,
        fc::Real   = NaN,
        epsc::Real = NaN,
        eta::Real  = 2.2,
        ft::Real   = NaN,
        GF::Real   = NaN,
        wc::Real   = NaN,
        ft_law     = :constant,
        fc_law     = :constant,
        beta::Real = 1.15,
    )
        @check E>0 "WillamWarnke: Young's modulus E must be > 0. Got $E."
        @check 0<=nu<0.5 "WillamWarnke: Poisson's ratio nu must be in the range [0, 0.5). Got $nu."
        @check 1<=beta<=1.5 "WillamWarnke: Factor beta must be in the range [1.0, 1.5]. Got $beta."
        @check eta>1 "WillamWarnke: Shape parameter eta must be > 1. Got $eta."

        wc, ft_law, ft_fun, status = setup_tensile_strength(ft, GF, wc, ft_law)
        failed(status) && throw(ArgumentError("WillamWarnke: " * status.message))

        fc_law, fc_fun, status = setup_compressive_strength(fc, epsc, fc_law)
        failed(status) && throw(ArgumentError("WillamWarnke: " * status.message))

        fc_fun = nothing
        if fc_law isa AbstractSpline
            fc_fun = fc_law
            fc_law = :custom
            fc     = fc_law(0.0)
        end

        # biaxial strength
        fb = beta*fc

        e  = (fc*fb - fc*ft + 3*fb*ft)/(2*fc*fb + ft*fc)
        @assert 0<e<=1


        this = new(E, nu, fc, epsc, eta, ft, wc, ft_law, ft_fun, fc_law, fc_fun, fb, e)
        return this
    end
end


mutable struct WillamWarnkeState<:ConstState
    ctx::Context
    σ  ::Vec6
    ε  ::Vec6
    εtp::Float64
    εcp::Float64
    εvp::Float64
    Δλ ::Float64
    h  ::Float64
    function WillamWarnkeState(ctx::Context)
        this     = new(ctx)
        this.σ   = zeros(Vec6)
        this.ε   = zeros(Vec6)
        this.εtp = 0.0
        this.εcp = 0.0
        this.Δλ  = 0.0
        this.h   = 0.0
        this
    end
end


# Type of corresponding state structure
compat_state_type(::Type{WillamWarnke}, ::Type{MechBulk}) = ctx.stress_state!=:plane_stress ? WillamWarnkeState : error("WillamWarnke: This model is not compatible with planestress")


function calc_θ(::WillamWarnke, σ::Vec6)
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


function calc_rθ(mat::WillamWarnke, σ::Vec6)
    e = mat.e
    θ = calc_θ(mat, σ)

    rnum   = 2*(1-e^2)*cos(θ) + (2*e-1)*√(4*(1-e^2)*cos(θ)^2 + 5*e^2 - 4*e)
    rden   = 4*(1-e^2)*cos(θ)^2 + (2*e-1)^2
    r      = rnum/rden

    return r
end


function calc_fc(mat::WillamWarnke, εcp::Float64)
    fc0 = 0.35*mat.fc
    fcr = 0.1*mat.fc
    return calc_compressive_strength(mat, fc0, fcr, εcp)
end


function calc_ft(mat::WillamWarnke, w::Float64)
    return calc_tensile_strength(mat, w)
end


function calc_ξ0_κ(mat::WillamWarnke, state::WillamWarnkeState, εtp::Float64, εcp::Float64)

    e  = mat.e
    w  = εtp*state.h

    # current values
    ft = calc_ft(mat, w)
    fc = calc_fc(mat, εcp)

    ξ0 = (e+1)*ft*ft/ ( (4*e+1)*fc - 3*e*ft )
    κ  = √2*fc/(fc + √3*ξ0)

    return ξ0, κ
end


function yield_func(mat::WillamWarnke, state::WillamWarnkeState, σ::AbstractArray, εtp::Float64, εcp::Float64)
    # f(σ) = ρ - rθ⋅(ξ0-ξ)⋅κ

    i1, j2 = tr(σ), J2(σ)

    ξ = i1/√3
    ρ = √(2*j2)
    rθ = calc_rθ(mat, σ)
    ξ0, κ = calc_ξ0_κ(mat, state, εtp, εcp)

    return ρ - rθ*(ξ0 - ξ)*κ
end


function yield_derivs(mat::WillamWarnke, state::WillamWarnkeState, σ::AbstractArray, εtp::Float64, εcp::Float64)

    e = mat.e
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

    ξ0, κ = calc_ξ0_κ(mat, state, εtp, εcp)


    # f derivatives
    dfdρ  = 1.0
    dfdξ  = κ*rθ
    dfdrθ = κ*(ξ0 - ξ)
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

    f_εcp  = εcp -> yield_func(mat, state, σ, εtp, εcp)
    dfdεcp = derive(f_εcp, εcp)

    f_εtp  = εtp -> yield_func(mat, state, σ, εtp, εcp)
    dfdεtp = derive(f_εtp, εtp)

    return dfdσ, dfdεtp, dfdεcp
end


function potential_derivs(mat::WillamWarnke, state::WillamWarnkeState, σ::AbstractArray, εtp::Float64, εcp::Float64)
    # g(σ) = ρ - (ξ0-ξ)⋅κ Drucker-Prager like surface

    ξ0, κ = calc_ξ0_κ(mat, state, εtp, εcp)

    # ᾱ = mat.ᾱ
    s = dev(σ)

    # g derivatives
    # dgdρ = 1.0
    # dgdξ = κ
    # dρdσ = s/norm(s)
    # dξdσ = √3/3*I2

    dgdσ = s/norm(s) + κ*√3/3*I2
    return dgdσ

end


function calcD(mat::WillamWarnke, state::WillamWarnkeState)
    De  = calcDe(mat.E, mat.ν, state.ctx.stress_state)

    state.Δλ==0.0 && return De

    dfdσ, dfdεtp, dfdεcp = yield_derivs(mat, state, state.σ, state.εtp, state.εcp)
    dgdσ = potential_derivs(mat, state, state.σ, state.εtp, state.εcp)

    Λ = eigvals(dgdσ, sort=false)
    # Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεcp*norm(min.(0.0, Λ)) - dfdεtp*norm(max.(0.0, Λ)))
    Dep = De - De*dgdσ*dfdσ'*De / (dfdσ'*De*dgdσ - dfdεcp*norm(min.(0.0, Λ)) - dfdεtp*maximum(max.(0.0, Λ)))
    
    return Dep
end


function calc_σ_εp_Δλ(mat::WillamWarnke, state::WillamWarnkeState, σtr::Vec6)
    maxits = 60
    tol    = 0.1
    tol    = 1.0
    dgdσ   = potential_derivs(mat, state, state.σ, state.εtp, state.εcp)
    De     = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    Δλ     = eps()

    σ  = σtr - Δλ*(De*dgdσ)

    εcp = state.εcp
    εtp = state.εtp

    f   = yield_func(mat, state, state.σ, εtp, εcp)
    η   = 1.0 # initial damping

    # iterative process
    for i in 1:maxits
        dfdσ, _ = yield_derivs(mat, state, σ, εtp, εcp)
        dgdσ    = potential_derivs(mat, state, σ, εtp, εcp)
        dfdΔλ   = -dfdσ'*De*dgdσ

        Δλ = Δλ - η*f/dfdΔλ
        if Δλ<0
            # Δλ = abs(Δλ)
            # @show Δλ
        end

        if isnan(Δλ)
            return state.σ, 0.0, 0.0, 0.0, failure("WillamWarnke: Δλ is NaN")
        end

        σ  = σtr - Δλ*(De*dgdσ)

        Λ   = eigvals(dgdσ, sort=false)
        εtp = state.εtp + Δλ*maximum(max.(0.0, Λ))
        εcp = state.εcp + Δλ*norm(min.(0.0, Λ))
        f   = yield_func(mat, state, σ, εtp, εcp)

        if abs(f) < tol
            Δλ < 0.0 && return σ, 0.0, 0.0, 0.0, 0.0, failure("WillamWarnke: negative Δλ")

            # @show Δλ
            # @show f
            return σ, εtp, εcp, Δλ, success()
        end

        # dumping
        i>10 && (η = 0.6)
        i>15 && (η = 0.3)
    end

    return state.σ, 0.0, 0.0, 0.0, failure("WillamWarnke: maximum iterations reached")
end


# function calc_σ_εpa_Δλ(mat::WillamWarnke, state::WillamWarnkeState, σtr::Vec6)

#     α  = mat.α
#     ᾱ  = mat.ᾱ

#     E, ν = mat.E, mat.ν
#     K, G  = E/(3.0*(1.0-2.0*ν)), E/(2.0*(1.0+ν))

#     maxits = 20
#     tol = 1e-3

#     # ftr = yield_func(mat, state, σtr, state.εpa)

#     # @show f
#     str   = dev(σtr)
#     n_str = norm(str)
#     ξtr   = tr(σtr)/√3
#     ρtr   = n_str
#     r     = calc_rθ(mat, σtr)
#     ξ0  = calc_ξmax(mat, state.εpa)

#     # iterative process since θ is unknown at step η+1
#     for i in 1:maxits

#         numΔλ = ρtr - α*r*(ξ0 - ξtr)
#         denΔλ = 2*G + 3*α*ᾱ*K*r

#         Δλ = numΔλ/denΔλ

#         # @assert Δλ >= 0.0
#         # todo: check if Δλ is negative and return to the apex

#         # estimative at η+1
#         s   = (1 - 2*G*Δλ/n_str)*str
#         ξ   = ξtr - 3*ᾱ*K*Δλ
#         σ   = ξ/√3*I2 + s
#         εpa = state.εpa + Δλ*√(1 + ᾱ^2)
#         ρ   = ρtr - 2*G*Δλ
#         r   = calc_rθ(mat, σ)
#         ξ0 = calc_ξmax(mat, εpa)

#         # yield function at η+1
#         f = ρ - r*α*(ξ0 - ξ)

#         if abs(f) < tol
#             @assert Δλ >= 0.0

#             return σ, εpa, Δλ, success()
#         end

#     end


#     return state.σ, 0.0, 0.0, failure("WillamWarnke: maximum iterations reached")

# end


# function calc_σ_εpa_Δλ_bis(mat::WillamWarnke, state::WillamWarnkeState, σtr::Vec6)
#     α  = mat.α
#     ᾱ  = mat.ᾱ

#     E, ν = mat.E, mat.ν
#     K, G  = E/(3.0*(1.0-2.0*ν)), E/(2.0*(1.0+ν))

#     # De = calcDe(mat.E, mat.ν, state.ctx.stress_state)

#     ξtr = tr(σtr)/√3
#     str = dev(σtr)
#     n_str = norm(str)
#     r     = calc_rθ(mat, σtr)
#     ξ0  = calc_ξmax(mat, state.εpa)

#     # estimative of Δλ
#     # dgdσ = potential_derivs(mat, state, σtr)
#     # Δλ0  = norm(σtr-state.σ)/norm(De*dgdσ)

#     # str   = dev(σtr)
#     ρtr   = n_str
#     numΔλ = ρtr - α*r*(ξ0 - ξtr)
#     denΔλ = 2*G + 3*α*ᾱ*K*r
#     Δλ0 = numΔλ/denΔλ

#     # function of Δλ
#     ff(Δλ)  = begin
#         # quantities at η+1
#         ρ = ρtr - 2*G*Δλ
#         if ρ>0
#             ξ = ξtr - 3*ᾱ*K*Δλ
#             s = (1 - 2*G*Δλ/n_str)*str  # todo: avoid to compute s
#             σ = ξ/√3*I2 + s
#         else
#             ξ = calc_ξmax(mat, state.εpa)
#             σ = ξ/√3*I2
#         end

#         εpa = state.εpa + Δλ*√(1 + ᾱ^2)
#         return yield_func(mat, state, σ, εpa)
#     end

#     a, b, status = findrootinterval(ff, 0.0, Δλ0)
#     failed(status) && return state.σ, 0.0, 0.0, status

#     Δλ, status = findroot(ff, a, b, ftol=1e-3, method=:bisection)
#     failed(status) && return state.σ, 0.0, 0.0, status
#     @assert Δλ >= 0.0

#     σ   = σtr - 2*G*Δλ*str/n_str - √3*K*Δλ*ᾱ*I2
#     εpa = state.εpa + Δλ*√(1 + ᾱ^2)

#     return σ, εpa, Δλ, success()

# end


function update_state(mat::WillamWarnke, state::WillamWarnkeState, Δε::AbstractArray)
    σini = state.σ
    De   = calcDe(mat.E, mat.ν, state.ctx.stress_state)
    σtr  = state.σ + De*Δε
    ftr  = yield_func(mat, state, σtr, state.εtp, state.εcp)

    Δλ  = 0.0
    tol = 1.0
    tol = 0.1

    if ftr < 1e-8
        # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        # plastic
        # σ, εpa, Δλ, status = calc_σ_εpa_Δλ(mat, state, σtr)
        # σ, εpa, Δλ, status = calc_σ_εpa_Δλ_bis(mat, state, σtr)
        state.σ, state.εtp, state.εcp, state.Δλ, status = calc_σ_εp_Δλ(mat, state, σtr)
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


function state_values(mat::WillamWarnke, state::WillamWarnkeState)
    σ, ε  = state.σ, state.ε
    ρ = √(2*J2(σ))
    ξ = tr(σ)/√3

    w  = state.εtp*state.h
    ft = calc_ft(mat, w)
    fc = calc_fc(mat, state.εcp)

    vals_d = stress_strain_dict(σ, ε, state.ctx.stress_state)
    vals_d[:εcp] = state.εcp
    vals_d[:εtp] = state.εtp
    vals_d[:ξ]   = ξ
    vals_d[:ρ]   = ρ
    vals_d[:fc]  = fc
    vals_d[:ft]  = ft

    return vals_d
end
