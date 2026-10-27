# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export TensionFieldSoftening

"""
    TensionFieldSoftening(; E, nu=0.0, fy, GF=NaN, beta=1.0)

Isotropic small-strain membrane plasticity with associated von Mises flow and
reversible tension-field relaxation (zero compressive principal stresses).
Without `GF`, the tensile response is perfectly plastic. Positive `GF` selects
`fy*max(0, 1-εpa/εcrit)^beta`, with `beta ≥ 1` and
`εcrit=(beta+1)*GF/(h*fy)`. MechMembrane sets `h=√(reference element area)`.
`GF` calibrates plastic dissipation per rupture area, not released elastic work.

Complete rupture is irreversible. Local failure requests global cutback; it
never marks rupture. Very steep softening may have no continuously accessible
return. The model uses neither viscosity nor local constitutive substepping.
"""
struct TensionFieldSoftening <: Constitutive
    E::Float64
    ν::Float64
    fy::Float64
    GF::Float64
    β::Float64

    function TensionFieldSoftening(; E::Real=NaN, nu::Real=0.0, fy::Real=0.0, GF::Real=NaN, beta::Real=1.0)
        @check isfinite(E) && E > 0.0 "TensionFieldSoftening: E must be finite and positive."
        @check isfinite(nu) && 0.0 <= nu < 0.5 "TensionFieldSoftening: nu must lie in [0, 0.5)."
        @check isfinite(fy) && fy > 0.0 "TensionFieldSoftening: fy must be finite and positive."
        @check isnan(GF) || (isfinite(GF) && GF > 0.0) "TensionFieldSoftening: GF must be positive when supplied."
        @check isfinite(beta) && beta >= 1.0 "TensionFieldSoftening: beta must be finite and at least one."
        return new(E, nu, fy, GF, beta)
    end
end

mutable struct TensionFieldSofteningState <: ConstState
    ctx::Context
    σ::Vec3
    ε::Vec3
    εp::Vec3
    ω::Vec3
    εpa::Float64
    Δλ::Float64
    Δγ::Float64
    h::Float64
    ruptured::Bool
    branch::Symbol

    function TensionFieldSofteningState(ctx::Context)
        return new(ctx, zeros(Vec3), zeros(Vec3), zeros(Vec3), zeros(Vec3),
            0.0, 0.0, 0.0, NaN, false, :initial)
    end
end

compat_state_type(::Type{TensionFieldSoftening}, ::Type{MechMembrane}) = TensionFieldSofteningState

function calcD(mat::TensionFieldSoftening, state::TensionFieldSofteningState)
    # The tensile-side elastic derivative is an admissible predictor at the
    # nonsmooth origin, without adding stiffness in accepted slack states.
    state.branch in (:initial, :taut_elastic) && return membrane_elastic_matrix(mat.E, mat.ν)
    state.branch in (:slack, :ruptured) && return zeros(3,3)
    E, ν = mat.E, mat.ν
    σ = state.σ
    q = sqrt(σ[1]^2-σ[1]*σ[2]+σ[2]^2+1.5σ[3]^2)
    gradient = Vec3(σ[1]-σ[2]/2, σ[2]-σ[1]/2, 1.5σ[3])/q
    # Recover the fixed-history trial strain from the accepted increment.
    e = state.ε-state.εp+state.Δλ*gradient
    F = eigen(Symmetric([e[1] e[3]/√2; e[3]/√2 e[2]]))
    V = F.vectors[:,end:-1:1]
    a, b = reverse(F.values)
    st = V'*[σ[1] σ[3]/√2; σ[3]/√2 σ[2]]*V
    s = [st[1,1],st[2,2]]
    εcrit = isnan(mat.GF) ? Inf : (mat.β+1)*mat.GF/(state.h*mat.fy)
    _, kp = tension_field_strength(mat,state.εpa,εcrit)
    if state.branch in (:wrinkled_elastic, :wrinkled_plastic)
        Et = state.Δλ == 0.0 ? E : E*kp/(E+kp)
        A = [Et 0.0; 0.0 0.0]
    else
        M = [1.0 -0.5; -0.5 1.0]
        P = inv([1.0 -ν; -ν 1.0]/E+state.Δγ*M)
        g = M*s/q
        ds = -P*M*s
        Rγ = (1-kp*state.Δγ)*dot(g,ds)-kp*q
        A = P-(1-kp*state.Δγ)/Rγ*(ds*(g'*P))
    end
    c, d = V[:,1]
    T = [c*c d*d √2*c*d; d*d c*c -√2*c*d; -√2*c*d √2*c*d c*c-d*d]
    scale = max(mat.fy/E, abs(a), abs(b))
    shear = abs(a-b) > 1e-10*scale ? (s[1]-s[2])/(a-b) : A[1,1]-A[1,2]
    principal = zeros(3,3)
    principal[1:2,1:2] = A
    principal[3,3] = shear
    return T'*principal*T
end

function tension_field_strength(mat::TensionFieldSoftening, εpa::Float64, εcrit::Float64)
    isinf(εcrit) && return mat.fy, 0.0
    εpa >= εcrit && return 0.0, 0.0
    z = 1.0-εpa/εcrit
    return mat.fy*z^mat.β, -mat.β*mat.fy/εcrit*z^(mat.β-1)
end

# Limiting plastic multiplier as stress tends to zero. The second case has
# positive limiting stresses in both principal directions; the first wrinkles.
function tension_field_terminal_strain(e::AbstractVector)
    a, b = reverse(LinearAlgebra.eigvals(Symmetric([e[1] e[3]/√2; e[3]/√2 e[2]])))
    a <= 0.0 && return 0.0
    b <= -a/2 && return a
    return sqrt(4/3*(a*a+a*b+b*b))
end

function update_state(mat::TensionFieldSoftening, state::TensionFieldSofteningState, cstate::TensionFieldSofteningState, Δε::AbstractVector)
    length(Δε) == 3 || throw(DimensionMismatch("TensionFieldSoftening: expected three Mandel strain components."))
    all(isfinite, Δε) || return zeros(Vec3), failure("TensionFieldSoftening: nonfinite strain increment.")
    ε = cstate.ε + Δε
    all(isfinite, ε) || return zeros(Vec3), failure("TensionFieldSoftening: nonfinite total strain.")
    εcrit = Inf
    if !isnan(mat.GF)
        isfinite(state.h) && state.h > 0.0 || return zeros(Vec3), failure("TensionFieldSoftening: a positive characteristic length h is required.")
        εcrit = (mat.β+1)*mat.GF/(state.h*mat.fy)
    end
    E, ν = mat.E, mat.ν
    εpa = cstate.εpa
    e = ε-cstate.εp
    F = eigen(Symmetric([e[1] e[3]/√2; e[3]/√2 e[2]]))
    V = F.vectors[:,end:-1:1]
    a, b = reverse(F.values)
    Δλ, Δγ = 0.0, 0.0
    stress = zeros(Vec3)
    plastic = cstate.εp
    wrinkle = zeros(Vec3)
    branch = :slack
    ruptured = cstate.ruptured
    tol = mat.fy*1e-10

    if ruptured
        branch = :ruptured
    elseif isfinite(εcrit) && tension_field_terminal_strain(e) >= εcrit-εpa
        # A zero-stress endpoint is not a root at finite Δγ. Certify that the
        # intact local map stays regular before accepting an endpoint crossing.
        # E/(2(1-ν)) is the smallest generalized eigenvalue of (M, Ce).
        _, kp = tension_field_strength(mat, εpa, εcrit)
        -kp < E/(2*(1-ν)) || return zeros(Vec3), failure("TensionFieldSoftening: rupture continuation is unresolved for this softening slope; reduce the increment.")
        e0 = cstate.ε-cstate.εp
        remaining = εcrit-εpa
        remaining > 0.0 && tension_field_terminal_strain(e0) < remaining ||
            return zeros(Vec3), failure("TensionFieldSoftening: no intact starting state for rupture continuation.")
        lo, hi = 0.0, 1.0
        for _ in 1:60
            mid = (lo+hi)/2
            if tension_field_terminal_strain(e0+mid*Δε) < remaining
                lo = mid
            else
                hi = mid
            end
        end
        endpoint = e0+hi*Δε
        ef = eigen(Symmetric([endpoint[1] endpoint[3]/√2; endpoint[3]/√2 endpoint[2]]))
        va, vb = reverse(ef.values)
        W = ef.vectors[:,end:-1:1]
        increments = vb <= -va/2 ? [remaining, -remaining/2] : [va, vb]
        ep = W*Diagonal(increments)*W'
        plastic = cstate.εp + Vec3(ep[1,1], ep[2,2], √2*ep[1,2])
        Δλ = remaining
        εpa = εcrit
        ruptured = true
        branch = :ruptured
    elseif a > 0.0
        Y, _ = tension_field_strength(mat, εpa, εcrit)
        # Uniaxial candidate; total strain minus committed plastic strain also
        # accounts for closure of the reversible wrinkles on reloading.
        if E*a > Y+tol
            f = λ -> E*(a-λ)-tension_field_strength(mat, εpa+λ, εcrit)[1]
            df = λ -> begin
                slope = -E-tension_field_strength(mat, εpa+λ, εcrit)[2]
                slope < 0.0 ? slope : NaN
            end
            Δλ, status = findroot_bracket_newton(f, df, 0.0;
                tol=1e-12*max(mat.fy/E, a), ftol=tol, maxits=80)
            failed(status) && return zeros(Vec3), failure("TensionFieldSoftening: wrinkled consistency solve failed; reduce the increment.")
            Δλ >= 0.0 && εpa+Δλ < εcrit || return zeros(Vec3), failure("TensionFieldSoftening: inadmissible wrinkled plastic root.")
        end
        s1 = E*(a-Δλ)
        isfinite(s1) && s1 > 0.0 || return zeros(Vec3), failure("TensionFieldSoftening: unresolved zero-stress wrinkled endpoint.")
        w = -ν*s1/E-Δλ/2-b
        s = [s1, 0.0]
        gp = [1.0, -0.5]
        if w >= 0.0
            _, kp = tension_field_strength(mat, εpa+Δλ, εcrit)
            Δλ == 0.0 || E+kp > 0.0 || return zeros(Vec3), failure("TensionFieldSoftening: singular wrinkled tangent.")
            branch = Δλ == 0.0 ? :wrinkled_elastic : :wrinkled_plastic
        else
            Ce = [1.0 -ν; -ν 1.0]/E
            M = [1.0 -0.5; -0.5 1.0]
            en = [a,b]
            trial = Ce\en
            qtrial = sqrt(dot(trial, M*trial))
            Δλ = 0.0
            w = 0.0
            if minimum(trial) >= 0.0 && qtrial <= Y+tol
                s = trial
                branch = :taut_elastic
            else
                # All stresses are eliminated analytically; Newton solves only
                # R(Δγ)=σVM-σmax(εpa+Δγ σVM)=0.
                function residual(γ)
                    P = inv(Ce+γ*M)
                    s = P*en
                    q = sqrt(dot(s,M*s))
                    g = M*s/q
                    ds = -P*M*s
                    strength, kp = tension_field_strength(mat, εpa+γ*q, εcrit)
                    R = q-strength
                    Rγ = (1-kp*γ)*dot(g,ds)-kp*q
                    return R, Rγ
                end
                f = γ -> residual(γ)[1]
                df = γ -> begin
                    slope = residual(γ)[2]
                    slope < 0.0 ? slope : NaN
                end
                Δγ, status = findroot_bracket_newton(f, df, 0.0;
                    tol=1e-12/E, ftol=tol, maxits=100)
                failed(status) && return zeros(Vec3), failure("TensionFieldSoftening: taut scalar consistency solve failed; reduce the increment.")
                Δγ >= 0.0 || return zeros(Vec3), failure("TensionFieldSoftening: negative auxiliary increment.")
                P = inv(Ce+Δγ*M)
                s = P*en
                q = sqrt(dot(s,M*s))
                Δλ = Δγ*q
                isfinite(q) && q > 0.0 && Δλ >= 0.0 && εpa+Δλ < εcrit && minimum(s) >= 0.0 ||
                    return zeros(Vec3), failure("TensionFieldSoftening: inadmissible taut plastic root.")
                gp = M*s/q
                R, Rγ = residual(Δγ)
                abs(R) <= tol && Rγ < 0.0 || return zeros(Vec3), failure("TensionFieldSoftening: singular or unconverged taut return.")
                branch = :taut_plastic
            end
        end
        st = V*Diagonal(s)*V'
        ep = V*Diagonal(Δλ*gp)*V'
        wt = V*Diagonal([0.0, max(0.0,w)])*V'
        stress = Vec3(st[1,1], st[2,2], √2*st[1,2])
        plastic = cstate.εp + Vec3(ep[1,1], ep[2,2], √2*ep[1,2])
        wrinkle = Vec3(wt[1,1], wt[2,2], √2*wt[1,2])
        εpa += Δλ
    else
        wrinkle = -e
    end

    # Assign only an accepted trial. Committed history is never mutated here.
    state.σ = stress
    state.ε = ε
    state.εp = plastic
    state.ω = wrinkle
    state.εpa = εpa
    state.Δλ = Δλ
    state.Δγ = Δγ
    state.ruptured = ruptured
    state.branch = branch
    return state.σ-cstate.σ, success()
end

function state_values(mat::TensionFieldSoftening, state::TensionFieldSofteningState)
    vals = membrane_state_values(state.σ, state.ε)
    vals[:εpxx] = state.εp[1]
    vals[:εpyy] = state.εp[2]
    vals[:εpxy] = state.εp[3]/√2
    vals[:εpzz] = -state.εp[1]-state.εp[2]
    vals[:εpa] = state.εpa
    vals[:ωxx] = state.ω[1]
    vals[:ωyy] = state.ω[2]
    vals[:ωxy] = state.ω[3]/√2
    vals[:ruptured] = Float64(state.ruptured)
    return vals
end
