export VonMises

"""
    VonMises(; E, nu=0.0, fy, H=0.0)

Linear-elastic constitutive model with Von Mises yield criterion and linear isotropic hardening.
Implements J2 (pressure-insensitive) plasticity with associated flow rule.

# Arguments
- `E::Float64`: Young’s modulus (must be > 0.0).
- `nu::Float64`: Poisson’s ratio (0.0 ≤ ν < 0.5).
- `fy::Float64`: Initial yield stress (must be > 0.0).
- `H::Float64`: Hardening modulus (≥ 0.0). A value of 0.0 corresponds to perfect plasticity.

# State Variables
Stored in `VonMisesState` (and its variants for reduced kinematics):
- `σ`: Stress tensor (full `Vec6` for 3D, reduced forms for plane stress, beam, and bar).
- `ε`: Strain tensor (same format as `σ`).
- `εpa::Float64`: Accumulated plastic strain.
- `Δλ::Float64`: Plastic multiplier increment.

# Variants
- `VonMisesState`: 3D continuum elements (full stress/strain in Voigt notation).
- `VonMisesPlaneStressState`: Plane stress elements.
- `VonMisesBeamState`: Beam elements (axial + bending stress/strain).
- `VonMisesBarState`: Truss elements (uniaxial stress/strain).
"""
mutable struct VonMises<:Constitutive
    E ::Float64
    ν ::Float64
    σy::Float64
    H ::Float64

    function VonMises(;
        E::Real=NaN,
        nu::Real=0.0,
        fy::Real=0.0,
        H::Real=0.0,
    )
        @check E > 0.0 "VonMises: Young's modulus E must be > 0.0. Got $E."
        @check nu >= 0.0 && nu < 0.5 "VonMises: Poisson's ratio nu must be in the range [0.0, 0.5). Got $nu."
        @check fy > 0.0 "VonMises: Initial yield stress fy must be > 0.0. Got $fy."
        @check H >= 0.0 "VonMises: Hardening modulus H must be >= 0.0. Got $H."
        return new(E, nu, fy, H)
    end

end


mutable struct VonMisesState<:IpState
    ctx::Context
    σ::Vec6
    ε::Vec6
    εpa::Float64
    Δλ::Float64
    αs::Float64
    function VonMisesState(ctx::Context)
        this = new(ctx)
        this.σ   = zeros(Vec6)
        this.ε   = zeros(Vec6)
        this.εpa = 0.0
        this.Δλ  = 0.0
        this.αs  = 1.0
        this
    end
end


mutable struct VonMisesBeamState<:IpState
    ctx::Context
    σ::Vec3
    ε::Vec3
    εpa::Float64
    Δλ::Float64
    αs::Float64
    function VonMisesBeamState(ctx::Context)
        this = new(ctx)
        this.σ   = zeros(Vec3)
        this.ε   = zeros(Vec3)
        this.εpa = 0.0
        this.Δλ  = 0.0
        this.αs  = 1.0
        this
    end
end

mutable struct VonMisesBarState<:IpState
    ctx::Context
    σ::Float64
    ε::Float64
    εpa::Float64
    Δλ::Float64
    function VonMisesBarState(ctx::Context; σ=0.0, ε=0.0)
        this = new(ctx)
        this.σ   = σ
        this.ε   = ε
        this.εpa = 0.0
        this.Δλ  = 0.0
        this
    end
end


compat_state_type(::Type{VonMises}, ::Type{MechBulk}) = VonMisesState
compat_state_type(::Type{VonMises}, ::Type{MechShell}) = VonMisesState

compat_state_type(::Type{VonMises}, ::Type{MechBeam}) = VonMisesBeamState
compat_state_type(::Type{VonMises}, ::Type{MechBar}) = VonMisesBarState
compat_state_type(::Type{VonMises}, ::Type{MechEmbBar}) = VonMisesBarState


# ❱❱❱ VonMises model for 3D and 2D bulk elements and shell elements

function yield_func(mat::VonMises, state::VonMisesState, σ::Vec6, εpa::Float64)
    if state.ctx.stress_state==:plane_stress || state.αs!=1.0
        # f = 1/2 σ*Psd*σ - 1/3 (fy + H εp)^2
        # f = J2D - 1/3 (fy + H εp)^2

        j2d = J2(σ)

        σy  = mat.σy
        H   = mat.H
        return j2d - 1/3*(σy + H*εpa)^2
    else
        j2d = J2(σ)
        σy  = mat.σy
        H   = mat.H
        return √(3*j2d) - σy - H*εpa
    end
end


function calcD(mat::VonMises, state::VonMisesState)
    if state.ctx.stress_state==:plane_stress || state.αs!=1.0
        αs  = state.αs
        De = calcDe(mat.E, mat.ν, :plane_stress, αs)
        state.Δλ==0.0 && return De
        σ = state.σ

        s     = SVector( 2/3*σ[1] - 1/3*σ[2], 2/3*σ[2] - 1/3*σ[1], -1/3*σ[1]-1/3*σ[2], σ[4], σ[5], σ[6] )
        dfdσ  = s
        dfdεp = -2/3*mat.H*(mat.σy + mat.H*state.εpa)

        return De - De*dfdσ*dfdσ'*De / (dfdσ'*De*dfdσ - norm(s)*dfdεp)

    else
        De  = calcDe(mat.E, mat.ν)
        state.Δλ==0.0 && return De

        j2d = J2(state.σ)
        @assert j2d>0

        s     = dev(state.σ)
        dfdσ  = √1.5*s/norm(s)
        dfdεp = -mat.H

        return De - De*dfdσ*dfdσ'*De / (dfdσ'*De*dfdσ - √1.5*dfdεp)
    end

end


function update_state(mat::VonMises, state::VonMisesState, Δε::Vector{Float64})
    if state.ctx.stress_state==:plane_stress || state.αs!=1.0
        σini = state.σ

        αs   = state.αs
        De   = calcDe(mat.E, mat.ν, :plane_stress, αs)
        σtr  = state.σ + De*Δε
        ftr  = yield_func(mat, state, σtr, state.εpa)
        tol  = 1e-8

        if ftr<tol
            # elastic
            state.Δλ = 0.0
            state.σ  = σtr
        else
            # plastic
            σ, εpa, Δλ, status = calc_σ_εpa_Δλ_plane_stress(mat, state, σtr)
            failed(status) && return state.σ, status

            state.σ, state.εpa, state.Δλ = σ, εpa, Δλ
        end

        state.ε += Δε
        Δσ     = state.σ - σini

        return Δσ, success()
    else
        σini = state.σ
        De   = calcDe(mat.E, mat.ν)
        σtr  = state.σ + De*Δε
        ftr  = yield_func(mat, state, σtr, state.εpa)
        tol  = 1e-8

        if ftr<tol
            # elastic
            state.Δλ = 0.0
            state.σ  = σtr
        else
            # plastic
            E, ν  = mat.E, mat.ν
            G     = E/(2*(1+ν))
            j2dtr = J2(σtr)

            state.Δλ = ftr/(3*G + √1.5*mat.H)
            √j2dtr - state.Δλ*√3*G >= 0.0 || return state.σ, failure("VonMisses: Negative value for √J2D")

            s = (1 - √3*G*state.Δλ/√j2dtr)*dev(σtr)
            state.σ = σtr - √6*G*state.Δλ*s/norm(s)
            state.εpa += state.Δλ
        end

        state.ε += Δε
        Δσ     = state.σ - σini
        return Δσ, success()
    end
end


function state_values(mat::VonMises, state::VonMisesState)
    σ, ε  = state.σ, state.ε
    # j1    = tr(σ)
    # srj2d = √J2(σ)

    stress_state = state.αs==1.0 ? state.ctx.stress_state : :plane_stress
    D = stress_strain_dict(σ, ε, stress_state)
    D[:εp]   = state.εpa

    return D
end


function calc_σ_εpa_Δλ_plane_stress(mat::VonMises, state::VonMisesState, σtr::Vec6)
    # Δλ estimative
    # De   = calcDe(mat.E, mat.ν, :plane_stress)
    # dfdσ = SVector( 2/3*σtr[1] - 1/3*σtr[2], 2/3*σtr[2] - 1/3*σtr[1], 0.0, σtr[4], σtr[5], σtr[6] )
    dfdσ = SVector( 2/3*σtr[1] - 1/3*σtr[2], 2/3*σtr[2] - 1/3*σtr[1], -1/3*σtr[1]-1/3*σtr[2], σtr[4], σtr[5], σtr[6] )

    # Δλ0  = norm(σtr-state.σ)/norm(De*dfdσ)
    Δλ0  = norm(σtr-state.σ)/(mat.E*norm(dfdσ))

    # find initial interval
    a = 0.0
    b = Δλ0

    σ, εpa = calc_σ_εpa_plane_stress(mat, state, σtr, a)
    fa     = yield_func(mat, state, σ, εpa)
    σ, εpa = calc_σ_εpa_plane_stress(mat, state, σtr, b)
    fb     = yield_func(mat, state, σ, εpa)

    # search for a valid interval
    if fa*fb>0
        maxits = 50
        for i in 1:maxits
            b  += Δλ0*(1.6)^i
            σ, εpa = calc_σ_εpa_plane_stress(mat, state, σtr, b)
            fb     = yield_func(mat, state, σ, εpa)
            fa*fb<0.0 && break

            i==maxits && return state.σ, 0.0, 0.0, failure("VonMises: Could not find interval for Δλ")
        end
    end

    ff(Δλ) = begin
        σ, εpa = calc_σ_εpa_plane_stress(mat, state, σtr, Δλ)
        yield_func(mat, state, σ, εpa)
    end

    tol = 10^-(8-log10(mat.σy))


    # findroot
    # Δλ, status = findroot(ff, a, b, tol)
    # failed(status) && return state.σ, 0.0, 0.0, status

    # σ, εpa = calc_σ_εpa_plane_stress(mat, state, σtr, Δλ)

    # bissection method
    local f, Δλ, σ, εpa
    σ0  = zeros(SVector{6}) # initial value

    tol    = 10^-(10-log10(mat.σy))
    maxits = 50

    for i in 1:maxits
        Δλ = (a+b)/2
        σ, εpa = calc_σ_εpa_plane_stress(mat, state, σtr, Δλ)
        f = yield_func(mat, state, σ, εpa)

        if fa*f<0
            b = Δλ
        else
            a  = Δλ
            fa = f
        end

        maximum(abs, σ-σ0) <= tol && break
        σ0 = σ

        i==maxits && return state.σ, 0.0, 0.0, failure("VonMises: could not find Δλ with NR/bissection (maxits reached, f=$f)")
    end

    return σ, εpa, Δλ, success()
end


function calc_σ_εpa_plane_stress(mat::VonMises, state::VonMisesState, σtr::Vec6, Δλ::Float64)
    E, ν = mat.E, mat.ν
    G    = state.αs*E/2/(1+ν)

    # σ at n+1
    den = E^2*Δλ^2 - 2*E*ν*Δλ + 4*E*Δλ - 3*ν^2 + 3
    m11 = (2*E*Δλ - E*ν*Δλ - 3*ν^2 + 3)/den
    m12 = (E*Δλ - 2*E*ν*Δλ)/den
    m66 = 1/(2*G*Δλ + 1)

    σ = SVector(
        m11*σtr[1] + m12*σtr[2],
        m12*σtr[1] + m11*σtr[2],
        0.0,
        m66*σtr[4],
        m66*σtr[5],
        m66*σtr[6]
    )

    dfdσ = SVector( 2/3*σ[1] - 1/3*σ[2], 2/3*σ[2] - 1/3*σ[1], -1/3*σ[1]-1/3*σ[2], σ[4], σ[5], σ[6] )

    εpa  = state.εpa + Δλ*norm(dfdσ)

    return σ, εpa
end


# ❱❱❱ VonMises model for beam elements

function yield_func(mat::VonMises, state::VonMisesBeamState, σ::Vec3, εpa::Float64)
    # Using Mendel's notation
    # f = √(3 J2) - fy - H εp
    # σ = [ σ1, √2*σ2, √2*σ3 ]
    # s = [ 2/3*σ1, -1/3*σ1, -1/3*σ1, 0.0, √2*σ2, √2*σ3 ]
    # σvm = √(σ1^2 + 3/2 (σ2^2 + σ3^2))

    σvm = √(σ[1]^2 + 3/2*(σ[2]^2 + σ[3]^2) )

    return σvm - mat.σy - mat.H*εpa
end


function calcD(mat::VonMises, state::VonMisesBeamState)
    E, ν = mat.E, mat.ν
    G    = state.αs*E/2/(1+ν)
    De = @SMatrix [ E    0.0  0.0
                    0.0  2*G  0.0
                    0.0  0.0  2*G ]

    state.Δλ==0.0 && return De

    σ   = state.σ
    σvm = √(σ[1]^2 + 3/2*(σ[2]^2 + σ[3]^2) )
    σvm < 1e-14 && return De

    Q       = Vec3( 1, 3/2, 3/2 )
    dfdσ    = 1/σvm * (Q .* σ)
    De_vec  = Vec3(E, 2G, 2G)
    De_dfdσ = De_vec .* dfdσ


    return De - (De_dfdσ*De_dfdσ') / (dot(dfdσ, De_dfdσ) + mat.H)

end


function update_state(mat::VonMises, state::VonMisesBeamState, Δε::Vector{Float64})
    σini = state.σ

    E, ν = mat.E, mat.ν
    G    = state.αs*E/2/(1+ν)
    De   = Vec3(E, 2*G, 2*G)

    σtr = state.σ + De.*Δε
    ftr = yield_func(mat, state, σtr, state.εpa)
    tol = 1e-8

    if ftr<tol
        # elastic
        state.Δλ = 0.0
        state.σ  = σtr
    else
        σ, εpa, Δλ, status = calc_σ_εpa_Δλ(mat, state, σtr)
        failed(status) && return state.σ, status

        state.σ, state.εpa, state.Δλ = σ, εpa, Δλ
    end

    state.ε += Δε
    Δσ     = state.σ - σini

    return Δσ, success()
end


function calc_σ_εpa_Δλ(mat::VonMises, state::VonMisesBeamState, σtr::Vec3)
    # Δλ estimative
    E, ν = mat.E, mat.ν
    G  = state.αs*E/2/(1+ν)
    De = Vec3(E, 2*G, 2*G)
    Q  = Vec3( 1, 1.5, 1.5 )

    σ      = σtr
    σvm_tr = √(σ[1]^2 + 1.5*(σ[2]^2 + σ[3]^2) )
    n_tr   = (Q.*σtr)/σvm_tr # freezing the plastic direction
    
    Δλ     = 0.0
    maxits = 20
    tol    = 1e-6*mat.σy
    for _ in 1:maxits
        σvm = √(σ[1]^2 + 1.5*(σ[2]^2 + σ[3]^2) )
        g = (Q.*σ)/σvm
        εpa = state.εpa + Δλ
        
        R = yield_func(mat, state, σ, εpa)
        abs(R) <= tol && return σ, εpa, Δλ, success()
        
        Rp = -dot(g, De.*n_tr) - mat.H
        Δλ = max(Δλ - R/Rp, 0.0)
        σ  = σtr - Δλ*(De.*n_tr)
    end

    return state.σ, 0.0, 0.0, failure("VonMises: Could not find Δλ")
end


function state_values(mat::VonMises, state::VonMisesBeamState)
    σ = state.σ
    σvm = √(σ[1]^2 + 1.5*(σ[2]^2 + σ[3]^2) )

    vals = OrderedDict{Symbol,Float64}(
        :σx´   => state.σ[1],
        :εx´   => state.ε[1],
        :εp    => state.εpa,
        :σx´y´ => state.σ[3]/SR2, # x´y´ component is the third one
        :σx´z´ => state.σ[2]/SR2, # x´z´ component (3d)
        :σvm   => σvm
    )

    return vals
end


# ❱❱❱ Von Mises for bar elements

function yield_func(mat::VonMises, state::VonMisesBarState, σ::Float64, εpa::Float64)
    return abs(σ) - (mat.σy + mat.H*εpa)
end


function calcD(mat::VonMises, state::VonMisesBarState)
    if state.Δλ == 0.0
        return mat.E
    else
        E, H = mat.E, mat.H
        return E*H/(E+H)
    end
end


function update_state(mat::VonMises, state::VonMisesBarState, Δε::Float64)
    E, H    = mat.E, mat.H
    σini    = state.σ
    σtr     = σini + E*Δε
    ftr     = yield_func(mat, state, σtr, state.εpa)

    if ftr<0
        state.Δλ = 0.0
        state.σ   = σtr
    else
        state.Δλ  = ftr/(E+H)
        Δεp       = state.Δλ*sign(σtr)
        state.εpa += state.Δλ
        state.σ   = σtr - E*Δεp
    end

    Δσ        = state.σ - σini
    state.ε  += Δε
    return Δσ, success()
end


function state_values(mat::VonMises, state::VonMisesBarState)
    return OrderedDict{Symbol,Float64}(
        :σx´  => state.σ,
        :εx´  => state.ε,
        :εp => state.εpa,
    )
end