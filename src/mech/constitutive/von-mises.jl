export VonMises

"""
    VonMises(; E, nu=0.0, fy, H=0.0)

Linear-elastic constitutive model with Von Mises yield criterion and linear isotropic hardening.
Implements J2 (pressure-insensitive) plasticity with associated flow rule.

# Arguments
- `E::Float64`: Young’s modulus (> 0.0).
- `nu::Float64`: Poisson’s ratio (0.0 ≤ ν < 0.5).
- `fy::Float64`: Initial yield stress (> 0.0).
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
- `VonMisesBeamState`: Beam elements (axial and transverse-shear stress/strain).
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


mutable struct VonMisesState<:ConstState
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


mutable struct VonMisesBeamState<:ConstState
    ctx::Context
    σ::Vec3
    ε::Vec3
    n::Vec3
    εpa::Float64
    Δλ::Float64
    αs::Float64
    function VonMisesBeamState(ctx::Context)
        this = new(ctx)
        this.σ   = zeros(Vec3)
        this.ε   = zeros(Vec3)
        this.n   = zeros(Vec3)
        this.εpa = 0.0
        this.Δλ  = 0.0
        this.αs  = 1.0
        this
    end
end


mutable struct VonMisesBarState<:ConstState
    ctx::Context
    σ::Float64
    ε::Float64
    εpa::Float64
    Δλ::Float64
    function VonMisesBarState(ctx::Context; σ::Float64=0.0)
        this = new(ctx)
        this.σ   = σ
        this.ε   = 0.0
        this.εpa = 0.0
        this.Δλ  = 0.0
        this
    end
end


compat_state_type(::Type{VonMises}, ::Type{MechSolid}) = VonMisesState
compat_state_type(::Type{VonMises}, ::Type{MechShell}) = VonMisesState

compat_state_type(::Type{VonMises}, ::Type{MechBeam}) = VonMisesBeamState
compat_state_type(::Type{VonMises}, ::Type{MechBar}) = VonMisesBarState
compat_state_type(::Type{VonMises}, ::Type{MechEmbBar}) = VonMisesBarState


# ❱❱❱ VonMises model for 3D and 2D bulk elements and shell elements

function yield_func(mat::VonMises, state::VonMisesState, σ::Vec6, εpa::Float64)
    j2d = J2(σ)
    σy  = mat.σy
    H   = mat.H
    return √(3*j2d) - σy - H*εpa
end


function calcD(mat::VonMises, state::VonMisesState)
    
    if state.ctx.stress_state==:plane_stress || state.αs!=1.0
        αs = state.αs
        De = calcDe(mat.E, mat.ν, :plane_stress, αs)
    else
        De = calcDe(mat.E, mat.ν)
    end
    
    state.Δλ==0.0 && return De

    j2d = J2(state.σ)
    @assert j2d>0
    
    σ = state.σ
    p = 1/3*(σ[1] + σ[2] + σ[3])
    s = SVector( σ[1]-p, σ[2]-p, σ[3]-p, σ[4], σ[5], σ[6] )

    dfdσ  = s*(√1.5/norm(s))
    # For the normalized von Mises flow in Mandel notation,
    # √(2/3)‖Δεp‖ = Δλ; hence dεpa/dλ = 1 for every stress state.
    dεpdλ = 1.0
    dfdεp = -mat.H

    return De - De*dfdσ*dfdσ'*De / (dfdσ'*De*dfdσ - dfdεp*dεpdλ)
end


function update_state(mat::VonMises, state::VonMisesState, cstate::VonMisesState, Δε::Vector{Float64})
    if state.ctx.stress_state==:plane_stress || state.αs!=1.0
        De = calcDe(mat.E, mat.ν, :plane_stress, state.αs)
    else
        De = calcDe(mat.E, mat.ν)
    end

    σtr  = cstate.σ + De*Δε
    ftr  = yield_func(mat, state, σtr, cstate.εpa)
    tol  = mat.σy*1e-8

    if ftr < tol
        state.Δλ = 0.0
        state.σ  = σtr
        Δσ = state.σ - cstate.σ
    else
        if state.ctx.stress_state==:plane_stress || state.αs!=1.0
            status = plastic_update_plane_stress(mat, state, cstate, σtr)
            failed(status) && return state.σ, status
    
            Δσ = state.σ - cstate.σ
    
            # # out-of-plane stress and strain components for plane stress condition TODO: Mandel
            # s33  = -1/3 * (state.σ[1] + state.σ[2])
            # Δε33 = -(mat.ν / mat.E) * (Δσ[1] + Δσ[2]) + state.Δγ * s33
            
            # update Δε
            # Δε = Vec6(Δε[1], Δε[2], Δε33, Δε[4], Δε[5], Δε[6])
    
            # return Δσ, success()
        else
            E, ν  = mat.E, mat.ν
            G     = E/(2*(1+ν))
            j2tr  = J2(σtr)

            Δλ = ftr/(3*G + mat.H)
            √j2tr - Δλ*√3*G >= 0.0 || return state.σ, failure("VonMisses: Negative value for √J2")

            s         = (1 - √3*G*Δλ/√j2tr)*dev(σtr)
            state.σ   = σtr - √6*G*Δλ*s/norm(s)
            # Accumulated equivalent plastic strain:
            # √(2/3)‖Δεp‖ = Δλ for the normalized flow direction.
            state.εpa = cstate.εpa + Δλ
            state.Δλ  = Δλ
    
            Δσ = state.σ - cstate.σ
        end
    end

    state.ε = cstate.ε + Δε
    
    return Δσ, success()
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


function plastic_update_plane_stress(mat::VonMises, state::VonMisesState, cstate::VonMisesState, σtr::Vec6)

    maxits = 50
    tol    = 1e-6*mat.σy
    Δλ     = 0.0 # plastic multiplier increment
    Δγ     = 0.0 # auxiliary variable
    c      = 2/3
    Δγmax  = mat.H > 0.0 ? 3/(2*mat.H) : Inf
    
    αs   = state.αs
    De   = calcDe(mat.E, mat.ν, :plane_stress, αs)
    invA = I4
    σ    = σtr
    εpa  = cstate.εpa

    for i in 1:maxits
        den = 1.0 - c*mat.H*Δγ
        den > 0.0 || return failure("VonMises: Invalid plane-stress plastic multiplier")

        invA = inv(I4 + Δγ*De*Psd)
        σ    = invA*σtr
        σvm = √(3*J2(σ))
        εpa = (cstate.εpa + c*mat.σy*Δγ)/den
        Δλ  = εpa - cstate.εpa

        R = yield_func(mat, state, σ, εpa)
        if abs(R) <= tol
            state.σ   = σ
            state.εpa = εpa
            state.Δλ  = Δλ
            return success()
        end
        
        p      = 1/3*(σ[1] + σ[2] + σ[3])
        s      = SVector( σ[1]-p, σ[2]-p, σ[3]-p, σ[4], σ[5], σ[6] )
        ∂σvm∂σ = √1.5*s/norm(s)
        ∂R∂εp  = -mat.H
        ∂σ∂Δγ  = -invA*De*s
        ∂εp∂Δγ = c*(mat.σy + mat.H*cstate.εpa)/den^2

        ∂R∂Δγ  = dot(∂σvm∂σ, ∂σ∂Δγ) + ∂R∂εp*∂εp∂Δγ
        isfinite(∂R∂Δγ) && ∂R∂Δγ!=0.0 || return failure("VonMises: Invalid plane-stress residual derivative")

        Δγnew = max(Δγ - R/∂R∂Δγ, 0.0)
        if Δγnew >= Δγmax
            Δγnew = 0.5*(Δγ + Δγmax)
        end
        Δγ = Δγnew

    end

    return failure("VonMises: plastic update failed")
end


function calc_σ_εpa_plane_stress(mat::VonMises, state::VonMisesState, cstate::VonMisesState, σtr::Vec6, Δγ::Float64)
    E, ν = mat.E, mat.ν
    G    = state.αs*E/2/(1+ν)

    # σ at n+1
    den = E^2*Δγ^2 - 2*E*ν*Δγ + 4*E*Δγ - 3*ν^2 + 3
    m11 = (2*E*Δγ - E*ν*Δγ - 3*ν^2 + 3)/den
    m12 = (E*Δγ - 2*E*ν*Δγ)/den
    m66 = 1/(2*G*Δγ + 1)

    σ = SVector(
        m11*σtr[1] + m12*σtr[2],
        m12*σtr[1] + m11*σtr[2],
        0.0,
        m66*σtr[4],
        m66*σtr[5],
        m66*σtr[6]
    )

    denε = 1.0 - 2/3*mat.H*Δγ
    denε > 0.0 || return σ, Inf
    εpa = (cstate.εpa + 2/3*mat.σy*Δγ)/denε

    return σ, εpa
end


# ❱❱❱ VonMises model for beam elements ❱❱❱

function yield_func(mat::VonMises, state::VonMisesBeamState, σ::Vec3, εpa::Float64)
    # Using Mandel's notation
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

    state.Δλ == 0.0 && return De

    σ   = state.σ
    σvm = √(σ[1]^2 + 3/2*(σ[2]^2 + σ[3]^2) )
    Q    = Vec3(1.0, 1.5, 1.5)
    n    = (Q.*σ)/σvm
    De_n = De*n

    return De - (De_n*De_n')/(dot(n, De_n) + mat.H)
end


function update_state(mat::VonMises, state::VonMisesBeamState, cstate::VonMisesBeamState, Δε::Vector{Float64})
    E, ν = mat.E, mat.ν
    G    = state.αs*E/2/(1+ν)
    De   = Vec3(E, 2*G, 2*G)
    
    σtr = cstate.σ + De.*Δε
    ftr = yield_func(mat, state, σtr, cstate.εpa)
    tol = 1e-8*mat.σy
    
    if ftr<tol
        state.Δλ = 0.0
        state.σ  = σtr
        state.n  = cstate.n
    else
        status = plastic_update_beam(mat, state, cstate, σtr)
        failed(status) && return state.σ, status
    end

    state.ε = cstate.ε + Δε
    Δσ      = state.σ - cstate.σ

    return Δσ, success()
end


function plastic_update_beam(mat::VonMises, state::VonMisesBeamState, cstate::VonMisesBeamState, σtr::Vec3)
    E, ν = mat.E, mat.ν
    G     = state.αs*E/2/(1+ν)
    De    = Vec3(E, 2*G, 2*G)
    Pd    = Vec3(2/3, 1.0, 1.0)
    Q     = De.*Pd

    c      = 2/3
    Δγ     = 0.0
    Δγmax  = mat.H > 0.0 ? 3/(2*mat.H) : Inf
    maxits = 50
    tol    = 1e-6*mat.σy

    σ    = σtr
    εpa  = cstate.εpa
    Δλ   = 0.0

    for i in 1:maxits
        den = 1.0 - c*mat.H*Δγ
        den > 0.0 || return failure("VonMises: Invalid beam plastic multiplier")

        A    = 1.0 .+ Δγ.*Q
        σ    = σtr./A
        σvm = √(σ[1]^2 + 1.5*(σ[2]^2 + σ[3]^2) )
        εpa = (cstate.εpa + c*mat.σy*Δγ)/den
        Δλ  = εpa - cstate.εpa

        R = yield_func(mat, state, σ, εpa)
        if abs(R) <= tol
            state.σ   = σ
            state.n   = Vec3(σ[1], 1.5*σ[2], 1.5*σ[3])/σvm
            state.εpa = εpa
            state.Δλ  = Δλ
            return success()
        end

        g       = Vec3(σ[1], 1.5*σ[2], 1.5*σ[3])/σvm
        ∂σ∂Δγ   = -(Q.*σ)./A
        ∂εp∂Δγ  = c*(mat.σy + mat.H*cstate.εpa)/den^2
        ∂R∂Δγ   = dot(g, ∂σ∂Δγ) - mat.H*∂εp∂Δγ
        isfinite(∂R∂Δγ) && ∂R∂Δγ!=0.0 || return failure("VonMises: Invalid beam residual derivative")

        Δγnew = max(Δγ - R/∂R∂Δγ, 0.0)
        if Δγnew >= Δγmax
            Δγnew = 0.5*(Δγ + Δγmax)
        end
        Δγ = Δγnew
    end

    return failure("VonMises: plastic update failed")
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


# ❱❱❱ Von Mises for bar elements ❱❱❱

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


function update_state(mat::VonMises, state::VonMisesBarState, cstate::VonMisesBarState, Δε::Float64)
    E, H = mat.E, mat.H
    σtr  = cstate.σ + E*Δε
    ftr  = yield_func(mat, state, σtr, cstate.εpa)

    if ftr<0
        state.Δλ = 0.0
        state.σ   = σtr
    else
        state.Δλ  = ftr/(E+H)
        Δεp       = state.Δλ*sign(σtr)
        state.εpa = cstate.εpa + state.Δλ
        state.σ   = σtr - E*Δεp
    end

    Δσ       = state.σ - cstate.σ
    state.ε  = cstate.ε + Δε
    return Δσ, success()
end


function state_values(mat::VonMises, state::VonMisesBarState)
    return OrderedDict{Symbol,Float64}(
        :σx´  => state.σ,
        :εx´  => state.ε,
        :εp => state.εpa,
    )
end
