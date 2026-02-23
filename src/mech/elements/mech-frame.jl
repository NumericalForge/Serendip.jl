# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechFrame

# MechFrame_params = [
#     FunInfo( :MechFrame, "A straight 2D frame (truss + simple beam) element"),
#     KwArgInfo( :E, "Young modulus",  cond=:(E>0.0)  ),
#     KwArgInfo( :A, "Section area",  cond=:(A>0.0)  ),
#     KwArgInfo( :I, "Moment of inertia",  cond=:(I>0.0)  ),
#     KwArgInfo( :gamma, "Specific weight", 0, cond=:(gamma>=0.0) ),
#     KwArgInfo( :rho, "Density", 0, cond=:(rho>=0.0)  ),
# ]
# @doc docstring(MechFrame_params) MechFrame


struct MechFrame<:MechFormulation
    E::Float64
    A::Float64
    I::Float64
    ρ::Float64
    γ::Float64

    function MechFrame(;
        E::Real=NaN,
        A::Real=NaN,
        I::Real=NaN,
        rho::Real=0.0,
        gamma::Real=0.0,
    )
        @check E > 0.0 "Young modulus must be positive"
        @check A > 0.0 "Section area must be positive"
        @check I > 0.0 "Moment of inertia must be positive"
        @check rho >= 0.0 "Density must be non-negative"
        @check gamma >= 0.0 "Specific weight must be non-negative"

        return new(E, A, I, rho, gamma)
    end
end


compat_role(::Type{MechFrame}) = :line
embedded_formulation(::Type{MechFrame}) = error("MechFrame: this element cannot be embedded")


function distributed_bc(elem::Element{MechFrame}, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return mech_line_distributed_forces(elem, t, key, val)
end


# function body_c(elem::Element{MechFrame}, key::Symbol, val::Union{Real,Symbol,Expr})
#     suitable_keys = (:qy,:qn, :wy)
#     key in suitable_keys || error("MechFrame: boundary condition $key is not applicable as distributed bc at element of type $(typeof(elem)). Suitable keys are $(string.(suitable_keys))")

#     return mech_line_distributed_forces(elem, 0.0, key, val)
# end

function elem_config_dofs(elem::Element{MechFrame})
    ndim = elem.ctx.ndim
    ndim==2 || error("MechFrame: Frame elements require ndim=2. Current ndim=$(ndim)")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :rz, :mz)
    end
end


function elem_map(elem::Element{MechFrame})
    keys =(:ux, :uy, :rz)
    return [ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]
end


dof_map(elem::Element{MechFrame}) = elem_map(elem)


function beam_shape_func(𝑥::Float64, ℓ::Float64)
    N = Array{Float64}(undef,6)
    N[1] = 1 - 𝑥/ℓ
    N[2] = 1 - 3*𝑥^2/ℓ^2 + 2*𝑥^3/ℓ^3
    N[3] = 𝑥 - 2*𝑥^2/ℓ + 𝑥^3/ℓ^2
    N[4] = 𝑥/ℓ
    N[5] = 3*𝑥^2/ℓ^2 - 2*𝑥^3/ℓ^3
    N[6] = 𝑥^3/ℓ^2 - 𝑥^2/ℓ

    return N
end


function body_c(elem::Element{MechFrame}, key::Symbol, val::Union{Real,Symbol,Expr})
    ndim  = elem.ctx.ndim

    # Check bcs
    !(key in (:qy, :qn)) && error("distributed_bc: boundary condition $key is not applicable as distributed bc at element with type $(typeof(elem))")

    nodes  = elem.nodes
    nnodes = length(nodes)
    t = 0.0

    # Force boundary condition
    nnodes = length(nodes)

    # Calculate the target coordinates matrix
    C = get_coords(nodes, ndim)
    L = norm(C[2,:]-C[1,:])

    # Calculate the nodal values
    F     = zeros(6)
    shape = LIN2
    ips   = get_ip_coords(shape)

    for i in 1:size(ips,1)
        R = ips[i].coord
        w = ips[i].w
        X = C'*LIN2.func(R)

        l = (C[2,:]-C[1,:])./L
        n = [-l[2], l[1]]

        if ndim==2
            x, y = X
            vip = evaluate(val, t=t, x=x, y=y)

            if key == :qy
                tl = vip*l[2]
                qn = vip*n[2]
            elseif key == :qn
                tl = 0.0
                qn = vip
            end
        else
            error("This beam element is for 1D only")
        end

        N = beam_shape_func(R[1]*L/2+L/2, L)
        Nl = [ N[1], 0, 0, N[4], 0, 0 ]
        Nn = [ 0, N[2], N[3], 0, N[5], N[6] ]

        F += (Nl*tl + Nn*qn)*L/2*w # F is a vector
    end

    # Rotation matrix
    c = (C[2,1] - C[1,1])/L
    s = (C[2,2] - C[1,2])/L
    T = [ c s 0  0 0 0
          -s c 0  0 0 0
          0 0 1  0 0 0
          0 0 0  c s 0
          0 0 0 -s c 0
          0 0 0  0 0 1 ]

    F = T'*F

    # generate a map
    keys = [:ux, :uy, :rz]
    map  = Int[ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]

    return F, map
end


function elem_stiffness(elem::Element{MechFrame})
    C  = get_coords(elem)
    ℓ  = norm(C[2,:]-C[1,:])
    ℓ2 = ℓ*ℓ
    ℓ3 = ℓ*ℓ*ℓ
    etype = elem.etype
    EA = etype.E*etype.A
    EI = etype.E*etype.I

    K0 = [ EA/ℓ     0         0         -EA/ℓ    0         0
           0       12*EI/ℓ3   6*EI/ℓ2    0     -12*EI/ℓ3   6*EI/ℓ2
           0        6*EI/ℓ2   4*EI/ℓ     0      -6*EI/ℓ2   2*EI/ℓ
          -EA/ℓ     0          0         EA/ℓ     0        0
           0      -12*EI/ℓ3  -6*EI/ℓ2    0      12*EI/ℓ3  -6*EI/ℓ2
           0        6*EI/ℓ2   2*EI/ℓ     0      -6*EI/ℓ2   4*EI/ℓ  ]


    # Rotation matrix
    c = (C[2,1] - C[1,1])/ℓ
    s = (C[2,2] - C[1,2])/ℓ

    T = [  c s 0  0 0 0
          -s c 0  0 0 0
           0 0 1  0 0 0
           0 0 0  c s 0
           0 0 0 -s c 0
           0 0 0  0 0 1 ]

    map = elem_map(elem)
    return T'*K0*T, map, map
end


function elem_mass(elem::Element{MechFrame})
    C  = get_coords(elem)
    ℓ  = norm(C[2,:]-C[1,:])
    ℓ2 = ℓ*ℓ
    mat = elem.cmodel


    M0 = mat.ρ*ℓ/420.0*[ 140   0      0      70    0      0
                         0     156    22*ℓ   0     54    -ℓ3*ℓ
                         0     22*ℓ   4*ℓ2   0     ℓ3*ℓ  ℓ3*ℓ2
                         70    0      0      140   0      0
                         0     54     ℓ3*ℓ   0     156   -22*ℓ
                         0    -ℓ3*ℓ  ℓ3*ℓ2   0    -22*ℓ   4*ℓ2 ]

    # Rotation matrix
    c = (C[2,1] - C[1,1])/ℓ
    s = (C[2,2] - C[1,2])/ℓ
    T = [  c s 0  0 0 0
          -s c 0  0 0 0
           0 0 1  0 0 0
           0 0 0  c s 0
           0 0 0 -s c 0
           0 0 0  0 0 1 ]

    map = elem_map(elem)
    return T'*M0*T, map, map
end


# function update_elem!(elem::Element{MechFrame}, U::Vector{Float64}, Δt::Float64)
#     K, map, map = elem_stiffness(elem)
#     dU = U[map]
#     dF = K*dU
#     return dF, map, success()
# end

function elem_internal_forces(elem::Element{MechFrame}, ΔU::Vector{Float64}=Float64[], Δt::Float64=0.0)
    K, map, map = elem_stiffness(elem)
    update = !isempty(ΔU)
    if update
        ΔF = K*ΔU
    else
        ΔF = zeros(length(map)) # TODO: use ips
    end
    return ΔF, map, success()
end
