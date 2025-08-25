# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechShell

mutable struct MechShell<:MechFormulation
    th::Float64
    κ::Float64
    ρ::Float64
    γ::Float64

    function MechShell(;thickness::Float64=0.0, kappa::Float64=1e-8, rho=0.0, gamma=0.0)
        return new(thickness, kappa, rho, gamma)
    end
end

compat_role(::Type{MechShell}) = :surface


function elem_init(elem::Element{MechShell})
    # check element dimension
    elem.shape.ndim==2 || throw(SerendipException("MechShell: Invalid element shape. Got $(elem.shape.name)"))

    # Compute nodal rotation matrices
    nnodes = length(elem.nodes)
    elem.cacheM = Vector{FixedSizeMatrix{Float64}}(undef, nnodes)
    C = get_coords(elem)

    for i in 1:nnodes
        Ri = elem.shape.nat_coords[i,:]
        dNdR = elem.shape.deriv(Ri)
        J = C'*dNdR

        V1 = J[:,1]
        V2 = J[:,2]
        V3 = cross(V1, V2)
        V2 = cross(V3, V1)
        normalize!(V1)
        normalize!(V2)
        normalize!(V3)

        elem.cacheM[i] = FixedSizeMatrix([V1'; V2'; V3'])
    end

    return nothing
end


function set_quadrature(elem::Element{MechShell}, n::Int=0)
    # Set integration points
    if n in (8, 18)
        n = div(n,2)
    end
    ip2d = get_ip_coords(elem.shape, n)
    ip1d = get_ip_coords(LIN2, 2)
    n = size(ip2d,1)

    resize!(elem.ips, 2*n)
    for k in 1:2
        for i in 1:n
            R = [ ip2d[i].coord[1:2]; ip1d[k].coord[1] ]
            w = ip2d[i].w*ip1d[k].w
            j = (k-1)*n + i
            elem.ips[j] = Ip(R, w)
            elem.ips[j].id = j
            elem.ips[j].state = compat_state_type(typeof(elem.pmodel), MechShell, elem.ctx)(elem.ctx)
            elem.ips[j].owner = elem
        end
    end

    # finding ips global coordinates
    C     = get_coords(elem)
    shape = elem.shape

    for ip in elem.ips
        R = [ ip.R[1:2]; 0.0 ]
        N = shape.func(R)
        ip.coord = C'*N

        dNdR = elem.shape.deriv(ip.R) # 3xn
        J = C'*dNdR
        No = normalize(cross(J[:,1], J[:,2]))
        ip.coord += elem.eform.th/2*ip.R[3]*No
    end

end


function distributed_bc(elem::Element{MechShell}, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return mech_boundary_forces(elem, facet, t, key, val)
end


function body_c(elem::Element{MechShell}, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_shell_body_forces(elem, 0.0, key, val)
end


function elem_config_dofs(elem::Element{MechShell})
    ndim = elem.ctx.ndim
    ndim==3 || error("MechShell: Shell elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
        add_dof(node, :rx, :mx)
        add_dof(node, :ry, :my)
        add_dof(node, :rz, :mz)
    end
end


# Rotation Matrix
function set_rot_x_xp(elem::Element{MechShell}, J::Matx, R::Matx)
    V1 = J[:,1]
    V2 = J[:,2]
    V3 = cross(V1, V2)
    V2 = cross(V3, V1)

    normalize!(V1)
    normalize!(V2)
    normalize!(V3)

    R[1,:] .= V1
    R[2,:] .= V2
    R[3,:] .= V3
end


function elem_map(elem::Element{MechShell})
    keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    # return [ node.dofdict[key].eq_id for node in elem.nodes for key in keys ]
    return [ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]
end


function setB(elem::Element{MechShell}, ip::Ip, N::Vect, L::Matx, dNdX::Matx, Rθ::Matx, Bil::Matx, Bi::Matx, B::Matx)
    nnodes = size(dNdX,1)
    th = elem.eform.th
    # Note that matrix B is designed to work with tensors in Mandel's notation

    ndof = 6
    for i in 1:nnodes
        ζ = ip.R[3]
        Rθ[1:3,1:3] .= L
        # Rθ[1:3,1:3] .= elem.cacheM[i]
        Rθ[4:5,4:6] .= elem.cacheM[i][1:2,:]

        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        Ni = N[i]

        Bil[1,1] = dNdx;                                                                                Bil[1,5] = dNdx*ζ*th/2
                             Bil[2,2] = dNdy;                           Bil[2,4] = -dNdy*ζ*th/2

                                                  Bil[4,3] = dNdy/SR2;  Bil[4,4] = -1/SR2*Ni
                                                  Bil[5,3] = dNdx/SR2;                                  Bil[5,5] = 1/SR2*Ni
        Bil[6,1] = dNdy/SR2; Bil[6,2] = dNdx/SR2;                       Bil[6,4] = -1/SR2*dNdx*ζ*th/2;  Bil[6,5] = 1/SR2*dNdy*ζ*th/2

        c = (i-1)*ndof
        @mul Bi = Bil*Rθ
        B[:, c+1:c+6] .= Bi
    end
end


function setB_dr(elem::Element{MechShell}, N::Vect, L::Matx, dNdX::Matx, Rθ_dr::Matx, Bil_dr::Matx, Bi_dr::Matx, B_dr::Matx)
    nnodes = size(dNdX,1)

    ndof = 6
    for i in 1:nnodes

        Rθ_dr[1:2,1:3] .= L[1:2,1:3]
        Rθ_dr[3,4:6] .= elem.cacheM[i][3,:]

        dNdx = dNdX[i,1]
        dNdy = dNdX[i,2]
        Ni   = N[i]


        Bil_dr[1,1] = +0.5*dNdy
        Bil_dr[1,2] = -0.5*dNdx
        Bil_dr[1,3] = Ni

        c = (i-1)*ndof
        @mul Bi_dr = Bil_dr*Rθ_dr

        B_dr[:, c+1:c+6] .= Bi_dr
    end
end


function setNN(elem::Element{MechShell}, ip::Ip, N::Vect, NNil::Matx, NNi::Matx, L::Matx, Rθ::Matx, NN::Matx)
    nnodes = length(N)
    ndof = 6
    th = elem.eform.th
    ζ = ip.R[3]

    for i in 1:nnodes

        Rθ[1:3,1:3] .= L
        Rθ[4:5,4:6] .= L[1:2,:]

        NNil[1,1] = N[i]
        NNil[2,2] = N[i]
        NNil[3,3] = N[i]
        NNil[1,5] = th/2*ζ*N[i]
        NNil[2,4] = -th/2*ζ*N[i]

        c = (i-1)*ndof
        @mul NNi = NNil*Rθ

        NN[:, c+1:c+6] .= NNi
    end
end


function elem_stiffness(elem::Element{MechShell})
    nnodes = length(elem.nodes)
    th     = elem.eform.th
    κ      = elem.eform.κ
    ndof   = 6
    nstr   = 6
    C      = get_coords(elem)
    K      = zeros(ndof*nnodes, ndof*nnodes)
    B      = zeros(nstr, ndof*nnodes)
    L      = zeros(3,3)
    Rθ     = zeros(5,ndof)
    Bil    = zeros(nstr,5)
    Bi     = zeros(nstr,ndof)

    B_dr   = zeros(1, ndof*nnodes)
    Bil_dr = zeros(1,3)
    Bi_dr  = zeros(1,ndof)
    Rθ_dr  = zeros(3,ndof)
    m      = div(length(elem.ips), 2) # half the number of integration points

    for (i,ip) in enumerate(elem.ips)
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J2D  = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′   = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)
        dNdR  = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        detJ′ = det(J′)
        @assert detJ′>0

        setB(elem, ip, N, L, dNdX′, Rθ, Bil, Bi, B)

        E  = elem.pmodel.E
        nu = elem.pmodel.ν
        G  = E/(2*(1+nu))

        coef  = detJ′*ip.w
        D     = calcD(elem.pmodel, ip.state, is_shell=true)
        K    += coef*B'*D*B
        # K    += coef*B'*S*D*B

        if i<=m # drilling stiffness (area integration)
            setB_dr(elem, N, L, dNdX′, Rθ_dr, Bil_dr, Bi_dr, B_dr)
            coef = κ*G*norm2(J2D)*th*ip.w
            @mul K += coef*B_dr'*B_dr
        end

    end

    map = elem_map(elem)
    return K, map, map
end


function elem_mass(elem::Element{MechShell})
    nnodes = length(elem.nodes)
    th     = elem.eform.th
    ndof   = 6 #6
    ρ      = elem.eform.ρ
    C      = get_coords(elem)
    M      = zeros(nnodes*ndof, nnodes*ndof)
    L      = zeros(3,3)
    Rθ   = zeros(5,ndof)

    NN     = zeros(3, nnodes*ndof)
    NNil    = zeros(3,5)
    NNi     = zeros(3,ndof)

    for ip in elem.ips
        # compute N matrix
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J2D  = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′   = [ L*J2D [ 0,0,th/2]  ]

        detJ′ = det(J′)
        @assert detJ′>0

        setNN(elem, ip, N, NNil, NNi, L, Rθ, NN)

        # compute M
        coef = ρ*detJ′*ip.w
        @mul M += coef*NN'*NN
    end

    map = elem_map(elem)
    return M, map, map
end


function elem_internal_forces(elem::Element{MechShell}, ΔUg::Vector{Float64}=Float64[], dt::Float64=0.0)
    nnodes = length(elem.nodes)
    th   = elem.eform.th
    ndof = 6

    map = elem_map(elem)
    ΔF  = zeros(ndof*nnodes)

    C = get_coords(elem)
    B = zeros(6, ndof*nnodes)

    L   = zeros(3,3)
    Rθ  = zeros(5,ndof)
    Bil = zeros(6,5)
    Bi  = zeros(6,ndof)

    update = !isempty(ΔUg)
    if update
        ΔU = ΔUg[map]
        Δε = zeros(6)
    end

    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R) # 3xn

        J2D = C'*dNdR
        set_rot_x_xp(elem, J2D, L)
        J′ = [ L*J2D [ 0,0,th/2]  ]
        invJ′ = inv(J′)

        dNdR = [ dNdR zeros(nnodes) ]
        dNdX′ = dNdR*invJ′

        setB(elem, ip, N, L, dNdX′, Rθ, Bil, Bi, B)

        if update
            @mul Δε = B*ΔU
            Δσ, status = update_state(elem.pmodel, ip.state, Δε)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        detJ′ = det(J′)
        coef  = detJ′*ip.w
        ΔF   += coef*B'*Δσ

    end

     return ΔF, map, success()
end

# function update_elem!(elem::Element{MechShell}, U::Array{Float64,1}, dt::Float64)
#     nnodes = length(elem.nodes)
#     th   = elem.eform.th
#     ndof = 6

#     map = elem_map(elem)
#     dU  = U[map]
#     dF  = zeros(length(dU))

#     C = get_coords(elem)
#     B = zeros(6, ndof*nnodes)

#     L    = zeros(3,3)
#     Rθ = zeros(5,ndof)
#     Bil  = zeros(6,5)
#     Bi   = zeros(6,ndof)
#     Δε   = zeros(6)

#     for ip in elem.ips
#         N = elem.shape.func(ip.R)
#         dNdR = elem.shape.deriv(ip.R) # 3xn

#         J2D = C'*dNdR
#         set_rot_x_xp(elem, J2D, L)
#         J′ = [ L*J2D [ 0,0,th/2]  ]
#         invJ′ = inv(J′)

#         dNdR = [ dNdR zeros(nnodes) ]
#         dNdX′ = dNdR*invJ′

#         setB(elem, ip, N, L, dNdX′, Rθ, Bil, Bi, B)
#         Δε = B*dU
#         Δσ, status = update_state(elem.pmodel, ip.state, Δε, is_shell=true)
#         failed(status) && return dF, map, failure("MechShell: Error at integration point $(ip.id)")

#         detJ′ = det(J′)
#         coef  = detJ′*ip.w
#         dF   += coef*B'*Δσ

#     end

#      return dF, map, success()
# end
