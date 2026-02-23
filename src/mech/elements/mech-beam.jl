export MechBeam

struct MechBeam<:MechFormulation
    b::Float64
    h::Float64
    r::Float64
    αs::Float64
    ρ::Float64
    γ::Float64

    function MechBeam(;r=0.0, b=0.0, h=0.0, A=0.0, rho=0.0, gamma=0.0)
        if r>0 || A>0 # circular section
            if A>0
                r = √(A/π)
            end
            αs = 9/10 # shear correction factor for circular sections
        else
            @check b>0.0 "MechBeam: b (width) must be > 0 for rectangular sections"
            @check h>0.0 "MechBeam: h (height) must be > 0 for rectangular sections"
            αs = 5/6 # shear correction factor for rectangular sections
        end

        @check rho >= 0.0
        @check gamma >= 0.0

        return new(b, h, r, αs, rho, gamma)
    end
end


compat_role(::Type{MechBeam}) = :line
embedded_formulation(::Type{MechBeam}) = error("MechBeam: this element cannot be embedded")


mutable struct MechBeamCache <: ElementCache
    Dlmn::Vector{FixedSizeMatrix{Float64}}  # rotation matrices at nodes
end


function elem_init(elem::Element{MechBeam})
    ndim = elem.ctx.ndim
    nnodes = length(elem.nodes)

    nnodes in (2,3) || error("MechBeam: only two-node and three-node beam elements are supported")

    Dlmn = Vector{FixedSizeMatrix{Float64}}(undef, nnodes)
    C = get_coords(elem)
    J = Matrix{Float64}(undef, ndim, 1)

    for i in 1:nnodes
        Ri = elem.shape.nat_coords[i,:]
        dNdR = elem.shape.deriv(Ri)
        @mul J = C'*dNdR
        L = zeros(ndim,ndim)
        set_rot_x_xp(elem, J, L)
        Dlmn[i] = FixedSizeMatrix(L)
    end
    elem.cache = MechBeamCache(Dlmn)

    # set the value of αs at integration points
    for ip in elem.ips
        ip.state.αs = elem.etype.αs
    end

    return nothing
end


function set_quadrature(elem::Element{MechBeam}, n::Int=0; state::NamedTuple=NamedTuple())
    ndim = elem.ctx.ndim

    if elem.etype.αs == 9/10 # circular section
        if ndim==3
            if n in (0,8)
                nl, nj, nk = 2, 2, 2
            elseif n==12
                nl, nj, nk = 3, 2, 2
            else 
                error("MechBeam: unsupported number of integration points for circular sections in 3D")
            end
        else
            if n in (0,4)
                nl, nj, nk = 2, 2, 1
            elseif n==6
                nl, nj, nk = 3, 2, 1
            else
                error("MechBeam: unsupported number of integration points for circular sections in 2D")
            end
        end

        # longitudinal
        ipL = get_ip_coords(LIN2, nl) 
        
        # transversal quadrature for circular sections
        if ndim==3
            ipT = [ QPoint( -0.5, -0.5,  0.0,  pi/4 ),
                    QPoint(  0.5, -0.5,  0.0,  pi/4 ),
                    QPoint(  0.5,  0.5,  0.0,  pi/4 ),
                    QPoint( -0.5,  0.5,  0.0,  pi/4 ) ]
        else
            ipT = [ QPoint( -0.5,  0.0,  0.0,  pi/4 ),
                    QPoint(  0.5,  0.0,  0.0,  pi/4 ) ]
        end

        nt = length(ipT)
        resize!(elem.ips, nl*nt)
        
        for i in 1:nl
            for j in 1:length(ipT)
                R = [ ipL[i].coord[1], ipT[j].coord[1], ipT[j].coord[2] ]
                w = ipL[i].w*ipT[j].w
                m = (i-1)*nt + j

                ipstate = compat_state_type(typeof(elem.cmodel), typeof(elem.etype))(elem.ctx; state...)
                elem.ips[m] = Ip(R, w, elem, ipstate)
            end
        end

    else # rectangular section
        if ndim==3
            if n in (0,8)
                nl, nj, nk = 2, 2, 2
            elseif n==2
                nl, nj, nk = 2, 1, 1
            elseif n==3
                nl, nj, nk = 3, 1, 1
            elseif n==12
                nl, nj, nk = 3, 2, 2
            else
                error("MechBeam: unsupported number of integration points for rectangular sections in 3D")
            end
        else
            if n in (0,4)
                nl, nj, nk = 2, 2, 1
            elseif n==2
                nl, nj, nk = 2, 1, 1
            elseif n==3
                nl, nj, nk = 3, 1, 1
            elseif n==6
                nl, nj, nk = 3, 2, 1
            else
                error("MechBeam: unsupported number of integration points for rectangular sections in 2D")
            end
        end

        ipL = get_ip_coords(LIN2, nl) # longitudinal
        ipT = get_ip_coords(LIN2, nj) # transversal

        resize!(elem.ips, nl*nj*nk)
        for i in 1:nl
            for j in 1:nj
                for k in 1:nk
                    if ndim==2
                        R = [ ipL[i].coord[1], ipT[j].coord[1], 0.0 ]
                        w = ipL[i].w*ipT[j].w
                    else
                        R = [ ipL[i].coord[1], ipT[j].coord[1], ipT[k].coord[1] ]
                        w = ipL[i].w*ipT[j].w*ipT[k].w
                    end
                    m = (i-1)*nj*nk + (j-1)*nk + k
                    ipstate = compat_state_type(typeof(elem.cmodel), typeof(elem.etype))(elem.ctx; state...)
                    elem.ips[m] = Ip(R, w, elem, ipstate)
                end
            end
        end

    end

    # finding ips global coordinates
    C     = get_coords(elem)
    shape = elem.shape

    s = 0.0
    for ip in elem.ips
        # @show ip.w ip.R
        s+= ip.w
    end
    # @show s

    for ip in elem.ips
        R = [ ip.R[1], 0.0, 0.0 ]
        # R = ip.R
        N = shape.func(R)
        ip.coord = C'*N
    end

end


function distributed_bc(elem::Element{MechBeam}, facet::Cell, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
end


function body_c(elem::Element{MechBeam}, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, key, val)
end


function elem_config_dofs(elem::Element{MechBeam})
    ndim = elem.ctx.ndim
    ndim in (2,3) || error("MechBeam: Beam elements do not work in $(ndim)d analyses")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :rz, :mz)
        if ndim==3
            add_dof(node, :uz, :fz)
            add_dof(node, :rx, :mx)
            add_dof(node, :ry, :my)
        end
    end
end


function elem_map(elem::Element{MechBeam})
    if elem.ctx.ndim==2
        keys =(:ux, :uy, :rz)
    else
        keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    end
    return [ get_dof(node,key).eq_id for node in elem.nodes for key in keys ]
end


dof_map(elem::Element{MechBeam}) = elem_map(elem)


# Rotation Matrix
function set_rot_x_xp(elem::Element{MechBeam}, J::Matx, R::Matx)
    ndim = elem.ctx.ndim
    V1 = normalize(vec(J))
    V1 = round.(V1, digits=14)

    if ndim==2
        V2 = [ -V1[2], V1[1] ]
        R[1,:] .= V1
        R[2,:] .= V2
    else
        
        if V1[1] == 0.0 && V1[2] == 0.0 # aligned with Z axis
            V2 = [ 0.0, 1.0, 0.0 ]
            V3 = [ -1.0, 0.0, 0.0 ]
        else
            Z = [ 0.0, 0.0, 1.0 ]
            V2 = cross(V1, Z)
            V3 = cross(V1, V2)
            normalize!(V2)
            normalize!(V3)
        end

        R[1,:] .= V1
        R[2,:] .= V2
        R[3,:] .= V3
    end

end


function get_ry_rz(elem::Element{MechBeam})
    ndim = elem.ctx.ndim
    if elem.etype.r > 0.0 # circular section
        r_y = elem.etype.r
        r_z = elem.etype.r
    elseif ndim==2 # height is in the y direction
        r_z = elem.etype.b/2
        r_y = elem.etype.h/2
    else # height is in the z direction
        r_z = elem.etype.h/2
        r_y = elem.etype.b/2
    end
    return r_y, r_z
end


function setB(elem::Element{MechBeam}, ip::Ip, L::Matx, N::Vect, dNdX::Matx, Rθ::Matx, Bil::Matx, Bi::Matx, B::Matx)
    ndim = elem.ctx.ndim
    ndof = ndim==2 ? 3 : 6
    nnodes = size(dNdX,1)
    Bil .= 0.0

    r_y, r_z = get_ry_rz(elem)

    # Note that matrix B is designed to work with tensors in Mandel's notation
    if ndim==2
        Rθ[3,3] = 1.0
        for i in 1:nnodes
            η = ip.R[2]
            Rθ[1:2,1:2] .= L

            # MITC: Mixed Interpolation of Tensorial Components Bathe (1985)
            Ni = nnodes>2 ? N[i] : 0.5 # MITC projectrion to reduce shear locking in two-node beams
            dNdx = dNdX[i]

            Bil[1,1] = dNdx;                           Bil[1,3] = -dNdx*η*r_y
                              Bil[3,2] = SR2/2*dNdx;   Bil[3,3] = -SR2/2*Ni

            c = (i-1)*ndof
            @mul Bi = Bil*Rθ
            B[:, c+1:c+ndof] .= Bi
        end
    else
        for i in 1:nnodes
            η = ip.R[2]
            ζ = ip.R[3]
            Rθ[1:3,1:3] .= L
            Rθ[4:6,4:6] .= elem.cache.Dlmn[i]
            # Rθ[1:3,1:3] .= elem.cache.Dlmn[i]

            Ni = nnodes>2 ? N[i] : 0.5 # MITC projectrion to reduce shear locking in two-node beams

            dNdx = dNdX[i]

            Bil[1,1] = dNdx;                                                                             Bil[1,5] = dNdx*ζ*r_z;  Bil[1,6] = -dNdx*η*r_y
                                                   Bil[2,3] = SR2/2*dNdx;  Bil[2,4] = SR2/2*dNdx*η*r_y;  Bil[2,5] = SR2/2*Ni
                            Bil[3,2] = SR2/2*dNdx;                         Bil[3,4] = -SR2/2*dNdx*ζ*r_z;                        Bil[3,6] = -SR2/2*Ni

            c = (i-1)*ndof
            @mul Bi = Bil*Rθ
            B[:, c+1:c+ndof] .= Bi
        end
    end
end


function elem_stiffness(elem::Element{MechBeam})
    ndim = elem.ctx.ndim
    nnodes = length(elem.nodes)
    ndof = ndim==2 ? 3 : 6
    nstr = 3

    r_y, r_z = get_ry_rz(elem)

    C   = get_coords(elem)
    K   = zeros(ndof*nnodes, ndof*nnodes)
    B   = zeros(nstr, ndof*nnodes)
    L   = zeros(ndim, ndim)
    Rθ  = zeros(ndof, ndof)
    Bil = zeros(nstr, ndof)
    Bi  = zeros(nstr, ndof)

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J1D  = C'*dNdR
        set_rot_x_xp(elem, J1D, L)
        dx′dξ = norm(J1D)
        dNdX′ = dNdR*inv(dx′dξ)
        D     = calcD(elem.cmodel, ip.state)

        setB(elem, ip, L, N, dNdX′, Rθ, Bil, Bi, B)

        detJ′ = dx′dξ*r_z*r_y
        w2d = ndim==2 ? 2.0 : 1.0
        coef = detJ′*ip.w*w2d

        K += coef*B'*D*B
    end

    map = elem_map(elem)
    return K, map, map
end


function elem_internal_forces(elem::Element{MechBeam}, ΔU::Vector{Float64}=Float64[], dt::Float64=0.0)
    ndim = elem.ctx.ndim
    nnodes = length(elem.nodes)
    ndof = ndim==2 ? 3 : 6
    nstr = 3

    r_y, r_z = get_ry_rz(elem)

    C   = get_coords(elem)
    B   = zeros(nstr, ndof*nnodes)
    L   = zeros(ndim,ndim)
    Rθ  = zeros(ndof,ndof)
    Bil = zeros(nstr,ndof)
    Bi  = zeros(nstr,ndof)

    map = elem_map(elem)
    ΔF  = zeros(ndof*nnodes)

    update = !isempty(ΔU)
    if update
        Δε = zeros(nstr)
    end

    for ip in elem.ips
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J1D  = C'*dNdR
        set_rot_x_xp(elem, J1D, L)
        dx′dξ = norm(J1D)
        dNdX′ = dNdR*inv(dx′dξ)
        setB(elem, ip, L, N, dNdX′, Rθ, Bil, Bi, B)
        if update
            @mul Δε = B*ΔU
            Δσ, status = update_state(elem.cmodel, ip.state, ip.cstate, Δε)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        detJ′ = dx′dξ*r_z*r_y
        w2d = ndim==2 ? 2.0 : 1.0
        coef = detJ′*ip.w*w2d
        ΔF += coef*B'*Δσ
    end

    return ΔF, map, success()
end


function elem_vals(elem::Element{MechBeam})
    # get ip average values
    ipvals = [ state_values(elem.cmodel, ip.state) for ip in elem.ips ]
    # sum  = merge(+, ipvals... )
    # nips = length(elem.ips)
    merger(x,y) = abs(x) > abs(y) ? x : y
    merged  = merge(merger, ipvals... )
    vals = OrderedDict( k=>v for (k,v) in merged)
    return vals
end


function elem_recover_nodal_values(elem::Element{MechBeam})
    ndim = elem.ctx.ndim
    r_y, r_z = get_ry_rz(elem)
    ndof = ndim==2 ? 3 : 6

    nnodes = length(elem.nodes)
    C = get_coords(elem) # global coordinates
    Ξ = elem.shape.nat_coords # natural coordinates

    # get local displacements from global
    if ndim==2
        keys =(:ux, :uy, :rz)
    else
        keys =(:ux, :uy, :uz, :rx, :ry, :rz)
    end

    U = [ get_dof(node, key).vals[key] for node in elem.nodes for key in keys ]
    U′ = similar(U)
    Rθ = zeros(ndof,ndof)

    # coefficient matrices to compute moment and shear
    Am = zeros(nnodes, nnodes)
    Av = zeros(nnodes, nnodes)

    for i in 1:nnodes
        if ndim==2
            Rθ[1:2,1:2] .= elem.cache.Dlmn[i]
            Rθ[3,3] = 1.0
        else
            Rθ[1:3,1:3] .= elem.cache.Dlmn[i]
            Rθ[4:6,4:6] .= elem.cache.Dlmn[i]
        end

        U′[(i-1)*ndof+1:i*ndof] = Rθ*U[(i-1)*ndof+1:i*ndof]
        ξ = Ξ[i]

        # first derivatives and jacobian
        dNdξ = elem.shape.deriv([ξ])
        jac  = norm(C'*dNdξ)

        # second derivatives
        d2Ndξ2 = elem.shape.deriv2([ξ])
        d2xdξ2 = dot(C'*d2Ndξ2, elem.cache.Dlmn[i][1,:])
        d2ξNdx2 = -d2xdξ2/jac^3

        for j in 1:nnodes
            dNdx = dNdξ[j]/jac
            Am[i,j] = dNdx
            d2Ndx2 = d2Ndξ2[j]*(1/jac)^2 + dNdξ[j]*d2ξNdx2
            Av[i,j] = d2Ndx2
        end

    end

    E = elem.cmodel.E
    if ndim==2
        θZ = U′[3:ndof:ndof*nnodes]
        # Izy = h*b^3/12
        Izy = 16*r_z*r_y^3/12

        # Bending moment
        Mxy = (Am*θZ).*(E*Izy)
        # Shear (negative to match convention)
        Vxy = -(Av*θZ).*(E*Izy)
        Vxy = ones(nnodes)*mean(Vxy)

        return OrderedDict(
            :Mx´y´ => Mxy,
            :Vx´y´ => Vxy
        )
    else
        θY = U′[5:ndof:ndof*nnodes-1]
        θZ = U′[6:ndof:ndof*nnodes]

        Iyz = 16*r_y*r_z^3/12
        Izy = 16*r_z*r_y^3/12

        # Bending moment
        Mxz = (Am*θY).*(E*Iyz)
        Mxy = (Am*θZ).*(E*Izy)
        # Shear (negative to match convention)
        Vxz = -(Av*θY).*(E*Iyz)
        Vxy = -(Av*θZ).*(E*Izy)
        Vxz = ones(nnodes)*mean(Vxz)
        Vxy = ones(nnodes)*mean(Vxy)

        return OrderedDict(
            :Mx´z´ => Mxz,
            :Mx´y´ => Mxy,
            :Vx´z´ => Vxz,
            :Vx´y´ => Vxy
        )
    end
end
