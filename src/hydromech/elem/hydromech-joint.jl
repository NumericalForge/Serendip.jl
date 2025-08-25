# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export HMJoint

struct HMJointProps<:ElemProperties
    function HMJointProps(; params...)
        return new()
    end    
end


mutable struct HMJoint<:Hydromech
    ctx::Context
    id    ::Int
    shape ::CellShape
    nodes ::Array{Node,1}
    ips   ::Array{Ip,1}
    tag   ::String
    mat   ::Material
    props ::HMJointProps
    active::Bool
    couplings::Array{Element,1}

    function HMJoint(props=HMJointProps())
        return new()
    end
end

# Return the shape family that works with this element
compat_role(::Type{HMJoint}) = JOINTCELL
compat_elem_props(::Type{HMJoint}) = HMJointProps



function elem_config_dofs(elem::HMJoint)
    nnodes = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints
    nlnodes = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            add_dof(node, :uw, :fw)  
        end
        if  i<=(dnlnodes)
            add_dof(node, :ux, :fx)
            add_dof(node, :uy, :fy)
            elem.ctx.ndim==3 && add_dof(node, :uz, :fz)         
        end
    end
end


function elem_init(elem::HMJoint)
    # Get linked elements
    e1 = elem.couplings[1]
    e2 = elem.couplings[2]

    # Volume from first linked element
    V1 = 0.0
    C1 = get_coords(e1)

    for ip in e1.ips
        dNdR = e1.shape.deriv(ip.R)
        J    = C1'*dNdR
        detJ = det(J)
        V1  += detJ*ip.w
    end

    # Volume from second linked element
    V2 = 0.0
    C2 = get_coords(e2)
    for ip in e2.ips
        dNdR = e2.shape.deriv(ip.R)
        J    = C2'*dNdR
        detJ = det(J)
        V2  += detJ*ip.w
    end

    # Area of joint element
    A = 0.0
    C = get_coords(elem)
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints 
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    C = C[1:nlnodes, :]
    fshape = elem.shape.facet_shape

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = fshape.deriv(ip.R)
        J    = C'*dNdR
        detJ = norm2(J)
        A += detJ*ip.w
    end

    # Calculate and save h at joint element's integration points
    h = (V1+V2)/(2.0*A)
    for ip in elem.ips
        ip.state.h = h
    end

    # Setting initial crack openning if available
    if hasfield(typeof(elem.pmodel), :w)
        if elem.pmodel.w > 0.0
            for ip in elem.ips
                ip.state.w[1] = elem.pmodel.w
                # ip.state.up = elem.pmodel.wc
                ip.state.up = 1e-10 # a value different from zero
                ip.state.Δλ = 1.0
            end 
        end
    end
end


function elem_stiffness(elem::HMJoint)
    ndim     = elem.ctx.ndim
    th       = elem.ctx.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints 
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape
    C        = get_coords(elem)[1:nlnodes,:]

    J        = Array{Float64}(undef, ndim, ndim-1)
    NN       = zeros(ndim, dnlnodes*ndim)
    Bu       = zeros(ndim, dnlnodes*ndim)
    DBu      = zeros(ndim, dnlnodes*ndim)
    K        = zeros(dnlnodes*ndim, dnlnodes*ndim)

    for ip in elem.ips
        # compute shape Jacobian
        N    = fshape.func(ip.R)
        dNdR = fshape.deriv(ip.R)

        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute Bu matrix
        T   = matrixT(J)
        NN .= 0.0  # NN = [ -N[]  N[] ]

        for i in 1:nlnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, nlnodes*ndim + (i-1)*ndim + dof] = N[i]
            end
        end

        @mul Bu = T*NN

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.pmodel, ip.state)
        @mul DBu = D*Bu
        @mul K  += coef*Bu'*DBu
    end
    
    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ node.dofdict[key].eq_id for node in elem.nodes[1:dnlnodes] for key in keys ]

    return K, map, map
end


function elem_coupling_matrix(elem::HMJoint) 
    ndim     = elem.ctx.ndim
    th       = elem.ctx.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape
    C        = get_coords(elem)[1:nlnodes,:]

    J        = Array{Float64}(undef, ndim, ndim-1)
    NN       = zeros(ndim, dnlnodes*ndim)
    Bu       = zeros(ndim, dnlnodes*ndim)
    mf       = [1.0, 0.0, 0.0][1:ndim]
    mfNf     = zeros(ndim, 3*nbsnodes)
    Cup      = zeros(dnlnodes*ndim, 3*nbsnodes) # u-p coupling matrix

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = fshape.deriv(ip.R)

        @mul J = C'*dNdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Np vector
        Np = elem.shape.base_shape.func(ip.R)
        N0 = 0*Np
        Nf = [N0' N0' Np']
        
        # compute Bu matrix
        N  = fshape.func(ip.R)
        T   = matrixT(J)
        NN .= 0.0

        for i in 1:nlnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof               ] = -N[i]
                NN[dof, nlnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        @mul Bu = T*NN
        # compute Cup
        coef = detJ*ip.w*th  
        mfNf = mf*Nf
        Cup -= coef*Bu'*mfNf
    end

    # map
    keys = (:ux, :uy, :uz)[1:ndim]
    map_u = [ node.dofdict[key].eq_id for node in elem.nodes[1:dnlnodes] for key in keys ]
    map_w = [ node.dofdict[:uw].eq_id for node in nodes_p ]

    return Cup, map_u, map_w
end


function elem_conductivity_matrix(elem::HMJoint)
    ndim     = elem.ctx.ndim
    th       = elem.ctx.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape

    C        = get_coords(elem)[1:nbsnodes,:]
    Cl       = zeros(nbsnodes, ndim-1)
    J        = Array{Float64}(undef, ndim, ndim-1)
    Jl       = zeros(ndim-1, ndim-1)
    Bp       = zeros(ndim-1, nbsnodes)

    H        = zeros(3*nbsnodes, 3*nbsnodes)

    nodes_p  = []
    
    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end

    for ip in elem.ips

        # compute shape Jacobian
        dNpdR = elem.shape.base_shape.deriv(ip.R)

        @mul J = C'*dNpdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Bp matrix  
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # new coordinate nodes

        @mul Jl = Cl'*dNpdR
        dNdX = dNpdR*inv(Jl)
        Bp = dNdX'
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 

        # compute NN matrix
        Np = elem.shape.base_shape.func(ip.R)  
        N0 = 0*Np
        Nb = [-Np' N0' Np']
        Nt = [N0' -Np' Np']

        # compute H
        coef  = detJ*ip.w*th*elem.pmodel.kt
        H -= coef*Nb'*Nb
        H -= coef*Nt'*Nt

         # compute crack aperture
        # if elem.pmodel.w == 0.0
            if ip.state.up == 0.0 || ip.state.w[1] <= 0.0  
                w = 0.0
            else
                w = ip.state.w[1]
            end
        # else
        #     if elem.pmodel.w >= ip.state.w[1]
        #         w = elem.pmodel.w
        #     else 
        #         w = ip.state.w[1]
        #     end
        # end    

        coef = detJ*ip.w*th*(w^3)/(12*elem.pmodel.η) 
        H -= coef*Bf'*Bf
    end
    
    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p  ]

    return H, map, map, nodes_p
end


function elem_compressibility_matrix(elem::HMJoint)
    ndim     = elem.ctx.ndim
    th       = elem.ctx.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape
    C        = get_coords(elem)[1:nbsnodes,:]

    J   = Array{Float64}(undef, ndim, ndim-1)
    Cpp = zeros(3*nbsnodes, 3*nbsnodes)

    nodes_p  = []
    
    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end 

    for ip in elem.ips
        # compute shape Jacobian
        dNpdR = elem.shape.base_shape.deriv(ip.R)
        @mul J = C'*dNpdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Np matrix
        Np = elem.shape.base_shape.func(ip.R)  
        N0 = 0*Np
        Nf = [N0' N0' Np']

        # compute crack aperture
        # if elem.pmodel.w == 0.0
            if ip.state.up == 0.0 || ip.state.w[1] <= 0.0  
                w = 0.0
            else
                w = ip.state.w[1]
            end
        # else
        #     if elem.pmodel.w >= ip.state.w[1]
        #         w = elem.pmodel.w
        #     else 
        #         w = ip.state.w[1]
        #     end
        # end    

        # compute Cpp
        coef = detJ*ip.w*elem.pmodel.β*w*th
        Cpp -= coef*Nf'*Nf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p  ]

    return Cpp, map, map
end


function elem_RHS_vector(elem::HMJoint)
    ndim     = elem.ctx.ndim
    th       = elem.ctx.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape
    C        = get_coords(elem)[1:nbsnodes,:]

    J        = Array{Float64}(undef, ndim, ndim-1)
    Cl       = zeros(nbsnodes, ndim-1)
    Jl       = zeros(ndim-1, ndim-1)
    Bp       = zeros(ndim-1, nbsnodes)
    Z        = zeros(ndim) 
    Z[end]   = 1.0
    bf       = zeros(ndim-1) 
    Q        = zeros(3*nbsnodes)

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end 
 
    for ip in elem.ips

        # compute shape Jacobian
        dNpdR = elem.shape.base_shape.deriv(ip.R)

        @mul J = C'*dNpdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Bp matrix
        T    = matrixT(J) #rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  #coordinate of new nodes

        @mul Jl = Cl'*dNpdR
        dNdX = dNpdR*inv(Jl)
        Bp = dNdX'
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 
        
        # compute Q

        # compute crack aperture
        # if elem.pmodel.w == 0.0
            if ip.state.up == 0.0 || ip.state.w[1] <= 0.0 
                w = 0.0
            else
                w = ip.state.w[1]
            end 
        # else
        #     if elem.pmodel.w >= ip.state.w[1]
        #         w = elem.pmodel.w
        #     else 
        #         w = ip.state.w[1]
        #     end
        # end    

        coef = detJ*ip.w*th*(w^3)/(12*elem.pmodel.η)   
        bf = T[(2:end), (1:end)]*Z*elem.ctx.γw
        
        @mul Q += coef*Bf'*bf
    end

    # map
    map = [  node.dofdict[:uw].eq_id for node in nodes_p  ]

    return Q, map
end

#=
function elem_internal_forces(elem::HMJoint, F::Array{Float64,1})
    ndim     = elem.ctx.ndim
    th       = elem.ctx.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end 

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes[1:dnlnodes] for key in keys ]
    map_w  = [ node.dofdict[:uw].eq_id for node in nodes_p ]

    dF     = zeros(dnlnodes*ndim)
    Bu     = zeros(ndim, dnlnodes*ndim)
    dFw    = zeros(3*nbsnodes)
    Bp     = zeros(ndim-1, nbsnodes)     

    J      = Array{Float64}(undef, ndim-1, ndim)
    NN     = zeros(ndim, dnlnodes*ndim)

    C      = get_coords(elem)[1:nbsnodes,:]
    Cl     = zeros(nbsnodes, ndim-1)
    Jl     = zeros(ndim-1, ndim-1)
    mf     = [1.0, 0.0, 0.0][1:ndim]    
    Bpuwf  = zeros(ndim-1)
    Z      = zeros(ndim) 
    Z[end] = 1.0


    for ip in elem.ips
        # compute shape Jacobian
        dNpdR = elem.shape.base_shape.deriv(ip.R)

        @mul J = dNpdR*C 
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # coordinate of new nodes
        @mul Jl = dNpdR*Cl
        Bp = inv(Jl)*dNpdR #dNdX
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 

        # compute bf vector
        bf = T[(2:end), (1:end)]*Z*elem.ctx.γw
        
        # compute Np vector
        Np = elem.shape.base_shape.func(ip.R)
        N0 = 0*Np
        Nb = [-Np' N0' Np']
        Nt = [N0' -Np' Np']
        Nf = [N0' N0'  Np']

        # compute NN matrix
        N    = fshape.func(ip.R)
        NN .= 0.0
        for i in 1:nlnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof               ] = -N[i]
                NN[dof, nlnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        @mul Bu = T*NN

        # internal force 
        uwf  = ip.state.uw[3]
        σ    = ip.state.σ[1:ndim] - mf*uwf # get total stress
        coef = detJ*ip.w*th
        @mul dF += coef*Bu'*σ

        # internal volumes dFw
        w  = ip.state.w[1:ndim]
        coef = detJ*ip.w*th
        mfw = mf'*w
        dFw-= coef*Nf'*mfw 

        coef = detJ*ip.w*elem.pmodel.β*th
        dFw -= coef*Nf'*uwf

        # longitudinal flow
        coef = detJ*ip.w*th  
        S = ip.state.S
        dFw -= coef*Bf'*S
         
        # transverse flow
        D = ip.state.D
        dFw -= coef*Nt'*D[1]
        dFw -= coef*Nb'*D[2]
    end

    F[map_u] += dF
    F[map_w] += dFw
end
=#

function update_elem!(elem::HMJoint, U::Array{Float64,1}, Δt::Float64)
    ndim     = elem.ctx.ndim
    th       = elem.ctx.thickness
    nnodes   = length(elem.nodes)
    nbsnodes = elem.shape.base_shape.npoints
    nlnodes  = Int((nnodes-nbsnodes)/2) 
    dnlnodes = nnodes-nbsnodes
    fshape   = elem.shape.facet_shape

    nodes_p  = []

    for (i, node) in enumerate(elem.nodes)
        if  i<=(nbsnodes) || (nlnodes+nbsnodes)>=i>(nlnodes) || i>(dnlnodes)
            push!(nodes_p, elem.nodes[i])
        end
    end 

    keys   = (:ux, :uy, :uz)[1:ndim]
    map_u  = [ node.dofdict[key].eq_id for node in elem.nodes[1:dnlnodes] for key in keys ]
    map_w  = [ node.dofdict[:uw].eq_id for node in nodes_p ]

    dU     = U[map_u] # nodal displacement increments
    dUw    = U[map_w] # nodal pore-pressure increments
    Uw     = [ node.dofdict[:uw].vals[:uw] for node in nodes_p ] # nodal pore-pressure at step n
    Uw    += dUw # nodal pore-pressure at step n+1

    dF     = zeros(dnlnodes*ndim)
    dFw    = zeros(3*nbsnodes)     

    J      = Array{Float64}(undef, ndim,  ndim-1)
    NN     = zeros(ndim, dnlnodes*ndim)
    Δω     = zeros(ndim)
    Δuw    = zeros(3)
    Bu     = zeros(ndim, dnlnodes*ndim)
    C      = get_coords(elem)[1:nbsnodes,:]
    Cl     = zeros(nbsnodes, ndim-1)
    Jl     = zeros(ndim-1, ndim-1)
    Bp     = zeros(ndim-1, nbsnodes)
    mf     = [1.0, 0.0, 0.0][1:ndim]
    BpUwf  = zeros(ndim-1)
    Z      = zeros(ndim) 
    Z[end] = 1.0

    for ip in elem.ips

        # compute shape Jacobian
        dNpdR = elem.shape.base_shape.deriv(ip.R)

        @mul J = C'*dNpdR
        detJ = norm2(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute Np vector
        Np = elem.shape.base_shape.func(ip.R)
        N0 = 0*Np
        Nb = [-Np' N0' Np']
        Nt = [N0' -Np' Np']
        Nf = [N0' N0'  Np']

        # compute Bp matrix
        T    = matrixT(J) # rotation matrix
        Cl   = C*T[(2:end), (1:end)]'  # coordinate of new nodes

        @mul Jl = Cl'*dNpdR
        dNdX = dNpdR*inv(Jl)
        Bp = dNdX'
        B0 = 0*Bp
        Bf = [B0 B0 Bp] 

        # compute bf vector
        bf = T[(2:end), (1:end)]*Z*elem.ctx.γw
        
        # compute NN matrix
        N    = fshape.func(ip.R)
        NN .= 0.0

        for i in 1:nlnodes
            for dof=1:ndim
                NN[dof, (i-1)*ndim + dof               ] = -N[i]
                NN[dof, nlnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        @mul Bu = T*NN

        # interpolation to the integ. point 
        Δuw  = [Np'*dUw[1:nbsnodes]; Np'*dUw[nbsnodes+1:2*nbsnodes]; Np'*dUw[2*nbsnodes+1:end]]
        G    = [ dot(Nt,Uw); dot(Nb,Uw)]
        BfUw = Bf*Uw + bf

        @mul Δω = Bu*dU
          
        # internal force dF
        Δσ, Vt, L, status = update_state(elem.pmodel, ip.state, Δω, Δuw, G, BfUw, Δt)
        failed(status) && return failure("HMJoint: error in update_elem!", status.message)
        Δσ -= mf*Δuw[3] # get total stress
        coef = detJ*ip.w*th
        @mul dF += coef*Bu'*Δσ

        # internal volumes dFw
        coef = detJ*ip.w*th
        mfΔω = mf'*Δω
        dFw -= coef*Nf'*mfΔω 

        # compute fluid compressibility
        # if elem.pmodel.w == 0.0
            if ip.state.up == 0.0 || ip.state.w[1] <= 0.0 
                w = 0.0
            else
                w = ip.state.w[1]
            end 
        # else
        #     if elem.pmodel.w >= ip.state.w[1]
        #         w = elem.pmodel.w
        #     else 
        #         w = ip.state.w[1]
        #     end
        # end    

        coef = detJ*ip.w*elem.pmodel.β*w*th
        dFw -= coef*Nf'*Δuw[3]

        # longitudinal flow
        coef = Δt*detJ*ip.w*th
        dFw -= coef*Bf'*L

        # transverse flow
        dFw += coef*Nt'*Vt[1]  
        dFw += coef*Nb'*Vt[2] 
    end

    return [dF; dFw], [map_u; map_w], success()
end


function elem_recover_nodal_values(elem::HMJoint)
    nips = length(elem.ips)

    keys = output_keys(elem.pmodel)
    vals = zeros(nips, length(keys))
    for (i,ip) in enumerate(elem.ips)
        dict = state_values(elem.pmodel, ip.state)
        vals[i,:] = [ dict[key] for key in keys ]
    end
    
    node_vals = OrderedDict{Symbol, Array{Float64,1}}()
    E = extrapolator(elem.shape.facet_shape, nips)
    for (i,key) in enumerate(keys)
        V = E*vals[:,i]
        node_vals[key] = [ V; V; V ]
    end

    return node_vals
end
