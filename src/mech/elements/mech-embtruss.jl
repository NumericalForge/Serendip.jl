# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

struct MechEmbBar<:MechFormulation
    A::Float64

    function MechEmbBar(;A=Nan)
        @check A > 0.0
        return new(A)
    end
end


compat_role(::Type{MechEmbBar}) = :line


mutable struct MechEmbBarCache <: ElementCache
    # NN matrix to map host displacements to embedded displacements
    NN::FixedSizeMatrix{Float64}
end


function elem_config_dofs(elem::Element{MechEmbBar})
    # The nodes of a MechEmbBar element are auxiliary nodes.
    # They are not used directly in the analysis but are used
    # to define the embedded element and for output.

    for node in elem.nodes
        add_dof(node, :ux, :fx)
        elem.ctx.ndim>=2 && add_dof(node, :uy, :fy)
        elem.ctx.ndim==3 && add_dof(node, :uz, :fz)
        node.aux = true # mark as auxiliary
    end
end


function elem_map(elem::Element{MechEmbBar})
    ndim = elem.ctx.ndim
    keys = (:ux, :uy, :uz)[1:ndim]
    solid = elem.couplings[1]
    return [ get_dof(node,key).eq_id for node in solid.nodes for key in keys ]
end


dof_map(elem::Element{MechEmbBar}) = elem_map(elem)


function _mountNN(elem::Element{MechEmbBar})
    ndim = elem.ctx.ndim
    solid = elem.couplings[1]
    n  = length(solid.nodes)
    m  = length(elem.nodes)
    NN = FixedSizeMatrix(zeros(ndim*m, ndim*n))
    Cs = get_coords(solid)

    for j in 1:m
        R = inverse_map(solid.shape, Cs, elem.nodes[j].coord)
        N = solid.shape.func(R)
        for i in 1:n
            for k in 1:ndim
                # NN[(i-1)*ndim+k, (j-1)*ndim+k] = N[i]
                NN[(j-1)*ndim+k, (i-1)*ndim+k] = N[i]
            end
        end
    end
    return NN
end


function elem_init(elem::Element{MechEmbBar})
    for node in elem.nodes
        get_dof(node, :ux).eq_id = -1
        get_dof(node, :uy).eq_id = -1
        elem.ctx.ndim==3 && (get_dof(node, :uz).eq_id = -1)
        # node.dofdict[:ux].eq_id = -1
        # node.dofdict[:uy].eq_id = -1
        # node.dofdict[:uz].eq_id = -1
    end

    elem.cache = MechEmbBarCache( _mountNN(elem) )
    # elem.cacheM = Array{Matrix}(undef, 1)
    # elem.cacheM = Vector{FixedSizeMatrix{Float64}}(undef, 1)

    # elem.cache.NN = _mountNN(elem)
    return nothing
end


function elem_displacements(elem::Element{MechEmbBar})
    ndim    = elem.ctx.ndim
    NN      = elem.cache.NN
    keys    = (:ux, :uy, :uz)[1:ndim]
    Uhost   = [ node.dofdict[key].vals[key] for node in elem.couplings[1].nodes for key in keys ]
    Ubar    = NN*Uhost
    nodemap = [ node.id for node in elem.nodes for key in keys ]
    dimmap  = [ i for node in elem.nodes for i in 1:ndim ]

    return Ubar, nodemap, dimmap
end


function elem_stiffness(elem::Element{MechEmbBar})
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    A = elem.etype.A
    C = get_coords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        E    = calcD(elem.cmodel, ip.state)
        coef = E*A*detJ*ip.w
        @mul K += coef*B'*B
    end

    NN = elem.cache.NN
    map = elem_map(elem)
    # @show "stiffness"
    return NN'*K*NN, map, map
end


function elem_internal_forces(elem::Element{MechEmbBar}, ΔU::Vector{Float64}=Float64[], dt::Float64=0.0)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    A      = elem.etype.A
    NN     = elem.cache.NN

    map = elem_map(elem)
    update = !isempty(ΔU)
    if update
        ΔUbar = NN*ΔU
    end

    ΔF = zeros(nnodes*ndim)
    C = get_coords(elem)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(undef, ndim, 1)

    for ip in elem.ips
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm(J)

        # mount B
        B .= 0.0
        for i in 1:nnodes
            for j in 1:ndim
                B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
            end
        end

        if update
            Δε = (B*ΔUbar)[1]
            Δσ, status = update_state(elem.cmodel, ip.state, ip.cstate, Δε)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        # σ = ip.state.σ
        coef = A*detJ*ip.w
        ΔF .+= coef*Δσ*vec(B')
    end

    return NN'*ΔF, map, success()
end


# function update_elem!(elem::Element{MechEmbBar}, U::Vector{Float64}, Δt::Float64)
#     ndim   = elem.ctx.ndim
#     nnodes = length(elem.nodes)
#     A      = elem.etype.A
#     NN     = elem.cache.NN

#     map = elem_map(elem)
#     dU  = U[map]
#     dUbar = NN*dU

#     dF = zeros(nnodes*ndim)
#     C  = get_coords(elem)
#     B  = zeros(1, nnodes*ndim)
#     J  = Array{Float64}(undef, ndim, 1)
#     for ip in elem.ips
#         dNdR = elem.shape.deriv(ip.R)
#         @mul J = C'*dNdR
#         detJ = norm(J)

#         # mount B
#         B .= 0.0
#         for i in 1:nnodes
#             for j in 1:ndim
#                 B[1,j+(i-1)*ndim] = dNdR[i,1]*J[j]/detJ^2.0
#             end
#         end

#         dε = (B*dUbar)[1]
#         dσ, _ = update_state(elem.cmodel, ip.state, dε)
#         coef = A*detJ*ip.w
#         @mul dF += coef*B'*dσ
#     end

#     # update nodal displacements
#     keys    = (:ux, :uy, :uz)[1:ndim]
#     Uhost   = [ get_dof(node, key).vals[key] for node in elem.couplings[1].nodes for key in keys ]
#     Ubar    = NN*Uhost + dUbar
#     for (i,node) in enumerate(elem.nodes)
#         for (j,key) in enumerate(keys)
#             get_dof(node, key).vals[key] = Ubar[(i-1)*ndim+j]
#         end
#     end

#     return NN'*dF, map, success()
# end


function elem_vals(elem::Element{MechEmbBar})
    # get area and average stress and axial force
    vals = OrderedDict(:A => elem.etype.A )
    σx´ = [ state_values(elem.cmodel, ip.state)[:σx´] for ip in elem.ips ]
    _, idx = findmax(abs, σx´)
    max_σx´ = σx´[idx]
    vals[:σx´] = max_σx´
    vals[:fx´] = elem.etype.A*max_σx´
    return vals
end


# function post_process(elem::Element{MechEmbBar})
#     ndim    = elem.ctx.ndim
#     NN      = elem.cache.NN

#     # update displacements
#     keys    = (:ux, :uy, :uz)[1:ndim]
#     Uhost   = [ node.vals[key] for node in elem.couplings[1].nodes for key in keys ]
#     Ubar    = NN*Uhost

#     for (i,node) in enumerate(elem.nodes)
#         for (j,key) in enumerate(keys)
#             node.vals[key] = Ubar[(i-1)*ndim+j]
#         end
#     end

# end
