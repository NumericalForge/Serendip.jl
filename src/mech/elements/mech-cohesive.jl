# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechCohesive

struct MechCohesive<:MechFormulation
end


# Return the shape family that works with this element
compat_role(::Type{MechCohesive}) = :cohesive


mutable struct MechCohesiveCache <: ElementCache
    mobilized::Bool
end


function elem_init(elem::Element{MechCohesive})
    # Computation of characteristic length 'h' for cohesive elements
    # and set it in the integration point state

    hasfield(typeof(elem.ips[1].state), :h) || return
    
    ndim = elem.ctx.ndim

    # Avg volume of linked elements
    V = 0.0
    for elem in elem.couplings
        V += cell_extent(elem)
    end
    V /= length(elem.couplings)

    # Area of cohesive element
    A = 0.0
    C = get_coords(elem)
    n = div(length(elem.nodes), 2)
    C = C[1:n, :]
    J = fzeros(ndim, ndim-1)

    for ip in elem.ips
        # compute shape Jacobian
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm2(J)
        detJ <= 0 && error("Invalid Jacobian norm for cohesive element")
        A += detJ*ip.w
    end

    # Calculate and save h at cohesive element's integration points
    h = V/A
    for ip in elem.ips
        ip.state.h = h
    end

    # Explicit cohesive elements already have duplicated/opposite-face nodes.
    # Intrinsic/extrinsic inserted interfaces usually start with the same node ids on both faces.
    mobilized = elem.nodes[1].id != elem.nodes[n+1].id
    elem.cache = MechCohesiveCache(mobilized)
end


function elem_stiffness(elem::Element{MechCohesive})
    elem.cache.mobilized == false && return zeros(0,0), Int[], Int[] # return empty arrays if the element is inactive

    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    nstr   = 3

    C = get_coords(elem)[1:hnodes,:]
    B = fzeros(nstr, nnodes*ndim)
    K = fzeros(nnodes*ndim, nnodes*ndim)

    DB = fzeros(nstr, nnodes*ndim)
    J  = fzeros(ndim, ndim-1)
    NN = fzeros(nstr, nnodes*ndim)

    for ip in elem.ips
        if elem.ctx.stress_state==:axisymmetric
            th = 2*pi*ip.coord.x
        end

        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute B matrix
        for i in 1:hnodes
            for dof in 1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end
        
        T = calc_interface_rotation(J)
        @mul B = T'*NN

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.cmodel, ip.state)
        
        @mul DB = D*B
        @mul K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = Int[ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_internal_forces(elem::Element{MechCohesive}, ΔU::Vector{Float64}=Float64[], Δt::Float64=0.0)
    if elem.cache.mobilized == false 
        # recover cohesive stress from neighboring bulk elements
        σ_arr = recover_cohesive_stress(elem)
        for (i, ip) in enumerate(elem.ips)
            ip.state.σ = cap_stress(elem.cmodel, σ_arr[i])
        end

        return zeros(0), Int[], success() # return empty arrays if the element is inactive
    end

    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    map    = dof_map(elem)
    nstr   = 3

    update = !isempty(ΔU)
    if update
        Δω = zeros(nstr)
    end

    ΔF = fzeros(nnodes*ndim)
    C  = get_coords(elem)[1:hnodes,:]
    B  = fzeros(nstr, nnodes*ndim)

    J  = fzeros(ndim, ndim-1)
    NN = fzeros(nstr, nnodes*ndim)

    for ip in elem.ips
        if elem.ctx.stress_state==:axisymmetric
            th = 2*pi*ip.coord.x
        end

        # compute shape Jacobian
        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        detJ = norm2(J)

        # compute B matrix
        for i in 1:hnodes
            for dof in 1:ndim
                NN[dof, (i-1)*ndim + dof              ] = -N[i]
                NN[dof, hnodes*ndim + (i-1)*ndim + dof] =  N[i]
            end
        end

        T = calc_interface_rotation(J)
        @mul B = T'*NN

        if update
            @mul Δω = B*ΔU
            
            Δσ, status = update_state(elem.cmodel, ip.state, ip.cstate, Δω)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        # internal force
        coef = detJ*ip.w*th
        @mul ΔF += coef*B'*Δσ
    end

    return ΔF, map, success()
end


function elem_recover_nodal_values(elem::Element{MechCohesive})
    nips = length(elem.ips)

    keys = output_keys(elem.cmodel)
    vals = zeros(nips, length(keys))
    for (i,ip) in enumerate(elem.ips)
        dict = state_values(elem.cmodel, ip.state)
        vals[i,:] = [ dict[key] for key in keys ]
    end

    node_vals = OrderedDict{Symbol, Vector{Float64}}()
    E = extrapolator(elem.shape, nips)
    for (i,key) in enumerate(keys)
        V = E*vals[:,i]
        node_vals[key] = [ V; V ]
    end

    return node_vals
end


"""
    recover_cohesive_stress(elem::Element{MechCohesive})

Calculates the stress tensor at each integration point of a cohesive element using  regression.

The stress is estimated by performing a least-squares regression on the stress values
from all integration points of the neighboring bulk elements. For each integration point of the
cohesive element, the resulting polynomial is evaluated at its coordinate. Finally, the stress
tensor is projected onto the element's local coordinate system to obtain normal and shear components.

# Arguments
- `elem::Element{MechCohesive}`: The cohesive element for which to calculate the stress.

# Returns
- A vector of vectors, where each inner vector contains the normal and shear stress components
  for one integration point in the element's local coordinate system.

"""
function recover_cohesive_stress(elem::Element{MechCohesive})
    ndim = elem.ctx.ndim
    nnodes = length(elem.nodes)
    nface_nodes = div(nnodes, 2)
    face_coords = get_coords(elem.nodes[1:nface_nodes], ndim)

    # Integration points from bulk elements coupled to the interface.
    neighbor_ips = [ip for bulk_elem in elem.couplings for ip in bulk_elem.ips]
    n_neighbor_ips = length(neighbor_ips)
    n_cohesive_ips = length(elem.ips)
    n_stress_comp = length(neighbor_ips[1].state.σ) # number of stress components from neighboring bulk elements

    @inline function reg_terms(x::Float64, y::Float64, nterms::Int)
        nterms == 4 && return (1.0, x, y, x*y)
        nterms == 3 && return (1.0, x, y)
        return (1.0,)
    end

    @inline function reg_terms(x::Float64, y::Float64, z::Float64, nterms::Int)
        nterms == 4 && return (1.0, x, y, z)
        return (1.0,)
    end

    nterms = if ndim == 3
        n_neighbor_ips >= 4 ? 4 : 1
    else
        n_neighbor_ips >= 4 ? 4 : n_neighbor_ips >= 3 ? 3 : 1
    end

    # Regression matrix at neighbor IP coordinates.
    Mreg = FMat(undef, n_neighbor_ips, nterms)
    for (i, ip) in enumerate(neighbor_ips)
        x, y, z = ip.coord
        Mreg[i, :] .= (ndim == 3) ? reg_terms(x, y, z, nterms) : reg_terms(x, y, nterms)
    end

    # Evaluation matrix at cohesive IP coordinates.
    Meval = FMat(undef, n_cohesive_ips, nterms)
    for (i, ip) in enumerate(elem.ips)
        x, y, z = ip.coord
        Meval[i, :] .= (ndim == 3) ? reg_terms(x, y, z, nterms) : reg_terms(x, y, nterms)
    end

    # Recovered continuum stress components at cohesive IPs.
    σ_eval = FMat(undef, n_stress_comp, n_cohesive_ips)
    σ_sample = FVec(undef, n_neighbor_ips)
    coeffs = FVec(undef, nterms)
    reg_solver = qr(Mreg)

    for icomp in 1:n_stress_comp
        for jip in 1:n_neighbor_ips
            σ_sample[jip] = neighbor_ips[jip].state.σ[icomp]
        end
        coeffs .= reg_solver \ σ_sample
        σ_eval[icomp, :] = Meval * coeffs
    end

    projected_stress = Vector{Float64}[]
    for (i, ip) in enumerate(elem.ips)
        dNdR = elem.shape.deriv(ip.R)
        J = face_coords' * dNdR
        T = calc_interface_rotation(J)

        n1 = T[:, 1]
        n2 = T[:, 2]
        σ = σ_eval[:, i]
        if ndim == 3
            n3 = T[:, 3]
            t1 = dott(σ, n1)
            σn = dot(t1, n1)
            τ1 = dot(t1, n2)
            τ2 = dot(t1, n3)
        else
            n1 = Vec3(n1[1], n1[2], 0.0)
            n2 = Vec3(n2[1], n2[2], 0.0)
            t1 = dott(σ, n1)
            σn = dot(t1, n1)
            τ1 = dot(t1, n2)
            τ2 = 0.0
        end
        push!(projected_stress, cap_stress(elem.cmodel, [σn, τ1, τ2]))
    end

    return projected_stress
end


function stress_strength_ratio(elem::Element{MechCohesive})
    η =  maximum( stress_strength_ratio(elem.cmodel, ip.state.σ) for ip in elem.ips )
    return η
end
