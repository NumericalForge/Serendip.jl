# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechCohesive

"""
    MechCohesive()

Zero-thickness cohesive interface formulation for mechanical analyses.
"""
struct MechCohesive<:MechFormulation
end


# Return the shape family that works with this element
compat_role(::Type{MechCohesive}) = :cohesive


mutable struct MechCohesiveCache <: ElementCache
    mobilized::Bool
    depth_list::FixedSizeVector{Float64}
end


function elem_init(elem::Element{MechCohesive})
    # Computation of characteristic length 'h' for cohesive elements
    # and set it in the integration point state

    ndim = elem.ctx.ndim
    state_has_h = hasfield(typeof(elem.ips[1].state), :h)

    # Cohesive depth from linked solids if applicable
    depth_list = if ndim == 2 && elem.ctx.stress_state != :axisymmetric
        length(elem.couplings) == 2 || error("MechCohesive: linked solid thickness requires exactly two linked owners.")
        nips = length(elem.ips)
        list = FixedSizeVector{Float64}(undef, nips)
        for (i, ip) in enumerate(elem.ips)
            th1 = evaluate(elem.couplings[1].etype.thickness, x=ip.coord.x, y=ip.coord.y, z=ip.coord.z)
            th2 = evaluate(elem.couplings[2].etype.thickness, x=ip.coord.x, y=ip.coord.y, z=ip.coord.z)
            list[i] = min(th1, th2)
        end

        list
    else
        list = FixedSizeVector{Float64}(undef, length(elem.ips))
        list .= elem.ctx.thickness
        list
    end

    # Average element size
    if state_has_h
        V = 0.0
        A = 0.0

        if ndim == 2 && elem.ctx.stress_state != :axisymmetric
            for owner in elem.couplings
                C = get_coords(owner)
                J = Array{Float64}(undef, ndim, ndim)
                V_owner = 0.0

                for ip in owner.ips
                    dNdR = owner.shape.deriv(ip.R)
                    @mul J = C'*dNdR
                    detJ = det(J)
                    detJ > 0.0 || error("Negative Jacobian determinant in cell $(owner.id)")
                    depth = evaluate(owner.etype.thickness, x=ip.coord.x, y=ip.coord.y, z=ip.coord.z)
                    V_owner += detJ*ip.w*depth
                end

                V += V_owner
            end
            V /= length(elem.couplings)

            C = get_coords(elem)
            n = div(length(elem.nodes), 2)
            C = C[1:n, :]
            J = fzeros(ndim, ndim-1)

            for (i, ip) in enumerate(elem.ips)
                dNdR = elem.shape.deriv(ip.R)
                @mul J = C'*dNdR
                detJ = norm2(J)
                detJ <= 0 && error("Invalid Jacobian norm for cohesive element")
                A += detJ*ip.w*depth_list[i]
            end
        else
            # Avg volume of linked elements
            for owner in elem.couplings
                V += cell_extent(owner)
            end
            V /= length(elem.couplings)

            # Area of cohesive element
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
        end

        # Calculate and save h at cohesive element's integration points
        h = V/A
        for ip in elem.ips
            ip.state.h = h
        end
    end

    # Explicit cohesive elements already have duplicated/opposite-face nodes.
    # Intrinsic/extrinsic inserted interfaces usually start with the same node ids on both faces.
    n = div(length(elem.nodes), 2)
    mobilized = elem.nodes[1].id != elem.nodes[n+1].id
    elem.cache = MechCohesiveCache(mobilized, depth_list)
end


function elem_stiffness(elem::Element{MechCohesive})
    elem.cache.mobilized == false && return zeros(0,0), Int[], Int[] # return empty arrays if the element is inactive

    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    nstr   = length(elem.ips[1].state.σ)

    C = get_coords(elem)[1:hnodes,:]
    B = fzeros(nstr, nnodes*ndim)
    K = fzeros(nnodes*ndim, nnodes*ndim)

    DB = fzeros(nstr, nnodes*ndim)
    J  = fzeros(ndim, ndim-1)
    NN = fzeros(nstr, nnodes*ndim)

    for (i, ip) in enumerate(elem.ips)
        depth = elem.ctx.stress_state == :axisymmetric ? 2*pi*ip.coord.x : elem.cache.depth_list[i]

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
        if size(T, 1) != nstr
            T = T[1:nstr, 1:nstr]
        end
        @mul B = T'*NN

        # compute K
        coef = detJ*ip.w*depth
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
    nnodes = length(elem.nodes)
    hnodes = div(nnodes, 2) # half the number of total nodes
    map    = dof_map(elem)
    nstr   = length(elem.ips[1].state.σ)

    update = !isempty(ΔU)
    if update
        Δω = zeros(nstr)
    end

    ΔF = fzeros(nnodes*ndim)
    C  = get_coords(elem)[1:hnodes,:]
    B  = fzeros(nstr, nnodes*ndim)

    J  = fzeros(ndim, ndim-1)
    NN = fzeros(nstr, nnodes*ndim)

    for (i, ip) in enumerate(elem.ips)
        depth = elem.ctx.stress_state == :axisymmetric ? 2*pi*ip.coord.x : elem.cache.depth_list[i]

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
        if size(T, 1) != nstr
            T = T[1:nstr, 1:nstr]
        end
        @mul B = T'*NN

        if update
            @mul Δω = B*ΔU
            
            Δσ, status = update_state(elem.cmodel, ip.state, ip.cstate, Δω)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        # internal force
        coef = detJ*ip.w*depth
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
