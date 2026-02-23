# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechBondSlip

"""
    MechBondSlip(; p)

Create a mechanical line-interface formulation for bond–slip between a rebar and surrounding material.

# Arguments
- `p::Real`: rebar perimeter (length).

# Behavior
- Pair with a bond–slip constitutive model.

# Returns
- A `MechBondSlip` object.
"""
struct MechBondSlip<:MechFormulation
    p::Float64

    function MechBondSlip(; p=NaN)
        @check p > 0.0
        return new(p)
    end
end

compat_role(::Type{MechBondSlip}) = :line_interface


function set_quadrature(elem::Element{MechBondSlip}, n::Int=0; state::NamedTuple=NamedTuple())

    bar   = elem.couplings[2]
    if n==0
        if bar.shape == LIN2
            n = 2
        else
            n = 3
        end
    end

    @check n in (2,3) "Only 2 or 3 integration points are allowed in MechBondSlip elements"

    ipc = get_ip_coords(bar.shape, n)
    n   = size(ipc,1)

    resize!(elem.ips, n)
    for i in 1:n
        R = ipc[i].coord
        w = ipc[i].w
        ipstate = compat_state_type(typeof(elem.cmodel), typeof(elem.etype))(elem.ctx; state...)
        elem.ips[i] = Ip(R, w, elem, ipstate)
    end

    # finding ips global coordinates
    C = get_coords(bar)
    for ip in elem.ips
        N = bar.shape.func(ip.R)
        ip.coord = extend!(C'*N, 3)
    end
end


function elem_init(elem::Element{MechBondSlip})
    host = elem.couplings[1]
    bar  = elem.couplings[2]
    Ch = get_coords(host)
    Ct = get_coords(bar)
    nips = length(elem.ips)
    elem.cacheM = Vector{FixedSizeMatrix{Float64}}(undef,nips) # B matrices
    elem.cacheV = [ FixedSizeVector{Float64}(undef,nips) ] # detJ values at slot 1
    for (i,ip) in enumerate(elem.ips)
        B, detJ = mountB(elem, ip.R, Ch, Ct)
        elem.cacheM[i] = FixedSizeMatrix(B)
        elem.cacheV[1][i] = detJ
    end

    return nothing
end


function mount_T(J::Matx)
    ndim = length(J)
    nJ   = norm(J)
    L1   = vec(J/nJ)

    if ndim==2
        L2 = [ -L1[2],  L1[1] ]
        return hcat(L1,L2)'
    end

    # Finding second vector
    if     abs(L1[1]) == 1.0; L2 = [0.0, 1.0, 0.0]
    elseif abs(L1[2]) == 1.0; L2 = [0.0, 0.0, 1.0]
    elseif abs(L1[3]) == 1.0; L2 = [1.0, 0.0, 0.0]
    else
        # Auxiliar vector L which must be different from L1
        L = [1.0, 0.0, 0.0]
        if norm(L-L1) < 1.0e-4; L = [0.0, 1.0, 0.0] end
        # Performing cross product to obtain a second vector
        L2  = normalize(cross(L1, L))
    end

    # Finding third vector
    L3 = normalize(cross(L1, L2))

    return [ L1'; L2'; L3' ]
end


function mountB(elem::Element{MechBondSlip}, R, Ch, Ct)
    # Calculates the matrix that relates nodal displacements with relative displacements


    # B = T* [NN*MM  -NN]      ndim x ndim*(m+n)

    # where
    # T is a direction cosines matrix
    # NN is a matrix containing truss element shape functions
    # evaluated at the point of interest R.
    # MM is a matrix containing tresspased element shape functions
    # evaluated at the n truss nodal points.

    #          [ M_11*I M_21*I ... M_m1*I]
    #     MM = [ M_12*I M_22*I ... M_m2*I]
    #          [ M_13*I M_23*I ... M_mn*I] ndim*n x ndim*m

    # where
    # M_12 is the first shape function from the tresspased solid element
    # evaluated at the second node of the truss element.
    # I is a ndim x ndim identity matrix


    ndim = elem.ctx.ndim
    host = elem.couplings[1]
    bar  = elem.couplings[2]
    # nnodes  = length(elem.nodes)
    nbnodes = length(bar.nodes)
    nsnodes = length(host.nodes)
    D = bar.shape.deriv(R)
    J = Ct'*D
    T = mount_T(J)

    # Mount NN matrix
    N = bar.shape.func(R)
    NN = hcat([ Ni*eye(ndim) for Ni in N  ]...)

    # Mount MM matrix
    stack = Array{Float64,2}[]
    for i in 1:nbnodes
        Xj = bar.nodes[i].coord
        R  = inverse_map(host.shape, Ch, Xj)
        M  = host.shape.func(R)
        for Mi in M
            push!(stack, Mi*eye(ndim))
        end
    end
    MM = hvcat(nsnodes, stack...)

    B = T*[ -NN*MM  NN ]
    detJ = norm(J)
    return B, detJ
end


function elem_stiffness(elem::Element{MechBondSlip})
    ndim = elem.ctx.ndim
    nnodes = length(elem.nodes)
    p = elem.etype.p

    K  = zeros(nnodes*ndim, nnodes*ndim)
    DB = zeros(ndim, nnodes*ndim)

    for (i,ip) in enumerate(elem.ips)
        B    = elem.cacheM[i]
        detJ = elem.cacheV[1][i]
        D    = calcD(elem.cmodel, ip.state)
        coef = p*detJ*ip.w
        @mul DB = D*B
        @mul K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ get_dof(node,key).eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_internal_forces(elem::Element{MechBondSlip}, ΔU::Vector{Float64}=Float64[], dt::Float64=0.0)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    p = elem.etype.p

    map  = dof_map(elem)

    update = !isempty(ΔU)
    if update
        Δu = zeros(ndim)
    end

    ΔF = zeros(nnodes*ndim)

    for (i,ip) in enumerate(elem.ips)
        B    = elem.cacheM[i]
        detJ = elem.cacheV[1][i]

        if update
            @mul Δu = B*ΔU
            Δσ, status = update_state(elem.cmodel, ip.state, ip.cstate, Δu)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end
        coef = p*detJ*ip.w
        @mul ΔF += coef*B'*Δσ
    end

    return ΔF, map, success()
end


function elem_recover_nodal_values(elem::Element{MechBondSlip})
    all_ip_vals = [ state_values(elem.cmodel, ip.state) for ip in elem.ips ]
    nips        = length(elem.ips)
    fields      = keys(all_ip_vals[1])
    nfields     = length(fields)

    # matrix with all ip values (nip x nvals)
    W = mapreduce(transpose, vcat, collect.(values.(all_ip_vals)))

    host = elem.couplings[1]
    bar  = elem.couplings[2]

    E = extrapolator(bar.shape, nips)
    N = E*W # (nbnodes x nfields)

    nhnodes = length(host.nodes)
    N = [ zeros(nhnodes, nfields); N ]

    # Filling nodal and elem vals
    node_vals = OrderedDict{Symbol, Vector{Float64}}(field => N[:,i] for (i,field) in enumerate(fields))

    return node_vals
end
