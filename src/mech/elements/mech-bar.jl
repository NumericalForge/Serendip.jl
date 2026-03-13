# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechBar

"""
    MechBar(; d=0.0, A=0.0)

Mechanical formulation for axial bar elements.
Defines uniaxial axial behavior along the element length for truss/rod analyses.

# Keyword arguments
- `d::Float64=0.0`: Circular section diameter (`d > 0`). When provided, it is used to compute the section area `A`.
- `A::Float64=0.0`: Cross-sectional area (`A > 0`) for direct area-based definition.

# Stored fields
- `A::Float64`: Cross-sectional area used by stiffness/mass routines.

"""
struct MechBar<:MechFormulation
    A::Float64

    function MechBar(; d::Float64=0.0, A::Float64=0.0)
        if d > 0.0
            @check d > 0.0 "MechBar: diameter d must be > 0 for circular sections"
            @check A == 0.0 "MechBar: if diameter d is specified, A must be 0"
            A = π*(d/2)^2
        else
            @check A > 0.0 "MechBar: cross-sectional area A must be > 0"
        end

        return new(A)
    end
end


compat_role(::Type{MechBar}) = :line
embedded_formulation(::Type{MechBar}) = MechEmbBar


function distributed_bc(elem::Element{MechBar}, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return mech_line_distributed_forces(elem, t, key, val)
end


function body_load(elem::Element{MechBar}, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_line_distributed_forces(elem, 0.0, key, val)
end


function elem_stiffness(elem::Element{MechBar})
    local E::Float64, A::Float64, coef::Float64, dNdR::Matrix{Float64}

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
    keys = [:ux, :uy, :uz][1:ndim]
    map  = [ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]
    return K, map, map
end


function elem_mass(elem::Element{MechBar})

    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    ρ = elem.etype.ρ
    A = elem.etype.A

    C = get_coords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, 1)
    N = zeros(ndim, ndim*nnodes)

    for ip in elem.ips

        dNdR = elem.shape.deriv(ip.R)
        Ni = elem.shape.func(ip.R)
        setNt(ndim,Ni,N)

        @mul J = C'*dNdR
        detJ = norm(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute M
        coef = ρ*A*detJ*ip.w
        @mul M += coef*N'*N

    end

    keys = [:ux, :uy, :uz][1:ndim]
    map  = [ get_dof(node,key).eq_id for node in elem.nodes for key in keys ]
    return M, map, map
end


function setNt(ndim::Int,Ni::Vect, N::Matx)
    nnodes = length(Ni)
    N .= 0.0

    if ndim==2
        for i in 1:nnodes
            j = i-1
            N[1,1+j*ndim] = Ni[i]
            N[2,2+j*ndim] = Ni[i]
        end
    elseif ndim==3
        for i in 1:nnodes
            j    = i-1
            N[1,1+j*ndim] = Ni[i]
            N[2,2+j*ndim] = Ni[i]
            N[3,3+j*ndim] = Ni[i]
       end
    else
        for i in 1:nodes
            j = i-1
            N[1,1+j*ndim] = Ni[i]
        end
    end

end


function elem_internal_forces(elem::Element{MechBar}, ΔU::Vector{Float64}=Float64[], dt::Float64=0.0)
    ndim   = elem.ctx.ndim
    nnodes = length(elem.nodes)
    A      = elem.etype.A
    map    = dof_map(elem)

    ΔF = zeros(nnodes*ndim)
    C = get_coords(elem)
    B = zeros(1, nnodes*ndim)
    J = Array{Float64}(undef, ndim, 1)

    update = !isempty(ΔU)

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
            Δε = (B*ΔU)[1]
            Δσ, status = update_state(elem.cmodel, ip.state, ip.cstate, Δε)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        coef = A*detJ*ip.w
        ΔF .+= coef*Δσ*vec(B')
    end

    return ΔF, map, success()
end


function elem_vals(elem::Element{MechBar})
    # get ip average values
    ipvals = [ state_values(elem.cmodel, ip.state) for ip in elem.ips ]
    merger(x,y) = abs(x) > abs(y) ? x : y
    merged  = merge(merger, ipvals... )
    vals = OrderedDict( k=>v for (k,v) in merged)

    # sum  = merge(+, ipvals... )
    # nips = length(elem.ips)
    # vals = OrderedDict( k=>v/nips for (k,v) in sum)
    return vals
end
