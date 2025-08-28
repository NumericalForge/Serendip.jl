# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechBulk


mutable struct MechBulk<:MechFormulation
    ρ::Float64
    γ::Float64

    function MechBulk(;rho=0.0, gamma=0.0)
        return new(rho, gamma)
    end
end


compat_role(::Type{MechBulk}) = :bulk


function elem_init(elem::Element{MechBulk})
    state_ty = typeof(elem.ips[1].state)
    if :h in fieldnames(state_ty)
        # Volume
        V = cell_extent(elem)*elem.ctx.thickness

        # Representative length size for the element
        nips = length(elem.ips)
        ndim = elem.ctx.ndim
        h = V^(1/ndim)

        for ip in elem.ips
            ip.state.h = h
        end
    end

    return nothing
end


function distributed_bc(elem::Element{MechBulk}, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return mech_boundary_forces(elem, facet, t, key, val)
end


function body_load(elem::Element{MechBulk}, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return mech_solid_body_forces(elem, key, val)
end


function setB(elem::Element, ip::Ip, dNdX::Matx, B::Matx)
    nnodes, ndim = size(dNdX)
    # Note that matrix B is designed to work with tensors in Mandel's notation

    if ndim==2
        for i in 1:nnodes
            j = i-1
            B[1,1+j*ndim] = dNdX[i,1]
            B[2,2+j*ndim] = dNdX[i,2]
            B[6,1+j*ndim] = dNdX[i,2]/SR2
            B[6,2+j*ndim] = dNdX[i,1]/SR2
        end
        if elem.ctx.stress_state==:axisymmetric
            N = elem.shape.func(ip.R)
            for i in 1:nnodes
                j = i-1
                r = ip.coord.x
                B[1,1+j*ndim] = dNdX[i,1]
                B[2,2+j*ndim] = dNdX[i,2]
                B[3,1+j*ndim] =    N[i]/r
                B[6,1+j*ndim] = dNdX[i,2]/SR2
                B[6,2+j*ndim] = dNdX[i,1]/SR2
            end
        end
    else
        for i in 1:nnodes
            dNdx = dNdX[i,1]
            dNdy = dNdX[i,2]
            dNdz = dNdX[i,3]
            j    = i-1
            B[1,1+j*ndim] = dNdx
            B[2,2+j*ndim] = dNdy
            B[3,3+j*ndim] = dNdz
            B[4,2+j*ndim] = dNdz/SR2;   B[4,3+j*ndim] = dNdy/SR2
            B[5,1+j*ndim] = dNdz/SR2;   B[5,3+j*ndim] = dNdx/SR2
            B[6,1+j*ndim] = dNdy/SR2;   B[6,2+j*ndim] = dNdx/SR2
        end
    end

end


function elem_stiffness(elem::Element{MechBulk})
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C = get_coords(elem)
    K = zeros(nnodes*ndim, nnodes*ndim)

    B = zeros(6, nnodes*ndim)
    DB = Array{Float64}(undef, 6, nnodes*ndim)
    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.ctx.stress_state==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        setB(elem, ip, dNdX, B)

        # compute K
        coef = detJ*ip.w*th
        D    = calcD(elem.cmodel, ip.state)
        @mul DB = D*B
        @mul K += coef*B'*DB
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = [ get_dof(node,key).eq_id for node in elem.nodes for key in keys ]

    return K, map, map
end


function elem_mass(elem::Element{MechBulk})
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    ρ = elem.eform.ρ
    C = get_coords(elem)
    M = zeros(nnodes*ndim, nnodes*ndim)
    N = zeros(ndim, nnodes*ndim)
    J = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.ctx.stress_state==:axisymmetric && (th = 2*pi*ip.coord.x)

        # compute N matrix
        Ni   = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)

        for i in 1:nnodes
            for j in 1:ndim
                N[j, (i-1)*ndim+j] = Ni[i]
            end
        end

        @mul J = C'*dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        # compute M
        coef = ρ*detJ*ip.w*th
        @mul M += coef*N'*N
    end

    keys = (:ux, :uy, :uz)[1:ndim]
    map  = Int[ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]

    return M, map, map
end


function elem_internal_forces(elem::Element{MechBulk}, ΔUg::Vector{Float64}=Float64[], dt::Float64=0.0)
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    keys   = (:ux, :uy, :uz)[1:ndim]
    map    = [ get_dof(node, key).eq_id for node in elem.nodes for key in keys ]

    dF = zeros(nnodes*ndim)
    B  = zeros(6, nnodes*ndim)

    J  = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    C = get_coords(elem)

    update = !isempty(ΔUg)
    if update
        ΔU = ΔUg[map]
        Δε = zeros(6)
    end

    for ip in elem.ips
        if elem.ctx.stress_state==:axisymmetric
            th = 2*pi*ip.coord.x
        end

        # compute B matrix
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C'*dNdR
        @mul dNdX = dNdR*inv(J)
        detJ = det(J)
        setB(elem, ip, dNdX, B)

        if update
            @mul Δε = B*ΔU
            Δσ, status = update_state(elem.cmodel, ip.state, Δε)
            failed(status) && return ΔF, map, status
        else
            Δσ = ip.state.σ
        end

        coef = detJ*ip.w*th
        @mul dF += coef*B'*Δσ
    end

    return dF, map, success()
end


function elem_vals(elem::Element{MechBulk})
    vals = OrderedDict{Symbol,Float64}()

    if haskey(state_values(elem.cmodel, elem.ips[1].state), :damt)

        mean_dt = mean( state_values(elem.cmodel, ip.state)[:damt] for ip in elem.ips )

        vals[:damt] = mean_dt
        mean_dc = mean( state_values(elem.cmodel, ip.state)[:damc] for ip in elem.ips )
        vals[:damc] = mean_dc
    end

    #vals = OrderedDict{String, Float64}()
    #keys = elem_vals_keys(elem)
#
    #dicts = [ state_values(elem.cmodel, ip.state) for ip in elem.ips ]
    #nips = length(elem.ips)
#
    #for key in keys
        #s = 0.0
        #for dict in dicts
            #s += dict[key]
        #end
        #vals[key] = s/nips
    #end

    return vals
end


