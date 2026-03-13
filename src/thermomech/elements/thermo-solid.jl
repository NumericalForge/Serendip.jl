# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export ThermoSolid


struct ThermoSolid <: ThermoMech
    ρ::Float64

    function ThermoSolid(; rho=0.0)
        rho >= 0.0 || error("ThermoSolid: `rho` must be nonnegative.")
        return new(rho)
    end
end


compat_role(::Type{ThermoSolid}) = :solid


function elem_config_dofs(elem::Element{ThermoSolid})
    for node in elem.nodes
        add_dof(node, :ut, :ft)
    end
end


@inline function elem_map_t(elem::Element{ThermoSolid})
    return Int[get_dof(node, :ut).eq_id for node in elem.nodes]
end


function distributed_bc(elem::Element{ThermoSolid}, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    key == :tq || error("distributed_bc: boundary condition $key is not applicable in a ThermoSolid element")

    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    nodes = facet.nodes
    nnodes = length(nodes)
    C = get_coords(nodes, ndim)

    F = zeros(nnodes)
    J = Array{Float64}(undef, ndim, ndim)
    ips = get_ip_coords(facet.shape)

    for i in 1:size(ips, 1)
        R = ips[i].coord
        w = ips[i].w
        N = facet.shape.func(R)
        D = facet.shape.deriv(R)
        @mul J = C' * D
        X = C' * N

        if ndim == 2
            x, y = X
            vip = evaluate(val, t=t, x=x, y=y)
            elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * X[1])
        else
            x, y, z = X
            vip = evaluate(val, t=t, x=x, y=y, z=z)
        end

        coef = vip * norm(J) * w * th
        F .+= coef .* N
    end

    map = Int[get_dof(node, :ut).eq_id for node in facet.nodes]
    return F, map
end


function elem_conductivity_matrix(elem::Element{ThermoSolid})
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C = get_coords(elem)
    H = zeros(nnodes, nnodes)
    Bt = zeros(ndim, nnodes)
    KBt = zeros(ndim, nnodes)

    J = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNdX = dNdR * inv(J)
        Bt .= dNdX'

        K = calcK(elem.cmodel, ip.state)
        coef = detJ * ip.w * th
        @mul KBt = K * Bt
        @mul H -= coef * Bt' * KBt
    end

    map = elem_map_t(elem)
    return H, map, map
end


function elem_mass_matrix(elem::Element{ThermoSolid})
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    ρ = elem.etype.ρ
    C = get_coords(elem)
    nnodes = length(elem.nodes)
    M = zeros(nnodes, nnodes)

    J = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        cv = calc_cv(elem.cmodel, ip.state.ut)
        coef = ρ * cv * detJ * ip.w * th
        M .-= coef .* (N * N')
    end

    map = elem_map_t(elem)
    return M, map, map
end


function elem_internal_forces(elem::Element{ThermoSolid}, F::Vector{Float64})
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    nnodes = length(elem.nodes)
    map_t = elem_map_t(elem)

    C = get_coords(elem)
    dFt = zeros(nnodes)
    Bt = zeros(ndim, nnodes)

    J = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNdX = dNdR * inv(J)

        Bt .= dNdX'

        cv = calc_cv(elem.cmodel, ip.state.ut)
        coef = detJ * ip.w * elem.etype.ρ * cv * th
        dFt .-= coef .* N .* ip.state.ut

        Q = ip.state.Q
        coef = detJ * ip.w * th
        @mul dFt += coef * Bt' * Q
    end

    F[map_t] .= dFt
end


function update_elem!(elem::Element{ThermoSolid}, DU::Vector{Float64}, Δt::Float64)
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    nnodes = length(elem.nodes)
    map_t = elem_map_t(elem)
    C = get_coords(elem)

    dUt = DU[map_t]
    Ut = [get_dof(node, :ut).vals[:ut] for node in elem.nodes]
    Ut .+= dUt

    Bt = zeros(ndim, nnodes)
    dFt = zeros(nnodes)

    J = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNdX = dNdR * inv(J)

        Bt .= dNdX'
        G = Bt * Ut
        Δut = N' * dUt
        q = update_state(elem.cmodel, ip.state, Δut, G, Δt)

        cv = calc_cv(elem.cmodel, ip.state.ut)
        coef = elem.etype.ρ * cv * detJ * ip.w * th
        dFt .-= coef .* N .* Δut

        coef = Δt * detJ * ip.w * th
        @mul dFt += coef * Bt' * q
    end

    return dFt, map_t, success()
end
