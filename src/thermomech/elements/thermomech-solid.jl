# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export TMSolid


struct TMSolid <: ThermoMech
    ρ::Float64
    γ::Float64
    cv::Float64
    α::Float64

    function TMSolid(; rho, cv, alpha, gamma=0.0)
        rho >= 0.0 || error("TMSolid: `rho` must be nonnegative.")
        gamma >= 0.0 || error("TMSolid: `gamma` must be nonnegative.")
        cv > 0.0 || error("TMSolid: `cv` must be positive.")
        0.0 <= alpha <= 1.0 || error("TMSolid: `alpha` must satisfy 0 <= alpha <= 1.")
        return new(rho, gamma, cv, alpha)
    end
end


compat_role(::Type{TMSolid}) = :solid


function elem_config_dofs(elem::Element{TMSolid})
    nbnodes = elem.shape.base_shape.npoints
    for (i, node) in enumerate(elem.nodes)
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        elem.ctx.ndim == 3 && add_dof(node, :uz, :fz)
        i <= nbnodes && add_dof(node, :ut, :ft)
    end
end


function distributed_bc(elem::Element{TMSolid}, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return mech_boundary_forces(elem, facet, t, elem.ctx.thickness, key, val)
end


function body_load(elem::Element{TMSolid}, key::Symbol, val::Union{Real,Symbol,Expr})
    return mech_solid_body_forces(elem, elem.ctx.thickness, key, val)
end


@inline function elem_map_u(elem::Element{TMSolid})
    keys = (:ux, :uy, :uz)[1:elem.ctx.ndim]
    return Int[get_dof(node, key).eq_id for node in elem.nodes for key in keys]
end


@inline function elem_map_t(elem::Element{TMSolid})
    nbnodes = elem.shape.base_shape.npoints
    return Int[get_dof(node, :ut).eq_id for node in elem.nodes[1:nbnodes]]
end


@inline function set_Bu(elem::Element{TMSolid}, ip::Ip, dNdX::Matx, B::Matx)
    setB(elem, ip, dNdX, B)
end


function elem_stiffness(elem::Element{TMSolid})
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C = get_coords(elem)
    K = zeros(nnodes * ndim, nnodes * ndim)
    Bu = zeros(6, nnodes * ndim)

    DBu = Array{Float64}(undef, 6, nnodes * ndim)
    J = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        @mul dNdX = dNdR * inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        set_Bu(elem, ip, dNdX, Bu)

        coef = detJ * ip.w * th
        D = calcD(elem.cmodel, ip.state)
        @mul DBu = D * Bu
        @mul K += coef * Bu' * DBu
    end

    map = elem_map_u(elem)
    return K, map, map
end


function elem_coupling_matrix(elem::Element{TMSolid})
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.base_shape.npoints
    E = elem.cmodel.E
    ν = elem.cmodel.ν
    α = elem.etype.α
    C = get_coords(elem)
    Bu = zeros(6, nnodes * ndim)
    Cut = zeros(nnodes * ndim, nbnodes)

    J = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    m = I2
    β = E * α / (1 - 2 * ν)
    if elem.ctx.stress_state == :plane_stress
        β = E * α / (1 - ν)
        m = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    end

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        @mul dNdX = dNdR * inv(J)
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        set_Bu(elem, ip, dNdX, Bu)

        Nt = elem.shape.base_shape.func(ip.R)
        coef = β * detJ * ip.w * th
        mNt = m * Nt'
        @mul Cut -= coef * Bu' * mNt
    end

    return Cut, elem_map_u(elem), elem_map_t(elem)
end


function elem_conductivity_matrix(elem::Element{TMSolid})
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.base_shape.npoints
    C = get_coords(elem)
    H = zeros(nbnodes, nbnodes)
    dNtdX = zeros(nbnodes, ndim)
    Bt = zeros(ndim, nbnodes)
    KBt = zeros(ndim, nbnodes)
    J = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        dNtdR = elem.shape.base_shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNtdX = dNtdR * inv(J)
        Bt .= dNtdX'

        K = calcK(elem.cmodel, ip.state)
        coef = detJ * ip.w * th
        @mul KBt = K * Bt
        @mul H -= coef * Bt' * KBt
    end

    map = elem_map_t(elem)
    return H, map, map
end


function elem_mass_matrix(elem::Element{TMSolid})
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    nbnodes = elem.shape.base_shape.npoints
    C = get_coords(elem)
    M = zeros(nbnodes, nbnodes)

    J = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        Nt = elem.shape.base_shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        coef = elem.etype.ρ * elem.etype.cv * detJ * ip.w * th
        M .-= coef .* (Nt * Nt')
    end

    map = elem_map_t(elem)
    return M, map, map
end


function elem_internal_forces(elem::Element{TMSolid}, F::Vector{Float64})
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    T0k = elem.ctx.T0 + 273.15
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.base_shape.npoints
    C = get_coords(elem)

    E = elem.cmodel.E
    ν = elem.cmodel.ν
    ρ = elem.etype.ρ
    α = elem.etype.α
    cv = elem.etype.cv
    β = E * α / (1 - 2 * ν)

    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

    m = I2
    if elem.ctx.stress_state == :plane_stress
        β = E * α / (1 - ν)
        m = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    end

    dF = zeros(nnodes * ndim)
    Bu = zeros(6, nnodes * ndim)
    dFt = zeros(nbnodes)
    Bt = zeros(ndim, nbnodes)

    J = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    dNtdX = Array{Float64}(undef, nbnodes, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in element $(elem.id)")
        invJ = inv(J)
        @mul dNdX = dNdR * invJ
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.base_shape.deriv(ip.R)
        @mul dNtdX = dNtdR * invJ
        Nt = elem.shape.base_shape.func(ip.R)
        Bt .= dNtdX'

        σ = ip.state.σ
        ut = ip.state.ut

        σ -= β * ut * m
        coef = detJ * ip.w * th
        @mul dF += coef * Bu' * σ

        ε = ip.state.ε
        εvol = dot(m, ε)
        coef = β * εvol * T0k * detJ * ip.w * th
        dFt .-= coef .* Nt

        coef = ρ * cv * detJ * ip.w * th
        dFt .-= coef .* Nt .* ut

        Q = ip.state.QQ
        coef = detJ * ip.w * th
        @mul dFt += coef * Bt' * Q
    end

    F[map_u] .= dF
    F[map_t] .= dFt
end


function update_elem!(elem::Element{TMSolid}, DU::Vector{Float64}, Δt::Float64)
    ndim = elem.ctx.ndim
    th = elem.ctx.thickness
    T0k = elem.ctx.T0 + 273.15
    nnodes = length(elem.nodes)
    nbnodes = elem.shape.base_shape.npoints
    C = get_coords(elem)

    E = elem.cmodel.E
    ν = elem.cmodel.ν
    ρ = elem.etype.ρ
    α = elem.etype.α
    cv = elem.etype.cv
    β = E * α / (1 - 2 * ν)

    map_u = elem_map_u(elem)
    map_t = elem_map_t(elem)

    dU = DU[map_u]
    dUt = DU[map_t]
    Ut = [get_dof(node, :ut).vals[:ut] for node in elem.nodes[1:nbnodes]]
    Ut .+= dUt
    m = I2
    if elem.ctx.stress_state == :plane_stress
        β = E * α / (1 - ν)
        m = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    end

    dF = zeros(nnodes * ndim)
    Bu = zeros(6, nnodes * ndim)
    dFt = zeros(nbnodes)
    Bt = zeros(ndim, nbnodes)

    J = Array{Float64}(undef, ndim, ndim)
    dNdX = Array{Float64}(undef, nnodes, ndim)
    dNtdX = Array{Float64}(undef, nbnodes, ndim)
    Δε = zeros(6)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in element $(elem.id)")
        invJ = inv(J)
        @mul dNdX = dNdR * invJ
        set_Bu(elem, ip, dNdX, Bu)

        dNtdR = elem.shape.base_shape.deriv(ip.R)
        @mul dNtdX = dNtdR * invJ
        Nt = elem.shape.base_shape.func(ip.R)

        @mul Δε = Bu * dU
        Δut = Nt' * dUt
        Bt .= dNtdX'
        G = Bt * Ut

        Δσ, q, status = update_state(elem.cmodel, ip.state, Δε, Δut, G, Δt)
        failed(status) && return [dF; dFt], [map_u; map_t], status

        Δσ -= β * Δut * m
        coef = detJ * ip.w * th
        @mul dF += coef * Bu' * Δσ

        Δεvol = dot(m, Δε)
        coef = β * Δεvol * T0k * detJ * ip.w * th
        dFt .-= coef .* Nt

        coef = ρ * cv * detJ * ip.w * th
        dFt .-= coef .* Nt .* Δut

        coef = Δt * detJ * ip.w * th
        @mul dFt += coef * Bt' * q
    end

    return [dF; dFt], [map_u; map_t], success()
end
