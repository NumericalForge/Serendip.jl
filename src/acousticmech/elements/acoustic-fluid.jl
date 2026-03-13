# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export AcousticFluid


struct AcousticFluid<:AcousticMech
    ρ::Float64
    c::Float64

    function AcousticFluid(; rho=0.0, c=0.0)
        ρ = Float64(rho)
        c = Float64(c)
        @check ρ >= 0.0
        @check c > 0.0
        return new(ρ, c)
    end
end

compat_role(::Type{AcousticFluid}) = :solid


function distributed_bc(elem::Element{AcousticFluid}, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return acoustic_mech_bc(elem, facet, t, key, val)
end


function elem_acoustic_stiffness(elem::Element{AcousticFluid})
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C      = get_coords(elem)
    K      = zeros(nnodes, nnodes)
    dNdX   = zeros(nnodes, ndim)
    J      = Array{Float64}(undef, ndim, ndim)

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")
        @mul dNdX = dNdR * inv(J)
        Bp = dNdX'

        coef = detJ * ip.w * th
        K .+= coef * Bp' * Bp
    end

    map = dof_map(elem)
    return K, map, map
end


function elem_acoustic_mass(elem::Element{AcousticFluid})
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nnodes = length(elem.nodes)
    C      = get_coords(elem)
    M      = zeros(nnodes, nnodes)
    J      = Array{Float64}(undef, ndim, ndim)
    c      = elem.etype.c

    for ip in elem.ips
        elem.ctx.stress_state == :axisymmetric && (th = 2 * pi * ip.coord.x)

        N    = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        @mul J = C' * dNdR
        detJ = det(J)
        detJ > 0.0 || error("Negative Jacobian determinant in cell $(elem.id)")

        coef = detJ * ip.w * th / c^2
        M .+= coef * N * N'
    end

    map = dof_map(elem)
    return M, map, map
end


function elem_internal_forces(elem::Element{AcousticFluid}, ΔU::Vector{Float64}=Float64[], Δt::Float64=0.0)
    map = dof_map(elem)
    isempty(ΔU) && return zeros(length(map)), map, success()

    dF = elem_acoustic_stiffness(elem)[1] * ΔU
    return dF, map, success()
end
