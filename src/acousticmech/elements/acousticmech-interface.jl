# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export AcousticMechInterface


struct AcousticMechInterface<:ElementFormulation
    ρ::Float64

    function AcousticMechInterface(; rho::Float64=0.0)
        @check rho >= 0.0
        return new(rho)
    end
end


compat_role(::Type{AcousticMechInterface}) = :contact


function elem_config_dofs(elem::Element{AcousticMechInterface})
    nnodes = length(elem.nodes)
    @check iseven(nnodes) "elem_config_dofs: AcousticMechInterface requires duplicated face nodes."
    return nothing
end


function interface_structural_map(elem::Element{AcousticMechInterface})
    ndim  = elem.ctx.ndim
    keys  = (:ux, :uy, :uz)[1:ndim]
    nface = div(length(elem.nodes), 2)

    structural_nodes = get_dof(elem.nodes[1], :up) === nothing ? elem.nodes[1:nface] : elem.nodes[nface+1:end]

    @check all(node -> get_dof(node, :ux) !== nothing, structural_nodes) "interface_structural_map: structural half nodes must contain displacement dofs."
    return Int[get_dof(node, key).eq_id for node in structural_nodes for key in keys]
end


function interface_pressure_map(elem::Element{AcousticMechInterface})
    nface = div(length(elem.nodes), 2)

    pressure_nodes = get_dof(elem.nodes[1], :up) !== nothing ? elem.nodes[1:nface] : elem.nodes[nface+1:end]

    @check all(node -> get_dof(node, :up) !== nothing, pressure_nodes) "interface_pressure_map: pressure half nodes must contain `:up` dofs."
    return Int[get_dof(node, :up).eq_id for node in pressure_nodes]
end


function dof_map(elem::Element{AcousticMechInterface})
    return [interface_structural_map(elem); interface_pressure_map(elem)]
end


function interface_fluid_element(elem::Element{AcousticMechInterface})
    fluid_elems = [coupling for coupling in elem.couplings if hasmethod(elem_acoustic_mass, (typeof(coupling),))]
    @check length(fluid_elems) == 1 "interface_fluid_element: expected exactly one coupled acoustic element."
    return fluid_elems[1]
end


function interface_reference_density(elem::Element{AcousticMechInterface})
    elem.etype.ρ > 0.0 && return elem.etype.ρ

    fluid_elem = interface_fluid_element(elem)
    hasproperty(fluid_elem.etype, :ρ) || error("interface_reference_density: coupled acoustic element must provide density.")
    return fluid_elem.etype.ρ
end


function interface_normal_operator(
    elem::Element{AcousticMechInterface},
    N::AbstractVector,
    J::AbstractMatrix,
    XΓ::AbstractVector,
    Xf::AbstractVector,
)
    ndim = elem.ctx.ndim
    nface = length(N)
    L = zeros(1, nface * ndim)

    n = if ndim == 2
        normalize([J[2], -J[1]])
    else
        normalize(cross(J[:, 1], J[:, 2]))
    end

    # Use the fluid centroid to orient the interface normal outward from the fluid domain.
    dot(n, XΓ - Xf) < 0.0 && (n .*= -1.0)

    for i in 1:nface
        for d in 1:ndim
            L[1, (i - 1) * ndim + d] = N[i] * n[d]
        end
    end

    return L
end


function elem_coupling_matrix(elem::Element{AcousticMechInterface})
    ndim   = elem.ctx.ndim
    th     = elem.ctx.thickness
    nface = div(length(elem.nodes), 2)
    first_half = elem.nodes[1:nface]
    last_half = elem.nodes[nface+1:end]
    pressure_first = get_dof(first_half[1], :up) !== nothing
    nodes_p = pressure_first ? first_half : last_half

    @check all(node -> get_dof(node, :up) !== nothing, nodes_p) "elem_coupling_matrix: pressure half nodes must contain `:up` dofs."
    nnodes = length(nodes_p)
    Cgeom  = get_coords(nodes_p, ndim)
    Ce     = zeros(nnodes, nnodes * ndim)
    fluid_elem = interface_fluid_element(elem)
    Xf = vec(sum(get_coords(fluid_elem), dims=1) ./ length(fluid_elem.nodes))

    for ip in elem.ips
        N = elem.shape.func(ip.R)
        dNdR = elem.shape.deriv(ip.R)
        J = Cgeom' * dNdR
        detJ = norm2(J)
        detJ > 0.0 || error("elem_coupling_matrix: invalid interface Jacobian for element $(elem.id)")

        coef = detJ * ip.w * th
        if elem.ctx.stress_state == :axisymmetric
            X = Cgeom' * N
            coef *= 2 * pi * X[1]
        end

        XΓ = vec(Cgeom' * N)
        L = interface_normal_operator(elem, N, J, XΓ, Xf)
        Ce .+= coef * (N * L)
    end

    return Ce, interface_pressure_map(elem), interface_structural_map(elem)
end


function elem_interface_inertia_matrix(elem::Element{AcousticMechInterface})
    Ce, pmap, umap = elem_coupling_matrix(elem)
    return interface_reference_density(elem) * Ce, pmap, umap
end


function elem_stiffness(elem::Element{AcousticMechInterface})
    return zeros(0, 0), Int[], Int[]
end


function elem_mass(elem::Element{AcousticMechInterface})
    return zeros(0, 0), Int[], Int[]
end


function elem_internal_forces(elem::Element{AcousticMechInterface}, ΔU::Vector{Float64}=Float64[], Δt::Float64=0.0)
    map = dof_map(elem)
    isempty(ΔU) && return zeros(length(map)), map, success()

    umap = interface_structural_map(elem)
    pmap = interface_pressure_map(elem)
    nu = length(umap)
    np = length(pmap)
    @check length(ΔU) == nu + np "elem_internal_forces: incompatible increment size for AcousticMechInterface."

    ΔUp = ΔU[nu+1:end]
    Ce, _, _ = elem_coupling_matrix(elem)

    ΔF = zeros(nu + np)
    ΔF[1:nu] .= Ce' * ΔUp
    return ΔF, map, success()
end
