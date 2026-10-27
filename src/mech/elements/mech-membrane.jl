# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

export MechMembrane

"""
    MechMembrane(; thickness, rho=0.0, gamma=0.0)

Small-displacement membrane embedded in 3D, with three translations per node.
Works with 2D area shapes providing interpolation, derivatives and surface
quadrature. `thickness` is positive; `rho` and `gamma` are mass density and
specific weight. There is no bending,
drilling stabilization or geometric stiffness. Constitutive components are
local in-plane Mandel components `(xx, yy, √2 xy)`.
"""
struct MechMembrane <: MechFormulation
    thickness::Float64
    ρ::Float64
    γ::Float64

    function MechMembrane(; thickness::Real=0.0, rho::Real=0.0, gamma::Real=0.0)
        @check isfinite(thickness) && thickness > 0.0 "MechMembrane: thickness must be finite and positive."
        @check isfinite(rho) && rho >= 0.0 "MechMembrane: rho must be finite and nonnegative."
        @check isfinite(gamma) && gamma >= 0.0 "MechMembrane: gamma must be finite and nonnegative."
        return new(thickness, rho, gamma)
    end
end

compat_role(::Type{MechMembrane}) = :surface

struct MechMembraneCache <: ElementCache
    B::Vector{Matrix{Float64}}
    N::Vector{Vector{Float64}}
    weights::Vector{Float64}
end

function elem_config_dofs(elem::Element{MechMembrane})
    elem.ctx.ndim == 3 || error("MechMembrane: requires a 3D analysis.")
    for node in elem.nodes
        add_dof(node, :ux, :fx)
        add_dof(node, :uy, :fy)
        add_dof(node, :uz, :fz)
    end
end

function elem_init(elem::Element{MechMembrane})
    elem.shape.ndim == 2 || error("MechMembrane: expected a 2D area shape, got $(elem.shape.kind).")
    C = get_coords(elem)
    blocks = Matrix{Float64}[]
    shapes = Vector{Float64}[]
    weights = Float64[]
    for ip in elem.ips
        J = C' * elem.shape.deriv(ip.R)
        a, b = J[:,1], J[:,2]
        area = norm(cross(a, b))
        isfinite(area) && area > eps(Float64)*norm(a)*norm(b) || error("MechMembrane: degenerate surface in element $(elem.id).")
        l = normalize(a)
        n = cross(a, b)/area
        m = cross(n, l)
        Rt = [l'; m']
        dNdX = elem.shape.deriv(ip.R) / (Rt*J)
        B = zeros(3, 3length(elem.nodes))
        for i in eachindex(elem.nodes)
            dx, dy = dNdX[i,:]
            B[:,3i-2:3i] .= [dx*l'; dy*m'; (dy*l'+dx*m')/√2]
        end
        push!(blocks, B)
        push!(shapes, elem.shape.func(ip.R))
        push!(weights, area*ip.w)
    end
    area = sum(weights)
    area > 0.0 || error("MechMembrane: nonpositive integrated area.")
    elem.cache = MechMembraneCache(blocks, shapes, weights)
    for ip in elem.ips
        if hasfield(typeof(ip.state), :h)
            ip.state.h = sqrt(area)
        end
        init_state(elem.cmodel, ip.state)
    end
    return nothing
end

function elem_stiffness(elem::Element{MechMembrane})
    cache = elem.cache
    K = zeros(3length(elem.nodes), 3length(elem.nodes))
    for (i, ip) in enumerate(elem.ips)
        B = cache.B[i]
        K += elem.etype.thickness*cache.weights[i] * (B' * calcD(elem.cmodel, ip.state) * B)
    end
    map = dof_map(elem)
    return K, map, map
end

function elem_mass(elem::Element{MechMembrane})
    cache = elem.cache
    M = zeros(3length(elem.nodes), 3length(elem.nodes))
    for i in eachindex(elem.ips)
        mass = elem.etype.ρ*elem.etype.thickness*cache.weights[i] * (cache.N[i]*cache.N[i]')
        for j in eachindex(elem.nodes), k in eachindex(elem.nodes), d in 1:3
            M[3j-3+d,3k-3+d] += mass[j,k]
        end
    end
    map = dof_map(elem)
    return M, map, map
end

function elem_internal_forces(elem::Element{MechMembrane}, ΔU::Vector{Float64}=Float64[], dt::Float64=0.0)
    cache = elem.cache
    F = zeros(3length(elem.nodes))
    map = dof_map(elem)
    for (i, ip) in enumerate(elem.ips)
        if isempty(ΔU)
            stress = ip.state.σ
        else
            stress, status = update_state(elem.cmodel, ip.state, ip.cstate, cache.B[i]*ΔU)
            failed(status) && return F, map, status
        end
        F += elem.etype.thickness*cache.weights[i] * (cache.B[i]'*stress)
    end
    return F, map, success()
end

function distributed_bc(elem::Element{MechMembrane}, facet::Cell, t::Float64, key::Symbol, val::Union{Real,Symbol,Expr,Symbolic})
    return mech_boundary_forces(elem, facet, t, 1.0, key, val)
end

function body_c(elem::Element{MechMembrane}, key::Symbol, val::Union{Real,Symbol,Expr})
    key in (:wx, :wy, :wz) || error("MechMembrane: invalid body-force component $key.")
    component = key == :wx ? :tx : key == :wy ? :ty : :tz
    face = Cell(elem.shape, :surface, elem.nodes, owner=elem)
    return mech_boundary_forces(elem, face, 0.0, elem.etype.thickness, component, val)
end

function elem_vals(elem::Element{MechMembrane})
    vals = OrderedDict{Symbol,Float64}()
    if hasfield(typeof(elem.ips[1].state), :ruptured)
        weights = elem.cache.weights
        # This cell field is an area fraction, not a recovered Boolean flag.
        vals[:ruptured] = sum(weights[i]*ip.state.ruptured for (i, ip) in enumerate(elem.ips))/sum(weights)
    end
    return vals
end

function recovery_samples(elem::Element{MechMembrane})
    return [
        (coord=ip.coord, values=surface_recovery_values(elem.cmodel, ip.state, ip.state.σ))
        for ip in elem.ips
    ]
end

recoverable_fields(elem::Element{MechMembrane}) = [:σ1′, :σ2′, :σvm′, :εp′]
