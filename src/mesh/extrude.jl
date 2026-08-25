# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

"""
    extrude(block; axis=[0,0,1], length=1.0, n=1, quiet=true)

Create a 3D `Block` by extruding a 2D quadrilateral `block` over the specified
`length`. The extrusion converts `QUAD4` blocks to `HEX8` and `QUAD8` blocks to
`HEX20`, with `n` divisions in the extrusion direction.

The `axis` direction is normalized, so its magnitude does not affect the
extrusion length. It can be given as a vector or as `:x`, `:y`, or `:z`.

# Arguments

- `block::Block`: Two-dimensional block to extrude.
- `axis=[0,0,1]`: Extrusion direction.
- `length::Number=1.0`: Total extrusion distance.
- `n::Int=1`: Number of divisions in the extrusion direction.
- `quiet=true`: Accepted for API consistency; this method does not print progress
  information.

# Returns

A new 3D `Block`; the input block is not modified.

# Example

Extrude a quadrilateral block through five divisions along the z-axis:

```julia
block2d = Block([0 0; 1 0; 1 1; 0 1]; nx=3, ny=4, shape=:quad4)
block3d = extrude(block2d; axis=:z, length=2.0, n=5)
```
"""
function extrude(block::Block; axis=[0,0,1], length::Number=1.0, n::Int=1, quiet=true)::Block

    if axis==:x
        axis = Vec3(1,0,0)
    elseif axis==:y
        axis = Vec3(0,1,0)
    elseif axis==:z
        axis = Vec3(0,0,1)
    end

    axis = Vec3(normalize(axis))
    Δl = length
    l = 0.0

    if block.shape==QUAD4
        newshape = HEX8
        ls = [l, l+Δl]
        nidx = [1:4;1:4]
        lidx = [1,1,1,1,2,2,2,2]
    elseif block.shape==QUAD8
        newshape = HEX20
        ls = [l, l+Δl/2, l+Δl]
        nidx = [1:4;1:4;5:8;5:8;1:4]
        lidx = [1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,2,2,2,2]
    else
        error("extrude: Block shape $(block.shape.kind) is not supported")
    end

    points = Point[]
    for (point,li) in zip(block.points[nidx], ls[lidx])
        coord = point.coord + li*axis
        push!(points, Point(coord))
    end

    return Block(points, nx=block.nx, ny=block.ny, nz=n, shape=newshape)

end


function extrude(blocks::AbstractArray; axis=[0,0,1], length=1.0::Number, n=1::Int, quiet=true)
    return extrude.(blocks, axis=axis, length=length, n=n, quiet=quiet)
end


# Generates a new mesh obtained by extrusion of a 2D mesh
"""
    extrude(mesh; length=1.0, n=1, axis=nothing, quiet=true, lagrangian=false)

Create a mesh one dimension higher by extruding each element of `mesh` through
`n` layers over the specified `length`. A 1D mesh produces a 2D mesh, while a
2D mesh produces a 3D mesh.

When `axis` is a vector, it is normalized and used as the extrusion direction.
When `axis=nothing`, the direction is computed from the local mesh normals.

# Arguments

- `length::Real=1.0`: Total extrusion distance.
- `n::Int=1`: Number of element layers in the extrusion direction. Must be positive.
- `axis=nothing`: Extrusion direction, given as a vector, or `nothing` to use the
  local mesh normals.
- `quiet=true`: Suppress progress information when `true`.
- `lagrangian=false`: Generate Lagrangian elements when available.

# Returns

A new `Mesh`; the input mesh is not modified.

# Example

Extrude one quadrilateral into two hexahedral layers:

```julia
coordinates = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
mesh2d = Mesh(coordinates, [[1, 2, 3, 4]], [:quad4])
mesh3d = extrude(mesh2d; axis=[0, 0, 1], length=2.0, n=2)
```
"""
function extrude(mesh::Mesh; length::Real=1.0, n::Int=1, axis=nothing, quiet=true, lagrangian=false)

    quiet || printstyled("Mesh extrude:\n", bold=true, color=:cyan)

    @check n>0

    if axis !== nothing
        if axis==:x
            axis = Vec3(1,0,0)
        elseif axis==:y
            axis = Vec3(0,1,0)
        elseif axis==:z
            axis = Vec3(0,0,1)
        else
            ax = Vec3(normalize(axis))
        end
    end

    # check cells
    for cell in mesh.elems
        celldim = cell.shape.ndim
        celldim==1 && mesh.ctx.ndim==3 && axis===nothing && error("extrude: cannot extrude cell of shape $(cell.shape.kind) in dimension 3 using normal")
        celldim==3 && error("extrude: cannot extrude cell of shape $(cell.shape.kind)")
    end

    # compute normals
    if axis===nothing
        nnodes  = Base.length(mesh.nodes)
        normals = zeros(nnodes, 3)
        counts  = zeros(nnodes)
        for cell in mesh.elems
            celldim = cell.shape.ndim
            # coords  = get_coords(cell, celldim)
            coords  = get_coords(cell)
            natcoords = cell.shape.nat_coords
            for (i,node) in enumerate(cell.nodes)
                R = natcoords[i,:]
                J = coords'*cell.shape.deriv(R)
                if celldim==1
                    N = Vec3(normalize([-J[2], J[1]]))
                else
                    N = Vec3(normalize(cross(J[:,1], J[:,2])))
                end
                normals[node.id, :] .+= N
                counts[node.id] += 1
            end
        end

        normals ./= counts
    end

    # generate new cells
    cells = Cell[]
    Δl    = length/n
    for l in range(0, step=Δl, length=n)
        for cell in mesh.elems

            if cell.shape==LIN2
                newshape = QUAD4
                ls = [l, l+Δl]
                nidx = [1,2,2,1]
                lidx = [1,1,2,2]
            elseif cell.shape==LIN3
                newshape = QUAD8
                ls = [l, l+Δl/2, l+Δl]
                nidx = [1,2,2,1,3,2,3,1]
                lidx = [1,1,3,3,1,2,3,2]
                if lagrangian
                    push!(nidx, 3)
                    push!(lidx, 2)
                end
            elseif cell.shape==LIN4
                newshape = QUAD8
                ls = [l, l+Δl*1/3, l+Δl*2/3, l+Δl]
                nidx = [1,2,2,1,3,4,2,2,4,3,1,1]
                lidx = [1,1,4,4,1,1,2,3,4,4,3,2]
            elseif cell.shape==TRI3
                newshape = WED6
                ls = [l, l+Δl]
                nidx = [1,2,3,1,2,3]
                lidx = [1,1,1,2,2,2]
            elseif cell.shape==QUAD4
                newshape = HEX8
                ls = [l, l+Δl]
                nidx = [1:4;1:4]
                lidx = [1,1,1,1,2,2,2,2]
            elseif cell.shape==TRI6
                newshape = WED15
                ls = [l, l+Δl/2, l+Δl]
                nidx = [1:3;1:3;4:6;4:6;1:3]
                lidx = [1,1,1,3,3,3,1,1,1,3,3,3,2,2,2]
            elseif cell.shape==QUAD8
                newshape = HEX20
                ls = [l, l+Δl/2, l+Δl]
                nidx = [1:4;1:4;5:8;5:8;1:4]
                lidx = [1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,2,2,2,2]
            elseif cell.shape==QUAD9
                newshape = HEX27
                ls = [l, l+Δl/2, l+Δl]
                nidx = [1:4;1:4;5:8;5:8;1:4; 8;6;5;7;9;9;9]
                lidx = [1,1,1,1,3,3,3,3,1,1,1,1,3,3,3,3,2,2,2,2, 2,2,2,2,1,3,2]
            else
                error("extrude: Cell shape $(cell.shape.kind) is not supported")
            end

            nodes = Node[]
            for (node,li) in zip(cell.nodes[nidx], ls[lidx])
                if axis===nothing
                    ax = Vec3(normals[node.id, :])
                end
                coord = node.coord + li*ax
                push!(nodes, Node(coord))
            end

            role = newshape.ndim==2 ? :surface : :solid
            newcell = Cell(newshape, role, nodes, tag=cell.tag)
            isinverted(newcell) && flip(newcell)
            push!(cells, newcell)
        end
    end

    # Merge coincident points by spatial position, not object identity.
    point_d = NodePosMap(n => n for c in cells for n in c.nodes)

    for cell in cells
        cell.nodes = Node[point_d[n] for n in cell.nodes]
    end
    nodes = collect(values(point_d))

    # new mesh
    newmesh = Mesh()
    newmesh.nodes = nodes
    newmesh.elems = cells
    synchronize(newmesh, sort=true)

    if !quiet
        @printf "  %5d points obtained\n" Base.length(newmesh.nodes)
        @printf "  %5d cells obtained\n" Base.length(newmesh.elems)
        @printf "  %5d faces obtained\n" Base.length(newmesh.faces)
        @printf "  %5d surface edges obtained\n" Base.length(newmesh.edges)
    end

    return newmesh

end
