
function get_facet_normal(face::AbstractCell)
    ndim = 1 + face.shape.ndim
    C = get_coords(face, ndim)

    if ndim==2
        C .+= [pi pi^1.1]
    else
        C .+= [pi pi^1.1 pi^1.2]
    end

    # calculate the normal
    I = ones(size(C,1))
    N = pinv(C)*I # best fit normal
    normalize!(N) # get unitary vector

    return N
end


function get_outline_edges(cells::Vector{<:AbstractCell}; angle=120, tol::Union{Nothing,Real}=nothing)

    boundary_faces = Dict{Tuple, Cell}()

    # Get faces
    for cell in cells
        cell.role in (:solid, :surface) || continue # only bulk cells
        if cell.shape.ndim==2
            key = _cell_key(cell, tol=tol)
            if haskey(boundary_faces, key)
                delete!(boundary_faces, key)
            else
                boundary_faces[key] = cell
            end
            continue
        end

        for face in get_facets(cell)
            key = _cell_key(face, tol=tol)
            if haskey(boundary_faces, key)
                delete!(boundary_faces, key)
            else
                boundary_faces[key] = face
            end
        end
    end

    faces = collect(values(boundary_faces))

    # Get normals
    normals = IdDict{Cell,Vector{Float64}}( f => get_facet_normal(f) for f in faces )
    pending_edges = Dict{Tuple,Cell}()
    outline_edges = CellEdge[]

    # Get edges with non-coplanar adjacent faces
    for face in faces
        n1 = normals[face]
        for edge in get_edges(face)
            key = _cell_key(edge, tol=tol)
            edge0 = get(pending_edges, key, nothing)
            if edge0===nothing
                pending_edges[key] = edge
            else
                delete!(pending_edges, key)
                n2 = normals[edge0.owner]
                α = 180 - acosd( abs(clamp(dot(n1,n2),-1,1)) )
                α = round(α, digits=2)
                α<=angle && push!(outline_edges, edge)
            end
        end
    end
    append!(outline_edges, values(pending_edges))

    return outline_edges
end
