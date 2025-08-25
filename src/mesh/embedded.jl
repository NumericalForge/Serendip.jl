# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

# Remove line-joint-cells from mesh and set up line cells as embedded cells
function generate_embedded_cells!(mesh::Mesh)

    newcells = []
    id = 0
    for cell in mesh.elems
        if cell.role==BONDSLIP_CELL
            # link solid cell to line cells
            solid, line = cell.couplings
            line.couplings = [ solid ]
        else
            # save non joint1D cells
            id += 1
            cell.id = id  # update id
            push!(newcells, cell)
        end
    end

    # update mesh
    mesh.elems = newcells
    return mesh
end
