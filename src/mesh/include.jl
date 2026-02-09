
include("mesh-context.jl")

export Cell, select
include("cell.jl")
include("collapse.jl")

export get_coords, get_node, get_nodes, get_facets, getfaces, getedges, get_patches, cell_extent, cell_quality, cell_aspect_ratio
include("partition.jl")

include("mesh.jl")
include("structured.jl")
# include("unstructured.jl")
include("gen-line-interface-elems.jl")
include("gen-mesh.jl")

include("io.jl")
export Mesh, save, get_outer_facets, threshold

include("refine.jl")
export hrefine, prefine

include("operators.jl")
export move, array, copy, mirror, rotate!, polar, scale, permute_coordinates

include("extrude.jl")
export extrude

include("revolve.jl")
export revolve

include("slice.jl")
export slice

include("smooth.jl")
export smooth!, laplacian_smooth!, fast_smooth!

include("gen-cohesive-elems.jl")
include("gen-contact-elems.jl")
export add_cohesive_elements, add_contact_elements, add_boundary_contact_elements, add_boundary_shell_elements

include("outline.jl")
export get_feature_edges, get_feature_mesh


