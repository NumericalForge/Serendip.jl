# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

using TranscodingStreams, CodecZlib

function save_vtk(mesh::AbstractDomain, filename::String; desc::String="")
    # Saves a UnstructuredGrid
    npoints = length(mesh.nodes)
    ncells  = length(mesh.elems)

    # Number of total connectivities
    nconns = 0
    for cell in mesh.elems
        nconns += 1 + length(cell.nodes)
    end

    # Open filename
    f = open(filename, "w")

    println(f, "# vtk DataFile Version 3.0")
    println(f, desc)
    println(f, "ASCII")
    println(f, "DATASET UNSTRUCTURED_GRID")
    println(f, "")
    println(f, "POINTS ", npoints, " float64")

    # Write nodes
    for (i,node) in enumerate(mesh.nodes)
        @printf f "%23.15e %23.15e %23.15e \n" node.coord.x node.coord.y node.coord.z
    end
    println(f)

    # Write connectivities
    println(f, "CELLS ", ncells, " ", nconns)
    for cell in mesh.elems
        print(f, length(cell.nodes), " ")
        for node in cell.nodes
            print(f, node.id-1, " ")
        end
        println(f)
    end
    println(f)

    # Write elem types
    println(f, "CELL_TYPES ", ncells)
    for cell in mesh.elems
        println(f, Int(cell.shape.vtk_type))
    end
    println(f)

    has_node_data = !isempty(mesh.node_data)
    has_elem_data = !isempty(mesh.elem_data)

    # Write node data
    if has_node_data
        println(f, "POINT_DATA ", npoints)
        for (field,D) in mesh.node_data
            isempty(D) && continue
            isfloat = eltype(D)<:AbstractFloat
            dtype = isfloat ? "float64" : "int"
            ncomps = size(D,2)
            if ncomps==1
                println(f, "SCALARS $field $dtype 1")
                println(f, "LOOKUP_TABLE default")
            else
                println(f, "VECTORS ", "$field $dtype")
            end
            for i in 1:npoints
                for j in 1:ncomps
                    if isfloat
                        @printf f "%23.10e" Float32(D[i,j])
                    else
                        # @show 30
                        @printf f "%10d" D[i,j]
                        # @show 40
                    end
                end
            end
            println(f)
        end
    end

    # Write cell data

    if has_elem_data
        # any( contains(field, "tag-") for field in keys(mesh.elem_data) ) && warn("save: Skipping tag string while saving mesh in .vtk legacy format")

        println(f, "CELL_DATA ", ncells)
        for (field,D) in mesh.elem_data
            isempty(D) && continue
            contains(field, "tag-") && continue # skip encoded tags

            isfloat = eltype(D)<:AbstractFloat

            dtype = isfloat ? "float64" : "int"
            ncomps = size(D,2)
            if ncomps==1
                println(f, "SCALARS $field $dtype 1")
                println(f, "LOOKUP_TABLE default")
            else
                println(f, "VECTORS ", "$field $dtype")
            end
            for i in 1:ncells
                for j in 1:ncomps
                    if isfloat
                        @printf f "%23.10e" Float32(D[i,j])
                    else
                        @printf f "%10d" D[i,j]
                    end
                end
            end
            println(f)
        end
    end

    close(f)

    return nothing
end


function get_array_node!(array::AbstractArray, name::String, compressed, buf)
    dtype = string(eltype(array))
    ncomps = size(array,2)
    if compressed # appends compressed data to buf
        xdata = XmlElement("DataArray", attributes=("type"=>dtype, "Name"=>name, "NumberOfComponents"=>"$ncomps", "format"=>"appended", "offset"=>"$(position(buf))"))
        level = 4 # compression level
        inipos = position(buf)

        array = collect(transpose(array))

        # temporary header
        write(buf, UInt64(1), UInt64(0), UInt64(0), UInt64(0))
        arr_size = length(array)*sizeof(eltype(array))

        # write compressed array
        zbuf = ZlibCompressorStream(buf, level=level)
        write(zbuf, array)
        write(zbuf, TranscodingStreams.TOKEN_END)
        flush(zbuf)
        TranscodingStreams.finalize(zbuf.codec) # finalize codec

        # rewrite header
        endpos = position(buf)
        comp_arr_size = endpos - inipos - 4*sizeof(UInt64(0)) # considering header size
        seek(buf, inipos)
        write(buf, UInt64(1), UInt64(arr_size), UInt64(arr_size), UInt64(comp_arr_size))
        seek(buf, endpos)
    else
        xdata = XmlElement("DataArray", attributes=("type"=>dtype, "Name"=>name, "NumberOfComponents"=>"$ncomps", "format"=>"ascii"))
        isfloat = eltype(array)<:AbstractFloat
        io = IOBuffer()
        nrows = size(array,1)
        for i in 1:nrows
            for j in 1:ncomps
                if isfloat
                    @printf io "%20.10e" Float32(array[i,j])
                else
                    print(io, array[i,j], "  ")
                end
            end
            i<nrows && print(io, "\n")
        end
        xdata.content = String(take!(io))
    end
    return xdata
end


function save_vtu(mesh::AbstractDomain, filename::String; desc::String="", compress=false)

    npoints = length(mesh.nodes)
    ncells  = length(mesh.elems)
    root = XmlElement("VTKFile", attributes=("type"=>"UnstructuredGrid", "version"=>"1.0",  "byte_order"=>"LittleEndian", "header_type"=>"UInt64", "compressor"=>"vtkZLibDataCompressor"))
    ugrid = XmlElement("UnstructuredGrid")
    piece = XmlElement("Piece", attributes=("NumberOfPoints"=>"$npoints", "NumberOfCells"=>"$ncells"))
    addchild!(ugrid, piece)
    addchild!(root, ugrid)

    buf = IOBuffer() # only for compressed data

    # Write coordinates
    xpoints = XmlElement("Points")
    coords = Float64[ node.coord[i] for node in mesh.nodes, i in 1:3 ]

    addchild!(xpoints, get_array_node!(coords, "Points", compress, buf))
    addchild!(piece, xpoints)

    xcells = XmlElement("Cells")

    # Write connectivities
    conn = Int32[]
    for cell in mesh.elems
        for node in cell.nodes
            push!(conn, node.id-1)
        end
    end
    addchild!(xcells, get_array_node!(conn, "connectivity", compress, buf))

    # Write offset
    offsets = Int32[]
    offset = 0
    for cell in mesh.elems
        offset += length(cell.nodes)
        push!(offsets, offset)
    end
    # append_compressed_array!(buf, offsets, level)
    addchild!(xcells, get_array_node!(offsets, "offsets", compress, buf))

    # Write cell types
    types = Int32[]
    for cell in mesh.elems
        if cell.role in (:contact, :cohesive, :line_interface, :tip)
            push!(types, Int32(VTK_POLY_VERTEX))
        else
            push!(types, Int32(cell.shape.vtk_type))
        end
    end
    addchild!(xcells, get_array_node!(types, "types", compress, buf))
    addchild!(piece, xcells)

    # Write node data
    has_node_data = !isempty(mesh.node_data)
    if has_node_data
        xpointdata = XmlElement("PointData")
        for (field, D) in mesh.node_data
            isempty(D) && continue
            addchild!(xpointdata, get_array_node!(D, field, compress, buf))
        end
        addchild!(piece, xpointdata)
        # push!(piece.children, xpointdata)
    end

    # Write cell data
    has_elem_data = !isempty(mesh.elem_data)
    if has_elem_data
        xcelldata = XmlElement("CellData")
        for (field, D) in mesh.elem_data
            isempty(D) && continue
            addchild!(xcelldata, get_array_node!(D, field, compress, buf))
        end
        addchild!(piece, xcelldata)
    end

    # Compression
    if compress
        xappended = XmlElement("AppendedData", attributes=("encoding"=>"raw",))
        xappended.content = take!(buf)
        addchild!(root, xappended)
    end

    fileatts = ("version"=>"1.0", "encoding"=>"UTF-8")
    doc = XmlDocument(fileatts)
    addchild!(doc, XmlComment(desc))
    addchild!(doc, root)
    save(doc, filename)
end


function add_extra_fields(mesh::AbstractDomain)
    ncells = length(mesh.elems)

    # Add one-based node-id, elem-id and cell type
    mesh.node_data["node-id"]   = Int[ node.id for node in mesh.nodes ]
    mesh.elem_data["elem-id"]   = Int[ elem.id for elem in mesh.elems ]
    mesh.elem_data["cell-type"] = [ cell.role in (:contact, :cohesive, :line_interface) ? Int32(VTK_POLY_VERTEX) : Int32(cell.shape.vtk_type) for cell in mesh.elems ]

    
    # Add field for interface elements
    interface_types = Set( c.role for c in mesh.elems if c.role in (:contact, :cohesive, :line_interface) )
    if length(interface_types)>0
        interface_data = zeros(Int, ncells, 5) # role, vtk_type, npoints, first linked cell, second linked cell
        mesh.elem_data["interface-data"] = interface_data

        for (i,cell) in enumerate(mesh.elems)
            if cell.role in (:contact, :cohesive)
                nlayers = 2 # Todo: Fix for 3 layers
                interface_data[i,1] = cell.role==:contact ? 1 : 2
                interface_data[i,2] = Int(cell.shape.vtk_type)
                interface_data[i,3] = length(cell.nodes)
                interface_data[i,4] = cell.couplings[1].id
                if length(cell.couplings)==2 # only for cohesive elements (boundary joints have only one linked cell)
                    interface_data[i,5] = cell.couplings[2].id
                end
            elseif cell.role==:line_interface
                interface_data[i,1] = 3
                interface_data[i,2] = Int(cell.shape.vtk_type)
                interface_data[i,3] = cell.shape.npoints
                interface_data[i,4] = cell.couplings[1].id
                interface_data[i,5] = cell.couplings[2].id
            elseif cell.role==:tip
                interface_data[i,1] = 4
                interface_data[i,2] = Int(cell.shape.vtk_type)
                interface_data[i,3] = cell.shape.npoints
                interface_data[i,4] = cell.couplings[1].id
                interface_data[i,5] = cell.couplings[2].id
            end
        end
    end

    # Add field for inset nodes
    # if :line_interface in interface_types
    #     for i in 1:ncells
    #         cell = mesh.elems[i]
    #         if cell.role==:line_interface
    #             interface_data[i,1] = 3
    #             interface_data[i,2] = Int(cell.shape.vtk_type)
    #             interface_data[i,3] = cell.shape.npoints
    #             interface_data[i,4] = cell.couplings[1].id
    #             interface_data[i,5] = cell.couplings[2].id
    #         end
    #     end
    # end

    # Add field for embedded elements
    if any( c.role==:line && length(c.couplings)>0 for c in mesh.elems )
        embedded_data = zeros(Int, ncells) # host cell id
        for i in 1:ncells
            cell = mesh.elems[i]
            if cell.role==:line && length(cell.couplings)>0
                embedded_data[i] = cell.couplings[1].id
            end
        end
        mesh.elem_data["embedded-data"] = embedded_data
    end

    # Add two UInt64 fields to enconde the tag of maximum 16 UTF units
    tags = collect(Set(elem.tag for elem in mesh.elems))
    if length(tags)>1 || tags[1] != ""
        tag_d = Dict( tag=>i for (i,tag) in enumerate(tags) )
        parts = Tuple{String,String}[]
        for tag in tags
            length(codeunits(tag))>16  && error("Mesh: tag '$tag' too long. Max length is 16 UTF units.")
            t1, t2 = safe_string_cut(tag, 8)
            push!(parts, (t1, t2))
        end

        Ts1 = zeros(UInt, length(mesh.elems))
        Ts2 = zeros(UInt, length(mesh.elems))
        for (i,elem) in enumerate(mesh.elems)
            t1, t2 = parts[tag_d[elem.tag]]
            Ts1[i] = encode_string_to_uint64(t1)
            Ts2[i] = encode_string_to_uint64(t2)
        end

        T = Int[ tag_d[elem.tag]-1 for elem in mesh.elems ]

        mesh.elem_data["tag-s1"] = Ts1
        mesh.elem_data["tag-s2"] = Ts2
        mesh.elem_data["tag"]    = T
    end
end


"""
    save(mesh, filename, quiet=true)

Saves a mesh object into a file. Available formats are vtu and vtk.
"""
function save(mesh::AbstractDomain, filename::String; compress=true, quiet=false)
    formats = (".vtk", ".vtu")
    _, format = splitext(filename)
    format in formats || error("save: Cannot save $(typeof(mesh)) to $filename. Available formats are $formats.")

    add_extra_fields(mesh)
    desc = "File generated by Serendip Finite Element Code"

    if format==".vtk"
        save_vtk(mesh, filename, desc=desc)
    elseif format==".vtu"
        save_vtu(mesh, filename, desc=desc, compress=compress)
    end
    quiet || printstyled( "  file $filename saved \e[K \n", color=:cyan)
    return nothing
end


function read_vtk(filename::String)
    # read nodal information
    alltext = read(filename, String)
    data    = split(alltext)

    npoints = 0
    ncells  = 0
    coords  = zeros(0,0)
    connects = Vector{Int}[]
    cell_types = Int[]

    node_data = OrderedDict{String,Array}()
    elem_data = OrderedDict{String,Array}()

    reading_node_data = false
    reading_elem_data  = false

    TYPES = Dict("float32"=>Float32, "float64"=>Float64, "int"=>Int64)

    idx = 1
    while idx<=length(data)
        if data[idx] == "DATASET"
            gridtype = data[idx+1]
            gridtype == "UNSTRUCTURED_GRID" || error("load_VTK_unstructured_grid: this reader only support files of VTK UNSTRUCTURED_GRID")
        end

        # read coords
        if data[idx] == "POINTS"
            npoints = parse(Int64, data[idx+1]) # read number of coords
            coords  = zeros(npoints,3)
            idx += 2
            for i in 1:npoints
                coords[i,1] = parse(Float64, data[idx+1])
                coords[i,2] = parse(Float64, data[idx+2])
                coords[i,3] = parse(Float64, data[idx+3])
                idx += 3
            end
        end

        # read cells connectivities
        if data[idx] == "CELLS"
            ncells = parse(Int64, data[idx+1])
            ncdata = parse(Int64, data[idx+2])
            idx += 2

            connects = Vector{Int}[]
            for i in 1:ncells
                npts = parse(Int64, data[idx+1])
                idx += 1
                conn = Int[]
                for j in 1:npts
                    idx += 1
                    id = parse(Int64, data[idx]) + 1
                    push!(conn, id)
                end
                push!(connects, conn)
            end
        end

        # read type of cells
        if data[idx] == "CELL_TYPES"
            idx += 1
            cell_types = Int[]
            for i in 1:ncells
                idx += 1
                vtk_shape = parse(Int64, data[idx])
                push!(cell_types, vtk_shape)
            end
        end

        if data[idx] == "POINT_DATA"
            idx += 1
            reading_node_data = true
            reading_elem_data  = false
        end

        if data[idx] == "CELL_DATA"
            idx += 1
            reading_elem_data  = true
            reading_node_data = false
        end

        if data[idx] == "VECTORS" && reading_node_data
            label = data[idx+1]
            ty = data[idx+2]
            dtype = TYPES[ty]
            idx += 2
            vectors = zeros(dtype, npoints,3)
            for i in 1:npoints
                vectors[i,1] = parse(dtype, data[idx+1])
                vectors[i,2] = parse(dtype, data[idx+2])
                vectors[i,3] = parse(dtype, data[idx+3])
                idx += 3
            end
            node_data[label] = vectors
        end

        if data[idx] == "VECTORS" && reading_elem_data
            label = data[idx+1]
            ty = data[idx+2]
            dtype = TYPES[ty]
            idx += 2
            vectors = zeros(dtype, ncells,3)
            for i in 1:ncells
                vectors[i,1] = parse(dtype, data[idx+1])
                vectors[i,2] = parse(dtype, data[idx+2])
                vectors[i,3] = parse(dtype, data[idx+3])
                idx += 3
            end
            elem_data[label] = vectors
        end

        if data[idx] == "SCALARS" && reading_node_data
            label = data[idx+1]
            ty = data[idx+2]
            idx += 5
            scalars = zeros(TYPES[ty], npoints)
            for i in 1:npoints
                idx += 1
                scalars[i] = parse(TYPES[ty], data[idx])
            end
            node_data[label] = scalars
        end

        if data[idx] == "SCALARS" && reading_elem_data
            label = data[idx+1]
            ty = data[idx+2]
            idx += 5
            scalars = zeros(TYPES[ty], ncells)
            for i in 1:ncells
                idx += 1
                scalars[i] = parse(TYPES[ty], data[idx])
            end
            elem_data[label] = scalars
        end

        idx += 1

    end

    return Mesh(coords, connects, cell_types, node_data, elem_data)
end


function read_vtu(filename::String)
    doc = XmlDocument(filename)

    # check if file is compressed
    if getchild(doc.root, "AppendedData")!==nothing
        raw = Vector{UInt8}(strip(getchild(doc.root, "AppendedData").content[2:end])) # remove first character _
        nodes = getallchildren(doc.root)

        # update content on xml nodes with offset att
        for node in nodes
            node isa XmlElement || continue
            node.name == "DataArray" && haskey(node.attributes, "offset") || continue
            offset = parse(Int, node.attributes["offset"])

            first  = offset + 1
            last   = offset + 4*sizeof(UInt64(0))

            header = Int.(reinterpret(UInt64, raw[first:last]))
            len    = header[4]

            first = offset + 4*sizeof(UInt64(0)) + 1
            last  = first + len - 1

            # fill xml nodes with decompressed data
            node.content = transcode(ZlibDecompressor, raw[first:last]) # decompression
        end
    end

    piece   = doc.root("UnstructuredGrid", "Piece")
    npoints = parse(Int, piece.attributes["NumberOfPoints"])
    ncells  = parse(Int, piece.attributes["NumberOfCells"])
    coords  = get_array(piece("Points", "DataArray"))

    xcells     = piece("Cells")
    conn       = get_array(xcells["Name"=>"connectivity"][1]) .+ 1
    offsets    = get_array(xcells["Name"=>"offsets"][1])
    cell_types = get_array(xcells["Name"=>"types"][1])

    connects = Vector{Int}[]
    pos = 1
    for off in offsets
        push!(connects, conn[pos:off])
        pos = off+1
    end

    node_data = OrderedDict{String,Array}()
    elem_data = OrderedDict{String,Array}()

    xpointdata = piece("PointData")
    if xpointdata!==nothing
        for array_data in xpointdata.children
            label = array_data.attributes["Name"]
            # label in ("node-id",) && continue # skip extra fields
            node_data[label] = get_array(array_data)
        end
    end

    xcelldata = piece("CellData")
    if xcelldata!==nothing
        for array_data in xcelldata.children
            label = array_data.attributes["Name"]
            # label in ("elem-id", "tag", "tag-s1", "tag-s2", "cell-type", "interface-data", "embedded-data") && continue # skip extra fields
            elem_data[label] = get_array(array_data)
        end
    end

    return Mesh(coords, connects, cell_types, node_data, elem_data)
end


function get_array(data_array::XmlElement)
    TYPES = Dict("Float32"=>Float32, "Float64"=>Float64, "Int32"=>Int32, "Int64"=>Int64, "UInt64"=>UInt64)

    ncomps = haskey(data_array.attributes, "NumberOfComponents") ? parse(Int, data_array.attributes["NumberOfComponents"]) : 1
    dtype  = TYPES[data_array.attributes["type"]]
    isbin  = data_array.content isa Vector{UInt8}
    if isbin
        if ncomps==1
            return reinterpret(dtype, data_array.content)
        else
            return transpose(reshape(reinterpret(dtype, data_array.content), ncomps, :))
        end
    else # string
        if ncomps==1
            return parse.(dtype, split(data_array.content))
        else
            return transpose(reshape(parse.(dtype, split(data_array.content)), ncomps, :))
        end
    end
end


# Setting a Mesh object
function Mesh(coords, connects, vtk_types, node_data, elem_data)

    npoints = size(coords,1)
    ncells  = length(connects)

    # Setting points
    nodes = Node[]
    for i in 1:npoints
        X = coords[i,:]
        node = Node(X)
        node.id = i
        push!(nodes, node)
    end

    # Mesh object
    ndim = get_ndim(nodes)

    # check for exceptional 3d cases
    if haskey(node_data, "rx") || haskey(node_data, "ry")
        ndim = 3
    end

    mesh = Mesh(ndim)
    mesh.nodes = nodes

    # Setting cells
    has_polyvertex = false

    for i in 1:ncells
        conn = mesh.nodes[ connects[i] ]
        vtk_shape = VTKCellType(vtk_types[i])
        if vtk_shape == VTK_POLY_VERTEX
            shape = POLYVERTEX
            has_polyvertex = true
            role = :undefined # to be fixed based on interface_data field
        else
            if vtk_shape==VTK_POLYGON
                shape = QUAD12
            else
                shape = Vtk2CellShape_dict[vtk_shape]
            end
            
            if shape.ndim==1
                role = :line
            else
                role = :bulk
            end

            if shape.ndim==2 && ndim==3
                role = :surface
            end
        end
        cell  = Cell(shape, role, conn)
        cell.id = i
        push!(mesh.elems, cell)
    end

    # update mesh and get faces and edges
    compute_facets(mesh)

    # Setting data
    mesh.node_data = node_data
    mesh.elem_data  = elem_data

    if haskey(mesh.elem_data, "interface-data")
        interface_data = mesh.elem_data["interface-data"]
        
        # Fix information for 1d and 2d/3d interface elements
        for (i,cell) in enumerate(mesh.elems)
            if cell.shape==POLYVERTEX && interface_data[i,1]>0
                # @show interface_data[i,1]
                idx = interface_data[i,1]
                cell.role  = idx==1 ? :contact : idx==2 ? :cohesive : idx==2 ? :line_interface : :tip
                vtk_shape  = VTKCellType(interface_data[i,2])
                cell.shape = Vtk2CellShape_dict[vtk_shape]

                if cell.role == :line_interface
                    cell.couplings = mesh.elems[interface_data[i, 4:5]]
                    cell.couplings[1].crossed = true # host cell is crossed
                elseif cell.role == :tip
                    cell.couplings = mesh.elems[interface_data[i, 4:5]]
                else
                    rng = interface_data[i,5]>0 ? (4:5) : (4:4)  # boundary contact has only one linked cell
                    cell.couplings = mesh.elems[interface_data[i, rng]]
                end
            end
        end
    end

    # check for remaining elements with role :undefined
    for cell in mesh.elems
        if cell.role==:undefined
            throw(SerendipException("Mesh: found element with undefined role (id: $(cell.id) vtk_shape: $(cell.shape.vtk_type))"))
        end
    end

    # Fix information for embedded elements
    if haskey(mesh.elem_data, "embedded-data")
        embedded_data = mesh.elem_data["embedded-data"]
        for (i,cell) in enumerate(mesh.elems)
            if cell.role==LINE_CELL && embedded_data[i]>0
                cell.embedded = true
                cell.couplings = Cell[ mesh.elems[embedded_data[i]] ]
            end
        end
    end

    # Flip cells
    for cell in mesh.elems
        isinverted(cell) && flip!(cell)
    end

    # Build tag if available
    if haskey(mesh.elem_data, "tag-s1") && haskey(mesh.elem_data, "tag-s2")
        Ts1 = mesh.elem_data["tag-s1"]
        Ts2 = mesh.elem_data["tag-s2"]
        for (i,cell) in enumerate(mesh.elems)
            cell.tag = decode_uint64_to_string(Ts1[i]) * decode_uint64_to_string(Ts2[i])
        end
    end

    synchronize(mesh)

    return mesh

end


"""
    Mesh(filename::String; sort=false, quiet=true) -> Mesh

Load a finite element mesh from file. Supports VTK legacy (`.vtk`) and 
VTK XML unstructured (`.vtu`) formats (preferred).

# Arguments
- `filename::String`: Path to the mesh file. Must have extension `.vtk` or `.vtu`.
- `sort::Bool=false`: If `true`, renumber nodes using a bandwidth-reduction
  algorithm after loading.
- `quiet::Bool=true`: If `true`, suppress console messages during loading.

# Returns
- `mesh::Mesh`: Mesh object containing nodes, elements, faces and edges.

# Example
```julia
mesh = Mesh("model.vtu"; sort=true, quiet=false)
```
"""
function Mesh(filename::String; sort=false, quiet=true)

    formats = (".vtk", ".vtu")

    quiet || printstyled("Mesh loading: filename $filename\n", bold=true, color=:cyan)

    _, format = splitext(filename)

    format in formats || error("Mesh: cannot read format \"$format\". Suitable formats are $formats.")

    if format==".vtk"
        quiet || print("  Reading VTK legacy format...\n")
        mesh = read_vtk(filename)
    elseif format==".vtu"
        quiet || print("  Reading VTU format...\n")
        mesh = read_vtu(filename)
    end

    quiet || printstyled( "  file $filename loaded \e[K \n", color=:cyan)

    # Reorder nodal numbering
    if sort
        quiet || print("  reordering points...\r")
        sort_mesh(mesh)
    end

    if !quiet
        npoints = length(mesh.nodes)
        ncells  = length(mesh.elems)
        nfaces  = length(mesh.faces)
        nedges  = length(mesh.edges)
        println("  ", mesh.ctx.ndim, "d                   ")
        @printf "  %5d points\n" npoints
        @printf "  %5d cells\n" ncells
        @printf "  %5d faces\n" nfaces
        @printf "  %5d surface edges\n" nedges
    end

    return mesh
end
