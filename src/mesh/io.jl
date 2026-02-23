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
        if cell.role in (:contact, :cohesive, :line_interface, :tip)
            println(f, Int(VTK_POLY_VERTEX))
        else
            println(f, Int(cell.shape.vtk_type))
        end
    end
    println(f)

    has_node_data = !isempty(mesh.node_fields)
    has_elem_data = !isempty(mesh.elem_fields)

    # Write node data
    if has_node_data
        println(f, "POINT_DATA ", npoints)
        for (field,D) in mesh.node_fields
            isempty(D) && continue
            isfloat = eltype(D)<:AbstractFloat
            T = isfloat ? "float64" : "int"
            ncomps = size(D,2)
            if ncomps==1
                println(f, "SCALARS $field $T 1")
                println(f, "LOOKUP_TABLE default")
            else
                println(f, "VECTORS ", "$field $T")
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
        # any( contains(field, "tag-") for field in keys(mesh.elem_fields) ) && warn("save: Skipping tag string while saving mesh in .vtk legacy format")

        println(f, "CELL_DATA ", ncells)
        for (field,D) in mesh.elem_fields
            isempty(D) && continue
            contains(field, "tag-") && continue # skip encoded tags

            isfloat = eltype(D)<:AbstractFloat

            T = isfloat ? "float64" : "int"
            ncomps = size(D,2)
            if ncomps==3
                println(f, "VECTORS ", "$field $T")
            else
                println(f, "SCALARS $field $T $ncomps")
                println(f, "LOOKUP_TABLE default")
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
    T = string(eltype(array))
    ncomps = size(array,2)
    if compressed # appends compressed data to buf
        xdata = XmlElement("DataArray", attributes=("type"=>T, "Name"=>name, "NumberOfComponents"=>"$ncomps", "format"=>"appended", "offset"=>"$(position(buf))"))
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
        xdata = XmlElement("DataArray", attributes=("type"=>T, "Name"=>name, "NumberOfComponents"=>"$ncomps", "format"=>"ascii"))
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
    root_atts = compress ?
        ("type"=>"UnstructuredGrid", "version"=>"1.0", "byte_order"=>"LittleEndian", "header_type"=>"UInt64", "compressor"=>"vtkZLibDataCompressor") :
        ("type"=>"UnstructuredGrid", "version"=>"1.0", "byte_order"=>"LittleEndian")
    root = XmlElement("VTKFile", attributes=root_atts)
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
    has_node_data = !isempty(mesh.node_fields)
    if has_node_data
        xpointdata = XmlElement("PointData")
        for (field, D) in mesh.node_fields
            isempty(D) && continue
            addchild!(xpointdata, get_array_node!(D, field, compress, buf))
        end
        addchild!(piece, xpointdata)
        # push!(piece.children, xpointdata)
    end

    # Write cell data
    has_elem_data = !isempty(mesh.elem_fields)
    if has_elem_data
        xcelldata = XmlElement("CellData")
        for (field, D) in mesh.elem_fields
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


function save_json(mesh::AbstractDomain, filename::String; desc::String="")

    # cells_obj = Vector{Dict{String,Any}}(undef, length(cells))

    # for (i, cell) in enumerate(mesh.elems)
    #     out[i] = Dict(
    #         "type" => cell.shape.name,
    #         "connectivity" => [n.id for n in cell.nodes]
    #     )
    # end

    nodes_obj = [ node.coord for node in mesh.nodes ]

    cells_obj = [
        Dict(
            "type" => cell.shape.name,
            "conn" => [n.id for n in cell.nodes]
        ) for cell in mesh.elems
    ]

    faces_obj = [
        Dict(
            "type" => cell.shape.name,
            "conn" => [n.id for n in cell.nodes]
        ) for cell in mesh.faces
    ]

    rounds = (V) -> eltype(V)<:AbstractFloat ? round.(V, sigdigits=4) : V

    node_data_obj = Dict{String,Array}( k => rounds(V) for (k,V) in mesh.node_fields if size(V,2)==1 )
    elem_data_obj = Dict{String,Array}( k => rounds(V) for (k,V) in mesh.elem_fields if size(V,2)==1 )

    mesh_obj = Dict(
        "mesh" => Dict(
            "nodes"    => nodes_obj,
            "elements" => cells_obj,
            "faces"    => faces_obj
        ),
        "nodal_fields"   => node_data_obj,
        "element_fields" => elem_data_obj
    )

    JSON.json(filename, mesh_obj, pretty=true)

    return nothing
end


function add_extra_fields(mesh::AbstractDomain)
    ncells = length(mesh.elems)

    # Add one-based node-id, elem-id and cell type
    mesh.node_fields["node-id"]   = Int[ node.id for node in mesh.nodes ]
    mesh.elem_fields["elem-id"]   = Int[ elem.id for elem in mesh.elems ]
    mesh.elem_fields["cell-type"] = [ cell.role in (:contact, :cohesive, :line_interface) ? Int32(VTK_POLY_VERTEX) : Int32(cell.shape.vtk_type) for cell in mesh.elems ]

    
    # Add field for interface elements
    interface_types = Set( c.role for c in mesh.elems if c.role in (:contact, :cohesive, :line_interface) )
    if length(interface_types)>0
        interface_data = zeros(Int, ncells, 5) # role, vtk_type, npoints, first linked cell, second linked cell
        mesh.elem_fields["interface-data"] = interface_data

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
        mesh.elem_fields["embedded-data"] = embedded_data
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

        mesh.elem_fields["tag-s1"] = Ts1
        mesh.elem_fields["tag-s2"] = Ts2
        mesh.elem_fields["tag"]    = T
    end
end


"""
    save(mesh, filename, quiet=true)

Saves a mesh object into a file. Available formats are vtu and vtk.
"""
function save(mesh::AbstractDomain, filename::String; compress=true, quiet=false)
    formats = (".vtk", ".vtu", ".json")
    _, format = splitext(filename)
    format in formats || error("save: Cannot save $(typeof(mesh)) to $filename. Available formats are $formats.")

    add_extra_fields(mesh)
    desc = "File generated by Serendip Finite Element Code"

    if format==".vtk"
        save_vtk(mesh, filename, desc=desc)
    elseif format==".vtu"
        save_vtu(mesh, filename, desc=desc, compress=compress)
    elseif format==".json"
        save_json(mesh, filename, desc=desc)
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

    node_fields = OrderedDict{String,Array}()
    elem_fields = OrderedDict{String,Array}()

    reading_node_data = false
    reading_elem_data = false

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
            T  = TYPES[ty]
            idx += 2
            vectors = zeros(T, npoints,3)
            for i in 1:npoints
                vectors[i,1] = parse(T, data[idx+1])
                vectors[i,2] = parse(T, data[idx+2])
                vectors[i,3] = parse(T, data[idx+3])
                idx += 3
            end
            node_fields[label] = vectors
        end

        if data[idx] == "VECTORS" && reading_elem_data
            label = data[idx+1]
            ty = data[idx+2]
            T = TYPES[ty]
            idx += 2
            vectors = zeros(T, ncells,3)
            for i in 1:ncells
                vectors[i,1] = parse(T, data[idx+1])
                vectors[i,2] = parse(T, data[idx+2])
                vectors[i,3] = parse(T, data[idx+3])
                idx += 3
            end
            elem_fields[label] = vectors
        end

        if data[idx] == "SCALARS" && reading_node_data
            label = data[idx+1]
            ty = data[idx+2]
            T        = TYPES[ty]
            idx += 5
            scalars = zeros(T, npoints)
            for i in 1:npoints
                idx += 1
                scalars[i] = parse(T, data[idx])
            end
            node_fields[label] = scalars
        end

        if data[idx] == "SCALARS" && reading_elem_data
            label    = data[idx+1]
            ty       = data[idx+2]
            T        = TYPES[ty]
            ncomps   = parse(Int, data[idx+3])
            idx     += 5
            
            if ncomps == 1
                scalars  = zeros(T, ncells)
                for i in 1:ncells
                    idx += 1
                    scalars[i] = parse(T, data[idx])
                end
                elem_fields[label] = scalars
            else
                vectors  = zeros(T, ncells, ncomps)
                for i in 1:ncells
                    for j in 1:ncomps
                        idx += 1
                        vectors[i,j] = parse(T, data[idx])
                    end
                end
                elem_fields[label] = vectors
            end
        end

        idx += 1

    end

    return Mesh(coords, connects, cell_types, node_fields, elem_fields)
end


function read_vtu(filename::String)
    doc = XmlDocument(filename)

    # Decode appended binary blocks (compressed or uncompressed) into DataArray contents.
    xapp = getchild(doc.root, "AppendedData")
    if xapp !== nothing
        # IMPORTANT:
        # For encoding="raw", XML text normalization can alter raw bytes.
        # Extract appended payload directly from file bytes instead of xapp.content.
        file_bytes = read(filename)
        app_tag = findfirst(UInt8[codeunits("<AppendedData")...], file_bytes)
        app_tag === nothing && error("read_vtu: <AppendedData> tag not found in file bytes")

        tag_close = findnext(==(UInt8('>')), file_bytes, first(app_tag))
        tag_close === nothing && error("read_vtu: Invalid <AppendedData> tag (missing '>').")

        marker = findnext(==(UInt8('_')), file_bytes, tag_close + 1)
        marker === nothing && error("read_vtu: AppendedData marker '_' not found.")

        app_end = findfirst(UInt8[codeunits("</AppendedData>")...], file_bytes)
        app_end === nothing && error("read_vtu: </AppendedData> closing tag not found.")

        raw_start = marker + 1
        raw_stop = first(app_end) - 1
        raw_start > raw_stop && error("read_vtu: Empty appended payload.")

        # Do not strip/trim: trailing whitespace bytes can belong to compressed blocks.
        raw = file_bytes[raw_start:raw_stop]
        nodes = getallchildren(doc.root)

        header_type = get(doc.root.attributes, "header_type", "UInt32")
        H = header_type == "UInt64" ? UInt64 :
            header_type == "UInt32" ? UInt32 :
            error("read_vtu: Unsupported header_type=$header_type")
        hbytes = sizeof(H)
        is_compressed = haskey(doc.root.attributes, "compressor")

        # Update XML node contents from appended payload using each DataArray offset.
        for node in nodes
            node isa XmlElement || continue
            node.name == "DataArray" && haskey(node.attributes, "offset") || continue
            offset = parse(Int, node.attributes["offset"])

            first = offset + 1
            first > length(raw) && error("read_vtu: Invalid offset $offset (out of payload bounds)")

            if is_compressed
                # VTK compressed layout:
                # [num_blocks, block_size, last_block_size, compressed_block_sizes...]
                nblocks = Int(reinterpret(H, raw[first:first+hbytes-1])[1])
                nblocks < 1 && error("read_vtu: Invalid number of compressed blocks: $nblocks")

                hlen = (3 + nblocks)*hbytes
                last = offset + hlen
                last > length(raw) && error("read_vtu: Compressed header exceeds payload at offset $offset")

                header = Int.(reinterpret(H, raw[first:last]))
                comp_sizes = header[4:end]
                comp_total = sum(comp_sizes)
                comp_first = offset + hlen + 1
                comp_last = comp_first + comp_total - 1
                comp_last > length(raw) && error("read_vtu: Compressed block exceeds payload at offset $offset")

                if nblocks == 1
                    node.content = transcode(ZlibDecompressor, raw[comp_first:comp_last])
                else
                    out = UInt8[]
                    pos = comp_first
                    for clen in comp_sizes
                        chunk_last = pos + clen - 1
                        chunk_last > comp_last && error("read_vtu: Invalid compressed chunk size at offset $offset")
                        append!(out, transcode(ZlibDecompressor, raw[pos:chunk_last]))
                        pos = chunk_last + 1
                    end
                    node.content = out
                end
            else
                # VTK uncompressed appended layout: [byte_length, payload...]
                last = first + hbytes - 1
                last > length(raw) && error("read_vtu: Uncompressed header exceeds payload at offset $offset")
                len = Int(reinterpret(H, raw[first:last])[1])
                data_first = last + 1
                data_last = data_first + len - 1
                data_last > length(raw) && error("read_vtu: Uncompressed payload exceeds bounds at offset $offset")
                node.content = raw[data_first:data_last]
            end
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

    node_fields = OrderedDict{String,Array}()
    elem_fields = OrderedDict{String,Array}()

    xpointdata = piece("PointData")
    if xpointdata!==nothing
        for array_data in xpointdata.children
            label = array_data.attributes["Name"]
            node_fields[label] = get_array(array_data)
        end
    end

    xcelldata = piece("CellData")
    if xcelldata!==nothing
        for array_data in xcelldata.children
            label = array_data.attributes["Name"]
            elem_fields[label] = get_array(array_data)
        end
    end

    return Mesh(coords, connects, cell_types, node_fields, elem_fields)
end


function get_array(data_array::XmlElement)
    TYPES = Dict("Float32"=>Float32, "Float64"=>Float64, "Int32"=>Int32, "Int64"=>Int64, "UInt64"=>UInt64)

    ncomps = haskey(data_array.attributes, "NumberOfComponents") ? parse(Int, data_array.attributes["NumberOfComponents"]) : 1
    T  = TYPES[data_array.attributes["type"]]
    isbin  = data_array.content isa Vector{UInt8}
    if isbin
        if ncomps==1
            return reinterpret(T, data_array.content)
        else
            return transpose(reshape(reinterpret(T, data_array.content), ncomps, :))
        end
    else # string
        if ncomps==1
            return parse.(T, split(data_array.content))
        else
            return transpose(reshape(parse.(T, split(data_array.content)), ncomps, :))
        end
    end
end


# Setting a Mesh object
function Mesh(coords, connects, vtk_types, node_fields, elem_fields)

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
    if haskey(node_fields, "rx") || haskey(node_fields, "ry")
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
                role = :cont
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
    mesh.node_fields = node_fields
    mesh.elem_fields  = elem_fields

    if haskey(mesh.elem_fields, "interface-data")
        interface_data = mesh.elem_fields["interface-data"]
        
        # Fix information for 1d and 2d/3d interface elements
        for (i,cell) in enumerate(mesh.elems)
            if cell.shape==POLYVERTEX && interface_data[i,1]>0
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
    if haskey(mesh.elem_fields, "embedded-data")
        embedded_data = mesh.elem_fields["embedded-data"]
        for (i,cell) in enumerate(mesh.elems)
            if cell.role==LINE_CELL && embedded_data[i]>0
                cell.embedded = true
                cell.couplings = Cell[ mesh.elems[embedded_data[i]] ]
            end
        end
    end

    # Flip cells
    for cell in mesh.elems
        isinverted(cell) && flip(cell)
    end

    # Build tag if available
    if haskey(mesh.elem_fields, "tag-s1") && haskey(mesh.elem_fields, "tag-s2")
        Ts1 = mesh.elem_fields["tag-s1"]
        Ts2 = mesh.elem_fields["tag-s2"]
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
