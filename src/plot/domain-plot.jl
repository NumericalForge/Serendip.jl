# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


"""
    DomainPlot(mesh; 
        size=(220,150), face_color=:aliceblue, warp=0.0,
        edge_width=0.3, edge_color=:auto, outline_width=0.4,
        line_element_width=0.6,
        field="", limits=[0.0,0.0], field_kind=:auto, field_mult=1.0,
        label="", colormap=:coolwarm, diverging=false,
        colorbar=:right, colorbar_length_factor=0.9, bins=6,
        font="NewComputerModern", font_size=7.0,
        interpolation=:linear,
        azimuth=30, elevation=30, distance=0.0,
        feature_edges=true, view_mode=:surface_with_edges,
        light_vector=[0,0,0],
        node_labels=false,
        axes=:none, axis_labels=String[],
        quiet=false)

Create a customizable domain plot for meshes and FE models.

# Arguments
- `mesh::Union{AbstractDomain,FEModel}`: object to plot.

# Keywords
- `size::Tuple{Int,Int}`: figure size in pt. (1 cm = 28.35 pt).
- `face_color::Symbol`: surface color for area/surface elements.
- `warp::Real`: displacement scale factor for warped views.
- `edge_width::Real`: internal edge width for area/surface cells.
- `edge_color::Union{Tuple,Symbol}`: edge color; `:auto` darkens face color.
- `outline_width::Real`: boundary outline width.
- `line_element_width::Real`: stroke width for line elements (bars/beams).
- `field::AbstractString`: scalar field name for coloring; empty disables.
- `limits::Vector{<:Real}`: field range `[min,max]`; `[0,0]` → auto.
- `field_kind::Symbol`: `:auto | :none | :node | :element`.
- `field_mult::Real`: multiplier applied to field values.
- `label::AbstractString`: colorbar label (alias: `colorbar_label`).
- `colormap::Union{Symbol,Colormap}`: e.g. `:viridis`, `:coolwarm`, or a `Colormap`.
- `diverging::Bool`: center colormap at zero.
- `colorbar::Symbol`: `:none | :right | :bottom`.
- `colorbar_length_factor::Real`: colorbar length scale (> 0).
- `bins::Int`: number of colorbar bins.
- `font::AbstractString`: font family.
- `font_size::Real`: font size (> 0).
- `interpolation::Symbol`: `:constant | :linear | :nonlinear` surface shading.
- `azimuth::Real`: 3D azimuth angle in degrees.
- `elevation::Real`: 3D elevation angle in degrees.
- `distance::Real`: camera distance (≥ 0).
- `feature_edges::Bool`: enhance feature lines (`:outline` forces on).
- `view_mode::Symbol`: `:surface_with_edges | :surface | :wireframe | :outline`.
- `light_vector::Vector{<:Real}`: light direction.
- `node_labels::Bool`: show node ids.
- `axes::Symbol`: axes widget location (see `_axes_widget_locations`).
- `axis_labels::Vector{<:AbstractString}`: custom axis labels; defaults to `["x","y","z"]`.
- `quiet::Bool`: suppress constructor log.

# Notes
- Use `save` to export the chart to a file.

# Example
```julia
plt = DomainPlot(model;
                 field="ux", field_kind=:node, warp=50.0,
                 colormap=:viridis, colorbar=:right, label="uₓ [mm]",
                 view_mode=:surface_with_edges, feature_edges=true)
save(plt, "model_plot.pdf")
```
"""
mutable struct DomainPlot<:Figure
    mesh::AbstractDomain
    canvas::FigureComponent
    colorbar::FigureComponent
    axes::FigureComponent
    nodes::Vector{Node}
    elems::Vector{AbstractCell}
    feature_edges_d::Dict
    values::Vector{Float64}
    outerpad::Float64
    shades::Vector{Float64}
    azimuth::Float64
    elevation::Float64
    distance::Float64

    interpolation::Symbol
    light_vector::Vector{Float64}

    width::Float64
    height::Float64

    edge_width::Float64
    edge_color::Union{Symbol,Tuple}
    outline_width::Float64
    line_elem_width::Float64
    face_color::Tuple
    field::String
    field_kind::Symbol
    limits::Vector{Float64}
    field_mult::Float64
    warp::Float64
    label::AbstractString
    colormap::Colormap
    diverging::Bool
    colorbar_location::Symbol
    colorbar_length_factor::Float64
    bins::Int
    font::String
    font_size::Float64
    show_feature_edges::Bool
    view_mode::Symbol
    node_labels::Bool
    axes_loc::Symbol
    axis_labels::Vector{AbstractString}
    quiet::Bool

    function DomainPlot(mesh;
        size::Tuple{Int,Int}=(220,150),
        face_color::Symbol=:aliceblue,
        warp::Real=0.0,
        edge_width::Real=0.3,
        edge_color::Symbol=:auto,
        outline_width::Real=0.4,
        line_elem_width::Real=0.6,
        field::AbstractString="",
        field_kind::Symbol=:auto,
        limits::Vector{Float64}=[0.0,0.0],
        field_mult::Real=1.0,
        label::AbstractString="",
        colormap::Union{Symbol,Colormap}=:coolwarm,
        diverging::Bool=false,
        colorbar::Symbol=:right,
        colorbar_length_factor::Real=0.9,
        bins::Int=6,
        font::AbstractString="NewComputerModern",
        font_size::Real=7.0,
        interpolation::Symbol=:linear,
        azimuth::Real=30.0,
        elevation::Real=30.0,
        distance::Real=0.0,
        feature_edges::Bool=true,
        view_mode::Symbol=:surface_with_edges,
        light_vector::Vector{Float64}=[0.0,0.0,0.0],
        node_labels::Bool=false,
        axes::Symbol=:none,
        axis_labels::Vector{<:AbstractString}=String[],
        quiet::Bool=false,
    )

        @check face_color in keys(_colors_dict) "face_color must be one of: $(collect(keys(_colors_dict))). Got $face_color"
        @check edge_color == :auto || edge_color in keys(_colors_dict) "edge_color must be :auto or one of: $(collect(keys(_colors_dict))). Got $edge_color"
        @check font in _available_fonts "font must be one of: $(_available_fonts). Got $(repr(font))"

        canvas = Canvas()
        the_colorbar = Colorbar()
        width, height = size
        axes_loc    = axes
        axis_labels = length(axis_labels)==0 ? [L"x", L"y", L"z"] : axis_labels
        axes_widget = AxisWidget(location=axes_loc, labels=axis_labels, font=font, font_size=font_size, azimuth=azimuth, elevation=elevation, arrow_length=20)

        mesh = mesh
        nodes = []
        elems = []
        feature_edges_d = Dict()
        values = []
        shades = []
        outerpad = 0.0

        face_color   = _colors_dict[face_color]
        field        = string(field)
        field_kind   = string(field) == "" ? :none : field_kind
        limits       = collect(limits)
        field_mult   = field_mult
        warp         = warp
        label        = label

        colormap = colormap isa Symbol ? Colormap(colormap) : colormap

        colorbar_location   = colorbar
        colorbar_length_factor = colorbar_length_factor

        font               = font
        font_size          = font_size
        azimuth             = azimuth
        elevation          = elevation
        distance           = distance
        show_feature_edges = (feature_edges || view_mode==:outline) && view_mode != :wireframe

        if !quiet
            printstyled("Domain plot\n", bold=true, color=:cyan)
            println("  size: $(width) x $(height) pt")
            println("  view mode: $(view_mode)")
            field != "" && println("  field: $(field)")
        end

        return new(mesh, canvas, the_colorbar, axes_widget, nodes, elems, 
            feature_edges_d, values, outerpad, shades,
            azimuth, elevation, distance,
            interpolation, light_vector,
            width, height,
            edge_width, edge_color, outline_width, 
            line_elem_width, face_color,
            field, field_kind, limits, field_mult,
            warp, label, colormap, diverging,
            colorbar_location, colorbar_length_factor,
            bins, font, font_size,
            show_feature_edges, view_mode,
            node_labels, axes_loc, axis_labels,
            quiet
        )
    end
end


function bezier_points(edge::AbstractCell)
    p1 = edge.nodes[1].coord[1:2]
    p4 = edge.nodes[2].coord[1:2]
    ξ1 = -1/3
    ξ2 = +1/3
    C = get_coords(edge.nodes, 2)
    p2 = C'*edge.shape.func([ξ1])
    p3 = C'*edge.shape.func([ξ2])

    cp2 = 1/6*(-5*p1+18*p2- 9*p3+2*p4)
    cp3 = 1/6*( 2*p1- 9*p2+18*p3-5*p4)

    return [p1, cp2, cp3, p4]
end


function project_to_2d!(nodes, azimuth, elevation, distance)
    # Find bounding box
    xmin, xmax = extrema( node.coord[1] for node in nodes)
    ymin, ymax = extrema( node.coord[2] for node in nodes)
    zmin, zmax = extrema( node.coord[3] for node in nodes)
    reflength = max(xmax-xmin, ymax-ymin, zmax-zmin)

    # Centralize
    center = 0.5*Vec3(xmin+xmax, ymin+ymax, zmin+zmax)
    for node in nodes
        node.coord = node.coord - center
    end

    # Rotation around z axis
    θ = -azimuth*pi/180
    R = Quaternion(cos(θ/2), 0, 0, sin(θ/2))
    for node in nodes
        node.coord = (R*node.coord*conj(R))[2:4]
    end

    # Rotation around y axis
    θ = elevation*pi/180
    R = Quaternion(cos(θ/2), 0, sin(θ/2), 0)
    for node in nodes
        node.coord = (R*node.coord*conj(R))[2:4]
    end

    # Set projection values
    distance==0 && (distance=reflength*3)
    distance = max(distance, reflength)
    focal_length = 0.1*distance

    # Make projection
    for node in nodes
        x´ = node.coord[1]
        y′ = node.coord[2]*focal_length/(distance-x´)
        z′ = node.coord[3]*focal_length/(distance-x´)
        node.coord = Vec3(y′, z′, distance-x´)
    end

    xmin, xmax = extrema( node.coord[1] for node in nodes)
    ymin, ymax = extrema( node.coord[2] for node in nodes)

    # normalize
    l = max(xmax-xmin, ymax-ymin)
    for node in nodes
        x = (node.coord[1]-xmin)/l
        y = (node.coord[2]-ymin)/l
        node.coord = Vec3(x, y, node.coord[3])
    end
end


function configure!(mplot::DomainPlot)
    mesh = copy(mplot.mesh)
    ndim = mesh.ctx.ndim

    # Change coords if warping
    if mplot.warp>0.0
        U = get(mesh.node_data, "U", nothing)
        if U === nothing
            alert("DomainPlot: Vector field U not found for warping.")
        else
            for node in mesh.nodes
                node.coord = node.coord + mplot.warp*U[node.id,:]
            end
        end
    end

    active_elems = select(mesh.elems, :active)
    if ndim==2
        areacells = [ elem for elem in active_elems if elem.role==:bulk ]
        linecells = [ cell for cell in active_elems if cell.role==:line ]
        feature_edges = get_outer_facets(areacells)

        elems    = [ areacells; linecells ]
        nodes    = get_nodes(elems)
    else
        # get surface cells and update
        volcells  = [ elem for elem in active_elems if elem.role==:bulk && elem.shape.ndim==3 ]
        areacells = [ elem for elem in active_elems if elem.role==:bulk && elem.shape.ndim==2 ]
        linecells = [ cell for cell in active_elems if cell.role==:line ]
        surfcells = get_outer_facets(volcells)
        feature_edges = get_feature_edges(surfcells)

        for c in surfcells
            c.id = c.owner.id # useful to get values but not shades
        end

        elems    = [ surfcells; areacells; linecells ]

        # unique by id
        nodes = unique( n -> n.id, [ node for elem in elems for node in elem.nodes ] )
    end

    node_data = mesh.node_data
    elem_data = mesh.elem_data

    # populate feature_edges_d
    for edge in feature_edges
        node_idxs = sort([ node.id for node in edge.nodes ])
        mplot.feature_edges_d[node_idxs] = edge
    end

    # 3D -> 2D projection
    if mesh.ctx.ndim==3
        # compute shades (before 2d projection)
        V = Vec3( cosd(mplot.elevation)*cosd(mplot.azimuth), cosd(mplot.elevation)*sind(mplot.azimuth), sind(mplot.elevation) ) # observer vector
        norm(mplot.light_vector)==0 && (mplot.light_vector = V)
        L = mplot.light_vector
        # shades = zeros(length(mesh.elems))
        shades = zeros(length(elems))
        for (i,elem) in enumerate(elems)
            elem.role==:surface || continue
            N = get_facet_normal(elem)
            dot(N,L)<0 && (N = -N)
            R = normalize(2*N*dot(L,N) - L)  # reflection vector
            dot(V,R)<0 && (R = -R)

            Ia = 0.68 # ambient
            Id = 0.34 # diffuse
            Is = 0.20 # specular
            sh = 2.00 # shininess

            shades[i] = Ia + Id*max(0,dot(L,N)) + Is*max(0,dot(V,R))^sh
        end
        mplot.shades = shades

        # compute projection
        project_to_2d!(nodes, mplot.azimuth, mplot.elevation, mplot.distance)
        # zmin, zmax = extrema(node.coord[3] for node in nodes)

        # distances = [ sum(node.coord[3] for node in elem.nodes)/length(elem.nodes)  for elem in mesh.elems ]
        # distances = [ minimum(node.coord[3] for node in elem.nodes) for elem in elems ]
        distances = [ 0.9*sum(node.coord[3] for node in elem.nodes)/length(elem.nodes) + 0.1*minimum(node.coord[3] for node in elem.nodes) for elem in elems ]
        perm = sortperm(distances, rev=true)
        elems = elems[perm]
        mplot.shades = mplot.shades[perm]
    else
        # sort elems by area (show largest first)
        areas = [ elem.role==:bulk ? cell_extent(elem) : 0 for elem in elems ]
        perm = sortperm(areas, rev=true)
        elems = elems[perm]
    end

    # Field
    has_field = mplot.field != ""

    if has_field
        field = string(mplot.field)
        mplot.label == ""  && (mplot.label = field)

        # available fields
        fields = [ string(field) for field in Iterators.flatten([keys(node_data), keys(elem_data)]) ]
        field in fields || error("mplot: field $field not found. Available fields are: $(join(fields, ", ", " and "))")

        if mplot.field_kind == :auto
            if haskey(node_data, field)
                mplot.field_kind = :node
            elseif haskey(elem_data, field)
                mplot.field_kind = :element
            else
                # fields = [ string(field) for field in Iterators.flatten([keys(node_data), keys(elem_data)]) ]
                fields = string.([ collect(keys(node_data)); collect(keys(elem_data)) ])
                error("mplot: field $field not found. Available fields are: $(join(fields, ", ", " and "))")
            end
        end

        if mplot.field_kind == :element
            fvals = elem_data[field].*mplot.field_mult
            fmax = maximum(fvals)
            fmin = minimum(fvals)
        else
            fvals = node_data[field].*mplot.field_mult
            fmax = maximum(fvals[node.id] for node in nodes)
            fmin = minimum(fvals[node.id] for node in nodes)
        end
        if fmin==fmax
            fmin -= 1
            fmax += 1
        end

        # Colormap
        mplot.values = fvals
        if mplot.limits==[0.0, 0.0]
            mplot.limits = [fmin, fmax]
        else
            fmin, fmax = mplot.limits
        end
        mplot.colormap = resize(mplot.colormap, fmin, fmax, diverging=mplot.diverging)

    end

    mplot.nodes = nodes
    mplot.elems = elems

    # Canvas
    # mplot.canvas = Canvas()
    width, height = mplot.width, mplot.height
    mplot.outerpad = 0.01*min(width, height)

    rpane = 0.0
    bpane = 0.0

    # Colorbar
    if has_field
        mplot.colorbar = Colorbar(;
            location      = mplot.colorbar_location,
            label         = mplot.label,
            length_factor = mplot.colorbar_length_factor,
            colormap      = mplot.colormap,
            limits        = mplot.limits,
            font_size     = mplot.font_size,
            font          = mplot.font,
            bins          = mplot.bins,
        )
        configure!(mplot, mplot.colorbar)

        if mplot.colorbar.location==:right
            rpane = mplot.colorbar.width
        elseif mplot.colorbar.location==:bottom
            bpane = mplot.colorbar.height
        end
    end

    # Canvas box
    canvas = mplot.canvas

    canvas.width = width - rpane - 2*mplot.outerpad
    canvas.height = height - bpane - 2*mplot.outerpad
    Xmin = mplot.outerpad
    Ymin = mplot.outerpad
    Xmax = width - rpane - mplot.outerpad
    Ymax = height - bpane - mplot.outerpad
    canvas.box = [ Xmin, Ymin, Xmax, Ymax ] # in user coordinates

    xmin, xmax = extrema( node.coord[1] for node in mplot.nodes)
    ymin, ymax = extrema( node.coord[2] for node in mplot.nodes)

    ratio = min((Xmax-Xmin)/(xmax-xmin), (Ymax-Ymin)/(ymax-ymin) )
    dX = 0.5*((Xmax-Xmin)-ratio*(xmax-xmin))
    dY = 0.5*((Ymax-Ymin)-ratio*(ymax-ymin))

    xmin = xmin - dX/ratio
    xmax = xmax + dX/ratio
    ymin = ymin - dY/ratio
    ymax = ymax + dY/ratio

    mplot.canvas.limits = [xmin, ymin, xmax, ymax] # in data coordinates

    # Axes widget
    if mplot.axes_loc != :none
        mplot.axes = AxisWidget(
            location  = mplot.axes_loc,
            font_size = mplot.font_size,
            font      = mplot.font,
            azimuth    = mplot.azimuth,
            elevation = mplot.elevation,
            labels    = mplot.axis_labels[1:ndim],
        )
        configure!(mplot.axes)
    end
end


function iscounterclockwise(points::Array{Node,1})
    val = 0.0
    n = length(points)
    for i in 1:n
        x1, y1 = points[i].coord
        x2, y2 = points[mod1(i+1, n)].coord
        val += (x1*y2) - (x2*y1)
    end
    return val > 0
end


function draw!(mplot::DomainPlot, ctx::CairoContext)
    # Cairo.push_group(ctx)
    # set_operator(ctx, Cairo.OPERATOR_SOURCE)

    set_line_join(ctx, Cairo.CAIRO_LINE_JOIN_ROUND)
    set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)
    font = get_font(mplot.font)
    select_font_face(ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
    set_font_size(ctx, mplot.font_size)

    # has_field = mplot.field != ""
    # is_nodal_field = has_field && haskey(mplot.mesh.node_data, mplot.field)

    Xmin, Ymin, Xmax, Ymax = mplot.canvas.box
    xmin, ymin, xmax, ymax = mplot.canvas.limits

    ratio = (Xmax-Xmin)/(xmax-xmin) # == (Ymax-Ymin)/(ymax-ymin)
    set_matrix(ctx, CairoMatrix([ratio, 0, 0, -ratio, Xmin-xmin*ratio, Ymax+ymin*ratio]...))

    # Draw elements
    for (i,elem) in enumerate(mplot.elems)

        elem.role in ( :line, :bulk, :surface ) || continue

        # culling back faces
        if elem.owner !== nothing && mplot.view_mode in (:surface, :surface_with_edges)
            if elem.owner.shape.ndim==3
                nnodes_basic = elem.shape.base_shape.npoints
                !iscounterclockwise( elem.nodes[1:nnodes_basic] ) && continue
            end
        end

        if elem.role==:line
            x, y = elem.nodes[1].coord
            move_to(ctx, x, y)
            color = Vec3(0.8, 0.2, 0.1)
            # if has_field && !is_nodal_field
            if mplot.field_kind == :element
                color = mplot.colormap(mplot.values[elem.id])
            end

            pts = bezier_points(elem)
            curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            set_line_width(ctx, mplot.line_elem_width)
            set_source_rgb(ctx, color...)
            stroke(ctx)
        else
            shade = mplot.mesh.ctx.ndim==3 ? mplot.shades[i] : 1.0
            draw_surface_cell!(ctx, mplot, elem, shade)
        end

        # point labels
        Cairo.save(ctx)
        set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
        if mplot.node_labels
            for node in elem.nodes
                x, y = data2user(mplot.canvas, node.coord[1], node.coord[2])
                set_source_rgb(ctx, _colors_dict[:blue]...)
                draw_text(ctx, x, y, string(node.id), halign="center", valign="center", angle=0)
            end
        end
        Cairo.restore(ctx)
    end

    # draw colorbar
    mplot.field_kind != :none && draw!(mplot, ctx, mplot.colorbar)

    # draw axes
    if mplot.axes_loc != :none
        if mplot.axes_loc == :bottomleft
            x = mplot.outerpad
            y = mplot.canvas.box[4] - mplot.axes.height
        elseif mplot.axes_loc == :topleft
            x = mplot.outerpad
            y = mplot.canvas.box[2]
        else
            error("DomainPlot: axes location $(mplot.axes_loc) not implemented")
        end
        move_to(ctx, x, y)
        draw!(ctx, mplot.axes)
    end

end


function draw_surface_cell!(ctx::CairoContext, mplot::DomainPlot, elem::AbstractCell, shade::Float64)
    set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)

    function get_edge_color(face_color, edge_color)
        if mplot.view_mode == :surface # disguise edges in surface view
            face_color
        elseif edge_color == :auto
            face_color.*0.55
        else
            _colors_dict[mplot.edge_color]
        end
    end

    function get_outline_color(face_color, edge_color)
        if edge_color == :auto
            face_color.*0.55
        else
            _colors_dict[mplot.edge_color].*0.55
        end
    end

    # constant_color = !is_nodal_field || mplot.interpolation==:constant
    constant_color = mplot.field_kind in (:element,:none) || mplot.interpolation==:constant

    # compute linear gradients
    if constant_color
        if mplot.field_kind == :node
            val = sum( mplot.values[node.id] for node in elem.nodes)/length(elem.nodes)
            color = mplot.colormap(val).*shade
        elseif mplot.field_kind == :element
            color = mplot.colormap(mplot.values[elem.id]).*shade
        else
            color = mplot.face_color.*shade
        end

        edge_color    = get_edge_color(color, mplot.edge_color)
        outline_color = get_outline_color(color, mplot.edge_color)
    else
        cm = mplot.colormap
        npoints = elem.shape.base_shape.npoints # use only corner nodes
        nodes = elem.nodes[1:npoints]

        # regression plane
        values = [ mplot.values[node.id] for node in nodes ]
        X = [ elem.nodes[i].coord[j] for i in 1:npoints, j in 1:3]

        X[:,3] .= 1.0 # for regression
        A = pinv(X)*values # regression plane coefficients
        V = normalize(A[1:2]) # gradient direction

        # distances along V and sorting
        d = [ dot(V, node.coord[1:2]) for node in nodes ]
        idxs = sortperm(d)
        nodes = nodes[idxs]
        values = values[idxs]

        # compute gradient extremes
        Xa = nodes[1].coord[1:2]
        Xb = nodes[end].coord[1:2]
        ℓ = norm(Xb - Xa)

        vmin = dot(A, [Xa[1], Xa[2], 1.0])
        vmax = dot(A, [Xb[1], Xb[2], 1.0])

        # compute stops
        stops = Float64[]
        for i in 1:npoints
            val  = values[i]
            stop = clamp(round((val-vmin)/(vmax-vmin), digits=8), 0.0, 1.0) # clamp is required
            push!(stops, stop)
        end

        pat = pattern_create_linear(0, 0, V[1], V[2])
        mat = CairoMatrix(1/ℓ, 0, 0, 1/ℓ,  -Xa[1]/ℓ, -Xa[2]/ℓ)
        Cairo.set_matrix(pat, mat) # linear gradients need transformation
        for (val, stop) in zip(values, stops)
            color = cm(val).*shade
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        # edges
        factor = 0.65
        if mplot.view_mode == :surface # disguise edges in surface view
            factor = 1.0
            # mplot.edge_width = 0.1
        end

        edge_pat = pattern_create_linear(0, 0, V[1], V[2])
        mat = CairoMatrix(1/ℓ, 0, 0, 1/ℓ,  -Xa[1]/ℓ, -Xa[2]/ℓ)
        Cairo.set_matrix(edge_pat, mat)
        for (val, stop) in zip(values, stops)
            # color = cm(val).*(shade*factor)
            color = get_edge_color(cm(val).*shade, mplot.edge_color)
            pattern_add_color_stop_rgb(edge_pat, stop, color...)
        end

        # outline
        outline_pat = pattern_create_linear(0, 0, V[1], V[2])
        mat = CairoMatrix(1/ℓ, 0, 0, 1/ℓ,  -Xa[1]/ℓ, -Xa[2]/ℓ)
        Cairo.set_matrix(outline_pat, mat)
        for (val, stop) in zip(values, stops)
            # color = cm(val).*(shade*0.55)
            color = get_outline_color(cm(val).*shade, mplot.edge_color)
            pattern_add_color_stop_rgb(outline_pat, stop, color...)
        end
    end


    # draw cells face
    center       = sum(get_coords(elem), dims=1)[1:2]/length(elem.nodes)
    edges        = get_edges(elem)
    show_surface = mplot.view_mode in (:surface, :surface_with_edges)
    show_edges   = mplot.view_mode in (:surface, :surface_with_edges, :wireframe) # edges are disguised in surface view
    if mplot.view_mode != :outline

        if constant_color
            # draw element
            x, y = edges[1].nodes[1].coord
            new_path(ctx)
            move_to(ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            close_path(ctx)

            if show_surface
                set_source_rgb(ctx, color...)
                fill_preserve(ctx)
            end

            if show_edges
                set_source_rgb(ctx, edge_color...)
                set_line_width(ctx, mplot.edge_width)
                stroke(ctx)
            end
        elseif mplot.interpolation==:linear

            # draw element
            x, y = edges[1].nodes[1].coord
            new_path(ctx)
            move_to(ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                # mplot.view_mode == :surface && expand_points!(center, pts)
                curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            close_path(ctx)

            if show_surface
                set_source(ctx, pat)
                fill_preserve(ctx)
            end

            if show_edges
                set_source(ctx, edge_pat)
                set_line_width(ctx, mplot.edge_width)
                stroke(ctx)
            end
        else # nonlinear
            # set pattern mesh for nonlinear gradient
            pattern = CairoPatternMesh()
            mesh_pattern_begin_patch(pattern)

            x, y = edges[1].nodes[1].coord
            mesh_pattern_move_to(pattern, x, y)
            nedges = length(edges)

            if length(edges)==3 # for triangles
                for edge in edges[1:2]
                    x, y = edge.nodes[2].coord
                    mesh_pattern_line_to(pattern, x, y)
                end
            else
                for edge in edges
                    pts = bezier_points(edge)
                    mesh_pattern_curve_to(pattern, pts[2]..., pts[3]..., pts[4]...)
                end
            end

            # elem colors
            color = mplot.face_color
            if mplot.field_kind == :element
                color = mplot.colormap(mplot.values[elem.id])
            end

            # set nodal colors
            for (i,node) in enumerate(elem.nodes[1:nedges])
                # if has_field && is_nodal_field
                if mplot.field_kind == :node
                    color = mplot.colormap(mplot.values[node.id])
                end
                scolor = color.*shade # apply shade
                mesh_pattern_set_corner_color_rgb(pattern, i-1, scolor...)
            end
            mesh_pattern_end_patch(pattern)

            # draw element
            x, y = edges[1].nodes[1].coord
            new_path(ctx)
            move_to(ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            if show_surface
                set_source(ctx, pattern)
                fill_preserve(ctx)
                # paint(ctx) # only if using pattern instead of a path
            end

            if show_edges
                if mplot.view_mode == :surface
                    set_source(ctx, pattern)
                else
                    set_source(ctx, edge_pat)
                end
                set_line_width(ctx, mplot.edge_width)
                stroke(ctx)
            end
        end
    end

    # draw feature edges
    set_line_cap(ctx, Cairo.CAIRO_LINE_CAP_ROUND)

    if mplot.show_feature_edges
        new_path(ctx) # clear path, e.g. when last command used preserve
        set_line_width(ctx, mplot.outline_width)
        if constant_color
            set_source_rgb(ctx, outline_color...)
        else
            set_source(ctx, outline_pat)
        end

        for edge in edges
            node_idxs = sort([ node.id for node in edge.nodes ])
            haskey(mplot.feature_edges_d, node_idxs) || continue

            x, y = edge.nodes[1].coord
            move_to(ctx, x, y)
            pts = bezier_points(edge)
            curve_to(ctx, pts[2]..., pts[3]..., pts[4]...)
            set_line_width(ctx, mplot.outline_width)
            stroke(ctx)
        end
    end

end

