# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


"""
    DomainPlot(;
        size=(220,150),
        colorbar=:right, colorbar_ratio=0.9,
        font="NewComputerModern", font_size=7.0,
        title="",
        light_vector=[0,0,0],
        azimuth=30, elevation=30, distance=0.0, up=:z,
        axes=:none, axis_labels=String[],
        quiet=false)

    DomainPlot(mesh;
        size=(220,150), face_color=:aliceblue, warp=0.0,
        edge_width=0.3, edge_color=:auto, outline_width=0.4,
        line_elem_width=0.6,
        field="", limits=Float64[], field_kind=:auto, field_mult=1.0,
        label="", colormap=:coolwarm, diverging=false,
        colorbar=:right, colorbar_ratio=0.9, bins=6,
        font="NewComputerModern", font_size=7.0,
        title="",
        interpolation=:linear,
        azimuth=30, elevation=30, distance=0.0, up=:z,
        feature_edges=true, view_mode=:surface_with_edges,
        light_vector=[0,0,0],
        mark=:none, mark_size=2.5, stage=nothing,
        node_labels=false,
        axes=:none, axis_labels=String[],
        quiet=false)

Create a customizable domain plot for meshes and FE models.

`DomainPlot(mesh; ...)` is the compatibility constructor for single-domain plots.
For layered compositions, construct an empty `DomainPlot(; ...)` and append
layers with `add_plot`.
"""

mutable struct DomainPlotLayer
    mesh::AbstractDomain
    nodes::Vector{Node}
    elems::Vector{AbstractCell}
    feature_edges_d::Dict
    values::Vector{Float64}
    shades::Vector{Float64}
    field_kind::Symbol
    node_marker_kinds::Dict{Int,Symbol}

    edge_width::Float64
    edge_color::Union{Symbol,Tuple}
    outline_width::Float64
    line_elem_width::Float64
    face_color::Tuple
    field::String
    auto_limits::Bool
    limits::Vector{Float64}
    field_mult::Float64
    warp::Float64
    label::String
    colormap::Colormap
    diverging::Bool
    colorbar_location::Symbol
    colorbar_ratio::Float64
    colorbar::Colorbar
    bins::Int
    interpolation::Symbol
    show_feature_edges::Bool
    view_mode::Symbol
    mark::Symbol
    mark_size::Float64
    stage::Any
    node_labels::Bool
end


struct DomainRenderElem
    layer::DomainPlotLayer
    elem::AbstractCell
    shade::Float64
    layer_index::Int
    index_in_layer::Int
end


mutable struct DomainPlot <: Figure
    mesh::Union{AbstractDomain,Nothing}
    figure_frame::Frame
    title_box::TextBox
    canvas::Canvas
    colorbar::Colorbar
    colorbars::Vector{Colorbar}
    axes::AxesWidget
    nodes::Vector{Node}
    elems::Vector{AbstractCell}
    feature_edges_d::Dict
    values::Vector{Float64}
    outerpad::Float64
    shades::Vector{Float64}
    azimuth::Float64
    elevation::Float64
    distance::Float64
    up::Symbol

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
    auto_limits::Bool
    field_mult::Float64
    warp::Float64
    label::AbstractString
    colormap::Colormap
    diverging::Bool
    colorbar_location::Symbol
    colorbar_ratio::Float64
    bins::Int
    font::String
    font_size::Float64
    show_feature_edges::Bool
    view_mode::Symbol
    mark::Symbol
    mark_size::Float64
    stage::Any
    node_marker_kinds::Dict{Int,Symbol}
    node_labels::Bool
    axes_loc::Symbol
    axis_labels::Vector{AbstractString}
    left_items::Vector{FigureComponent}
    right_items::Vector{FigureComponent}
    top_items::Vector{FigureComponent}
    bottom_items::Vector{FigureComponent}
    layers::Vector{DomainPlotLayer}
    render_elems::Vector{DomainRenderElem}
    ndim::Int
    quiet::Bool

    function DomainPlot(;
        size::Tuple{<:Real,<:Real}=(220,150),
        colorbar::Symbol=:right,
        colorbar_ratio::Real=0.9,
        font::AbstractString="NewComputerModern",
        font_size::Real=7.0,
        title::AbstractString="",
        light_vector::Vector{<:Real}=[0.0,0.0,0.0],
        azimuth::Real=30.0,
        elevation::Real=30.0,
        distance::Real=0.0,
        up::Symbol=:z,
        axes::Symbol=:none,
        axis_labels::Vector{<:AbstractString}=String[],
        outerpad=0.0,
        quiet::Bool=false,
    )
        @check !isempty(font) "font must be a non-empty string"
        @check font_size > 0 "font_size must be positive"
        @check colorbar_ratio > 0 "colorbar_ratio must be positive"
        @check distance >= 0 "distance must be non-negative"
        @check up in (:x, :y, :z) "up must be one of :x, :y, :z. Got $up"

        width, height = float.(size)
        axes_loc = axes
        axis_labels = length(axis_labels) == 0 ? AbstractString["x", "y", "z"] : AbstractString[string(label) for label in axis_labels]
        axis_arrow_length = max(28.0, 4.0*font_size)
        title_box = TextBox(title)
        axes_widget = AxesWidget(
            location=axes_loc,
            labels=axis_labels,
            font=font,
            font_size=font_size,
            azimuth=azimuth,
            elevation=elevation,
            distance=distance,
            up=up,
            arrow_length=axis_arrow_length,
        )

        if !quiet
            printstyled("Domain plot\n", bold=true, color=:cyan)
            println("  size: $(width) x $(height) pt")
            title != "" && println("  title: $(title)")
        end

        return new(
            nothing,
            Frame(0.0, 0.0, width, height),
            title_box,
            Canvas(),
            Colorbar(location=:none),
            Colorbar[],
            axes_widget,
            Node[],
            AbstractCell[],
            Dict(),
            Float64[],
            float(outerpad),
            Float64[],
            float(azimuth),
            float(elevation),
            float(distance),
            up,
            :linear,
            Float64[float(v) for v in light_vector],
            width,
            height,
            0.3,
            :auto,
            0.3,
            0.6,
            _colors_dict[:aliceblue],
            "",
            :auto,
            [0.0, 0.0],
            true,
            1.0,
            0.0,
            "",
            Colormap(:coolwarm),
            false,
            colorbar,
            float(colorbar_ratio),
            6,
            string(font),
            float(font_size),
            true,
            :surface_with_edges,
            :none,
            2.5,
            nothing,
            Dict{Int,Symbol}(),
            false,
            axes_loc,
            axis_labels,
            FigureComponent[],
            FigureComponent[],
            FigureComponent[],
            FigureComponent[],
            DomainPlotLayer[],
            DomainRenderElem[],
            0,
            quiet,
        )
    end
end


function DomainPlot(mesh;
    size::Tuple{<:Real,<:Real}=(220,150),
    face_color::Symbol=:aliceblue,
    warp::Real=0.0,
    edge_width::Real=0.3,
    edge_color::Union{Symbol,Tuple}=:auto,
    outline_width::Real=0.3,
    line_elem_width::Real=0.6,
    field::AbstractString="",
    field_kind::Symbol=:auto,
    limits::AbstractVector{<:Real}=Float64[],
    field_mult::Real=1.0,
    label::AbstractString="",
    colormap::Union{Symbol,Colormap}=:coolwarm,
    diverging::Bool=false,
    colorbar::Symbol=:right,
    colorbar_ratio::Real=0.9,
    bins::Int=6,
    font::AbstractString="NewComputerModern",
    font_size::Real=7.0,
    title::AbstractString="",
    light_vector::Vector{<:Real}=[0.0,0.0,0.0],
    interpolation::Symbol=:linear,
    azimuth::Real=30.0,
    elevation::Real=30.0,
    distance::Real=0.0,
    up::Symbol=:z,
    feature_edges::Bool=true,
    view_mode::Symbol=:surface_with_edges,
    mark::Symbol=:none,
    mark_size::Real=2.5,
    stage=nothing,
    node_labels::Bool=false,
    axes::Symbol=:none,
    axis_labels::Vector{<:AbstractString}=String[],
    outerpad=0.0,
    quiet::Bool=false,
)
    mplot = DomainPlot(
        size=size,
        colorbar=colorbar,
        colorbar_ratio=colorbar_ratio,
        font=font,
        font_size=font_size,
        title=title,
        light_vector=light_vector,
        azimuth=azimuth,
        elevation=elevation,
        distance=distance,
        up=up,
        axes=axes,
        axis_labels=axis_labels,
        outerpad=outerpad,
        quiet=quiet,
    )

    add_plot(
        mplot,
        mesh;
        face_color=face_color,
        warp=warp,
        edge_width=edge_width,
        edge_color=edge_color,
        outline_width=outline_width,
        line_elem_width=line_elem_width,
        field=field,
        field_kind=field_kind,
        limits=limits,
        field_mult=field_mult,
        label=label,
        colormap=colormap,
        diverging=diverging,
        colorbar=colorbar,
        colorbar_ratio=colorbar_ratio,
        bins=bins,
        interpolation=interpolation,
        feature_edges=feature_edges,
        view_mode=view_mode,
        mark=mark,
        mark_size=mark_size,
        stage=stage,
        node_labels=node_labels,
    )

    return mplot
end


function _domain_plot_title_font_size(mplot::DomainPlot)
    return 1.2 * mplot.font_size
end


const _domain_mark_list = (:none, :auto, :support)
const _domain_translation_keys = (:ux, :uy, :uz)
const _domain_rotation_keys = (:rx, :ry, :rz)
const _domain_regular_marker_kind = :regular



function _domain_edge_key(edge::AbstractCell)
    return Tuple(sort(Int[node.id for node in edge.nodes]))
end


function _domain_layer_has_field(layer::DomainPlotLayer)
    return !isempty(layer.field) && layer.field_kind != :none
end


function _domain_layer_field_kind(layer::DomainPlotLayer)
    return _domain_layer_has_field(layer) ? layer.field_kind : :none
end


function _domain_sync_legacy_layer_fields!(mplot::DomainPlot, layer::DomainPlotLayer)
    mplot.mesh = layer.mesh
    mplot.feature_edges_d = layer.feature_edges_d
    mplot.values = layer.values
    mplot.shades = layer.shades

    mplot.edge_width = layer.edge_width
    mplot.edge_color = layer.edge_color
    mplot.outline_width = layer.outline_width
    mplot.line_elem_width = layer.line_elem_width
    mplot.face_color = layer.face_color
    mplot.field = layer.field
    mplot.field_kind = layer.field_kind
    mplot.limits = layer.limits
    mplot.auto_limits = layer.auto_limits
    mplot.field_mult = layer.field_mult
    mplot.warp = layer.warp
    mplot.label = layer.label
    mplot.colormap = layer.colormap
    mplot.diverging = layer.diverging
    mplot.bins = layer.bins
    mplot.interpolation = layer.interpolation
    mplot.colorbar = layer.colorbar
    mplot.show_feature_edges = layer.show_feature_edges
    mplot.view_mode = layer.view_mode
    mplot.mark = layer.mark
    mplot.mark_size = layer.mark_size
    mplot.stage = layer.stage
    mplot.node_marker_kinds = layer.node_marker_kinds
    mplot.node_labels = layer.node_labels
end


"""
    add_plot(plot::DomainPlot, mesh; kwargs...)

Append one mesh layer to a `DomainPlot`.

Use this with the composable constructor:

```julia
plot = DomainPlot(size=(8cm, 6cm))
add_plot(plot, mesh1; face_color=:aliceblue)
add_plot(plot, model2; warp=1.0, view_mode=:outline, edge_color=:red)
save(plot, "out.pdf")
```

# Arguments
- `plot::DomainPlot`: Target plot figure (mutated).
- `mesh::Union{AbstractDomain,FEModel}`: Source mesh or FE model for this layer.

# Layer keywords
- `face_color::Symbol`: surface color for area/surface cells.
- `warp::Real`: displacement scale factor applied from nodal field `"U"`.
- `edge_width::Real`: internal edge width for area/surface cells.
- `edge_color::Union{Symbol,Tuple}`: edge color; `:auto` darkens the face color.
- `outline_width::Real`: boundary outline width.
- `line_elem_width::Real`: stroke width for line elements.
- `field::AbstractString`: scalar field name for coloring; empty disables coloring.
- `field_kind::Symbol`: `:auto | :none | :node | :element`.
- `limits::Vector{<:Real}`: field range `[min,max]`; empty vector enables auto range.
- `field_mult::Real`: multiplier applied to field values.
- `label::AbstractString`: colorbar label for the field-bearing layer.
- `colormap::Union{Symbol,Colormap}`: layer colormap.
- `diverging::Bool`: center the colormap at zero.
- `colorbar::Symbol`: `:none | :left | :right | :top | :bottom`; defaults to the figure setting.
- `colorbar_ratio::Real`: colorbar length factor for this layer; defaults to the figure setting.
- `bins::Int`: number of colorbar bins.
- `interpolation::Symbol`: `:constant | :linear | :nonlinear`.
- `feature_edges::Bool`: enhance feature lines (`:outline` forces on).
- `view_mode::Symbol`: `:surface_with_edges | :surface | :wireframe | :outline`.
- `mark::Symbol`: node marker mode. Use `:none`, `:auto`, or `:support`.
- `mark_size::Real`: node marker size in points.
- `stage`: optional analysis stage whose boundary conditions feed semantic markers.
- `node_labels::Bool`: show node ids for this layer.

# Notes
- All layers in one `DomainPlot` must have the same dimension.
- Shared figure options such as camera, axes, title, font, and `light_vector`
  belong to `DomainPlot`.
- Figure-level `colorbar` and `colorbar_ratio` act as defaults for layers unless
  overridden in `add_plot`.

# Returns
- The appended `DomainPlotLayer`.
"""
function add_plot(
    mplot::DomainPlot,
    mesh;
    face_color::Symbol=:aliceblue,
    warp::Real=0.0,
    edge_width::Real=0.3,
    edge_color::Union{Symbol,Tuple}=:auto,
    outline_width::Real=0.3,
    line_elem_width::Real=0.6,
    field::AbstractString="",
    field_kind::Symbol=:auto,
    limits::AbstractVector{<:Real}=Float64[],
    field_mult::Real=1.0,
    label::AbstractString="",
    colormap::Union{Symbol,Colormap}=:coolwarm,
    diverging::Bool=false,
    colorbar::Union{Symbol,Nothing}=nothing,
    colorbar_ratio::Union{Real,Nothing}=nothing,
    bins::Int=6,
    interpolation::Symbol=:linear,
    feature_edges::Bool=true,
    view_mode::Symbol=:surface_with_edges,
    mark::Symbol=:none,
    mark_size::Real=2.5,
    stage=nothing,
    node_labels::Bool=false,
)
    @check face_color in keys(_colors_dict) "face_color must be one of: $(collect(keys(_colors_dict))). Got $face_color"
    @check edge_color == :auto || edge_color in keys(_colors_dict) "edge_color must be :auto or one of: $(collect(keys(_colors_dict))). Got $edge_color"
    @check length(limits) in (0, 2) "limits must have length 2 (manual) or be empty (auto)"
    @check mark in _domain_mark_list "mark must be one of $(_domain_mark_list). Got $mark"
    @check mark_size > 0 "mark_size must be positive"
    @check mark != :support || stage !== nothing "DomainPlot: mark=:support requires a stage"
    @check stage === nothing || hasproperty(stage, :bcs) "DomainPlot: stage must provide a bcs field"
    @check view_mode in (:surface_with_edges, :surface, :wireframe, :outline) "DomainPlot: invalid view_mode $(repr(view_mode))"
    @check interpolation in (:constant, :linear, :nonlinear) "DomainPlot: invalid interpolation $(repr(interpolation))"
    colorbar === nothing || @check colorbar in (:none, :left, :right, :top, :bottom) "DomainPlot: invalid colorbar location $(repr(colorbar))"
    colorbar_ratio === nothing || @check colorbar_ratio > 0 "DomainPlot: colorbar_ratio must be positive"

    if !isempty(mplot.layers)
        @check mesh.ctx.ndim == mplot.layers[1].mesh.ctx.ndim "DomainPlot: all layers must have the same dimension"
    end
    layer_colorbar = something(colorbar, mplot.colorbar_location)
    layer_colorbar_ratio = float(something(colorbar_ratio, mplot.colorbar_ratio))

    layer = DomainPlotLayer(
        mesh,
        Node[],
        AbstractCell[],
        Dict(),
        Float64[],
        Float64[],
        field_kind,
        Dict{Int,Symbol}(),
        float(edge_width),
        edge_color,
        float(outline_width),
        float(line_elem_width),
        _colors_dict[face_color],
        string(field),
        length(limits) == 0,
        length(limits) == 0 ? [0.0, 0.0] : collect(float.(limits)),
        float(field_mult),
        float(warp),
        string(label),
        colormap isa Symbol ? Colormap(colormap) : colormap,
        diverging,
        layer_colorbar,
        layer_colorbar_ratio,
        Colorbar(location=:none),
        bins,
        interpolation,
        (feature_edges || view_mode == :outline) && view_mode != :wireframe,
        view_mode,
        mark,
        float(mark_size),
        stage,
        node_labels,
    )

    push!(mplot.layers, layer)
    _domain_sync_legacy_layer_fields!(mplot, layer)
    return layer
end


function _domain_bc_target_nodes(bc)
    if bc.kind == :node
        return bc.target
    elseif bc.kind in (:face, :edge)
        return [node for facet in bc.target for node in facet.nodes]
    else
        return Node[]
    end
end


function _domain_is_zero_bc(cond, node::Node)
    cond isa Real || return false
    x, y, z = node.coord
    return isapprox(evaluate(cond, x=x, y=y, z=z, t=0.0), 0.0; atol=1e-12)
end


function _domain_node_bc_classes(stage, ndim::Int)
    class_data = Dict{Int,NamedTuple}()
    stage === nothing && return class_data

    translations = Set(_domain_translation_keys[1:ndim])
    rotations = Set(_domain_rotation_keys)

    for bc in stage.bcs
        bc.kind == :body && continue
        for node in _domain_bc_target_nodes(bc)
            data = get(class_data, node.id, (
                natural = false,
                essential = false,
                zero_translations = Set{Symbol}(),
                zero_rotation = false,
            ))

            natural = data.natural
            essential = data.essential
            zero_translations = copy(data.zero_translations)
            zero_rotation = data.zero_rotation

            for (key, cond) in bc.conds
                dof = get_dof(node, key)
                if dof !== nothing && dof.name == key
                    essential = true
                    if key in translations && _domain_is_zero_bc(cond, node)
                        push!(zero_translations, key)
                    elseif key in rotations && _domain_is_zero_bc(cond, node)
                        zero_rotation = true
                    end
                elseif dof !== nothing && dof.natname == key
                    natural = true
                elseif bc.kind in (:face, :edge)
                    natural = true
                end
            end

            class_data[node.id] = (
                natural = natural,
                essential = essential,
                zero_translations = zero_translations,
                zero_rotation = zero_rotation,
            )
        end
    end

    return class_data
end


function _domain_resolve_node_marker_kinds(stage, ndim::Int)
    marker_kinds = Dict{Int,Symbol}()
    translations = Set(_domain_translation_keys[1:ndim])

    for (node_id, data) in _domain_node_bc_classes(stage, ndim)
        if data.zero_rotation
            marker_kinds[node_id] = :clamp
        elseif data.zero_translations == translations
            marker_kinds[node_id] = :full_translation
        elseif !isempty(data.zero_translations)
            marker_kinds[node_id] = :zero_translation
        elseif data.essential
            marker_kinds[node_id] = :essential
        elseif data.natural
            marker_kinds[node_id] = :natural
        end
    end

    return marker_kinds
end


function _domain_node_marker_style(mark_kind::Symbol)
    if mark_kind == :natural
        return (:circle, :lightgray, true)
    elseif mark_kind == :essential
        return (:triangle, :lightgray, true)
    elseif mark_kind == :zero_translation
        return (:triangle, :black, true)
    elseif mark_kind == :full_translation
        return (:square, :black, true)
    elseif mark_kind == :clamp
        return (:hash, :black, true)
    else
        return (:circle, :black, false)
    end
end


function _domain_draw_node_marker!(ctx::RenderContext, mplot::DomainPlot, layer::DomainPlotLayer, node::Node)
    layer.mark == :none && return

    mark_kind = get(layer.node_marker_kinds, node.id, _domain_regular_marker_kind)
    layer.mark == :support && mark_kind == _domain_regular_marker_kind && return

    mark, color, has_stroke = _domain_node_marker_style(mark_kind)

    cairo_ctx = ctx.cairo_ctx
    x, y = data2user(mplot.canvas, node.coord[1], node.coord[2])
    set_line_width(cairo_ctx, max(0.3, 0.1*layer.mark_size) * ctx.width_scale)
    draw_mark(cairo_ctx, x, y, mark, layer.mark_size, Color(color), Color(:black); draw_stroke=has_stroke)
end


function _domain_draw_node_markers!(ctx::RenderContext, mplot::DomainPlot, layer::DomainPlotLayer, nodes)
    layer.mark == :none && return
    cairo_ctx = ctx.cairo_ctx

    Cairo.save(cairo_ctx)
    reset_matrix!(ctx)
    for node in nodes
        _domain_draw_node_marker!(ctx, mplot, layer, node)
    end
    Cairo.restore(cairo_ctx)
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


function project_to_2d!(nodes, azimuth, elevation, distance, up::Symbol=:z)
    xmin, xmax = extrema(node.coord[1] for node in nodes)
    ymin, ymax = extrema(node.coord[2] for node in nodes)
    zmin, zmax = extrema(node.coord[3] for node in nodes)
    reflength = max(xmax-xmin, ymax-ymin, zmax-zmin)

    center = 0.5*Vec3(xmin+xmax, ymin+ymax, zmin+zmax)
    for node in nodes
        node.coord = project_view_point(node.coord, azimuth, elevation, distance, up=up, center=center, reflength=reflength)
    end

    xmin, xmax = extrema(node.coord[1] for node in nodes)
    ymin, ymax = extrema(node.coord[2] for node in nodes)

    l = max(xmax-xmin, ymax-ymin)
    for node in nodes
        x = (node.coord[1]-xmin)/l
        y = (node.coord[2]-ymin)/l
        node.coord = Vec3(x, y, node.coord[3])
    end
end


function _domain_configure_layer_field!(layer::DomainPlotLayer, mesh::AbstractDomain, nodes)
    has_field = _domain_layer_has_field(layer)
    if !has_field
        layer.values = Float64[]
        return false
    end

    field = layer.field
    node_fields = mesh.node_fields
    elem_fields = mesh.elem_fields

    layer.label == "" && (layer.label = field)

    fields = [string(name) for name in Iterators.flatten([keys(node_fields), keys(elem_fields)])]
    field in fields || error("DomainPlot: field $field not found. Available fields are: $(join(fields, ", ", " and "))")

    if layer.field_kind == :auto
        if haskey(node_fields, field)
            layer.field_kind = :node
        elseif haskey(elem_fields, field)
            layer.field_kind = :element
        else
            all_fields = string.([collect(keys(node_fields)); collect(keys(elem_fields))])
            error("DomainPlot: field $field not found. Available fields are: $(join(all_fields, ", ", " and "))")
        end
    end

    if layer.field_kind == :element
        fvals = collect(float.(elem_fields[field] .* layer.field_mult))
        fmax = maximum(fvals)
        fmin = minimum(fvals)
    elseif layer.field_kind == :node
        fvals = collect(float.(node_fields[field] .* layer.field_mult))
        fmax = maximum(fvals[node.id] for node in nodes)
        fmin = minimum(fvals[node.id] for node in nodes)
    else
        error("DomainPlot: invalid field_kind $(repr(layer.field_kind)) for field $(repr(field))")
    end

    if fmin == fmax
        fmin -= 1
        fmax += 1
    end

    layer.values = fvals
    if layer.auto_limits
        layer.limits = [fmin, fmax]
    else
        fmin, fmax = layer.limits
    end

    layer.colormap = resize(layer.colormap, fmin, fmax, diverging=layer.diverging)
    return true
end


function _domain_prepare_layer!(mplot::DomainPlot, layer::DomainPlotLayer)
    mesh = copy(layer.mesh)
    ndim = mesh.ctx.ndim

    if layer.warp > 0.0
        U = get(mesh.node_fields, "U", nothing)
        if U === nothing
            alert("DomainPlot: Vector field U not found for warping.")
        else
            for node in mesh.nodes
                node.coord = node.coord + layer.warp*U[node.id,:]
            end
        end
    end

    active_elems = select(mesh.elems, :active)
    if ndim == 2
        areacells = [elem for elem in active_elems if elem.role == :solid]
        linecells = [elem for elem in active_elems if elem.role == :line]
        feature_edges = get_outer_facets(areacells)

        layer.elems = AbstractCell[areacells; linecells]
        layer.nodes = get_nodes(layer.elems)
        layer.shades = ones(Float64, length(layer.elems))
    else
        volcells = [elem for elem in active_elems if elem.role == :solid && elem.shape.ndim == 3]
        areacells = [elem for elem in active_elems if elem.role == :solid && elem.shape.ndim == 2]
        linecells = [elem for elem in active_elems if elem.role == :line]
        surfcells = get_outer_facets(volcells)
        feature_edges = get_feature_edges(surfcells)

        for cell in surfcells
            cell.id = cell.owner.id
        end

        layer.elems = AbstractCell[surfcells; areacells; linecells]
        layer.nodes = unique(node -> node.id, [node for elem in layer.elems for node in elem.nodes])

        Vview = Vec3(cosd(mplot.elevation)*cosd(mplot.azimuth), cosd(mplot.elevation)*sind(mplot.azimuth), sind(mplot.elevation))
        V = view_to_world(Vview, up=mplot.up)
        L = norm(mplot.light_vector) == 0 ? V : Vec3(mplot.light_vector)
        shades = zeros(Float64, length(layer.elems))
        for (i, elem) in enumerate(layer.elems)
            elem.role == :surface || continue
            N = get_facet_normal(elem)
            dot(N, L) < 0 && (N = -N)
            R = normalize(2*N*dot(L, N) - L)
            dot(V, R) < 0 && (R = -R)

            Ia = 0.68
            Id = 0.34
            Is = 0.20
            sh = 2.00

            shades[i] = Ia + Id*max(0, dot(L, N)) + Is*max(0, dot(V, R))^sh
        end
        layer.shades = shades
    end

    empty!(layer.feature_edges_d)
    for edge in feature_edges
        layer.feature_edges_d[_domain_edge_key(edge)] = edge
    end

    _domain_configure_layer_field!(layer, mesh, layer.nodes)
    layer.node_marker_kinds = layer.mark in (:auto, :support) ? _domain_resolve_node_marker_kinds(layer.stage, ndim) : Dict{Int,Symbol}()
    return nothing
end


function _domain_render_depth_key(render::DomainRenderElem)
    zavg = sum(node.coord[3] for node in render.elem.nodes)/length(render.elem.nodes)
    zmin = minimum(node.coord[3] for node in render.elem.nodes)
    depth = 0.9*zavg + 0.1*zmin
    return (-depth, render.layer_index, render.index_in_layer)
end


function _domain_render_2d_key(render::DomainRenderElem)
    role_rank = render.elem.role == :line ? 1 : 0
    area_key = role_rank == 0 ? -cell_extent(render.elem) : 0.0
    return (render.layer_index, role_rank, area_key, render.index_in_layer)
end


function _domain_reset_render_state!(mplot::DomainPlot)
    mplot.nodes = Node[]
    mplot.elems = AbstractCell[]
    mplot.render_elems = DomainRenderElem[]
    mplot.feature_edges_d = Dict()
    mplot.values = Float64[]
    mplot.shades = Float64[]
end


function _domain_prepare_layers!(mplot::DomainPlot)
    field_layers = DomainPlotLayer[]

    for layer in mplot.layers
        _domain_prepare_layer!(mplot, layer)
        append!(mplot.nodes, layer.nodes)
        _domain_layer_has_field(layer) && push!(field_layers, layer)
    end

    return field_layers
end


function _domain_build_render_elems!(mplot::DomainPlot)
    for (layer_index, layer) in enumerate(mplot.layers)
        for (index_in_layer, elem) in enumerate(layer.elems)
            shade = mplot.ndim == 3 ? layer.shades[index_in_layer] : 1.0
            push!(mplot.render_elems, DomainRenderElem(layer, elem, shade, layer_index, index_in_layer))
        end
    end

    if mplot.ndim == 3
        sort!(mplot.render_elems, by=_domain_render_depth_key)
    else
        sort!(mplot.render_elems, by=_domain_render_2d_key)
    end

    mplot.elems = AbstractCell[render.elem for render in mplot.render_elems]
    mplot.shades = Float64[render.shade for render in mplot.render_elems]
    return nothing
end


function _domain_title_metrics(mplot::DomainPlot)
    isempty(mplot.title_box.text) && return (0.0, 0.0)

    surf = CairoImageSurface(4, 4, Cairo.FORMAT_ARGB32)
    cc = CairoContext(surf)
    title_font_size = _domain_plot_title_font_size(mplot)
    select_font_face(cc, get_font(mplot.font), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(cc, title_font_size)
    title_height = getsize(cc, mplot.title_box.text, title_font_size)[2]
    title_gap = 0.6 * mplot.font_size
    return (title_height, title_gap)
end


function _domain_reset_side_items!(mplot::DomainPlot)
    mplot.colorbars = Colorbar[]
    mplot.left_items = FigureComponent[]
    mplot.right_items = FigureComponent[]
    mplot.top_items = FigureComponent[]
    mplot.bottom_items = FigureComponent[]
    return nothing
end


function _domain_push_side_item!(mplot::DomainPlot, side::Symbol, item::FigureComponent)
    if side == :left
        push!(mplot.left_items, item)
    elseif side == :right
        push!(mplot.right_items, item)
    elseif side == :top
        push!(mplot.top_items, item)
    elseif side == :bottom
        push!(mplot.bottom_items, item)
    else
        error("DomainPlot: invalid side $(repr(side))")
    end
    return nothing
end


function _domain_configure_colorbars!(mplot::DomainPlot)
    _domain_reset_side_items!(mplot)

    for layer in mplot.layers
        if !_domain_layer_has_field(layer) || layer.colorbar_location == :none
            layer.colorbar = Colorbar(location=:none)
            continue
        end

        layer.colorbar = Colorbar(
            location=layer.colorbar_location,
            label=layer.label,
            length_factor=layer.colorbar_ratio,
            colormap=layer.colormap,
            limits=layer.limits,
            font_size=mplot.font_size,
            font=mplot.font,
            bins=layer.bins,
        )
        configure!(mplot, layer.colorbar)
        push!(mplot.colorbars, layer.colorbar)
        _domain_push_side_item!(mplot, layer.colorbar.location, layer.colorbar)
    end

    mplot.colorbar = isempty(mplot.colorbars) ? Colorbar(location=:none) : mplot.colorbars[1]
    return nothing
end


function _domain_side_pane_size(items::Vector{FigureComponent}, side::Symbol)
    isempty(items) && return 0.0
    if side in (:left, :right)
        return maximum(item.width for item in items)
    else
        return maximum(item.height for item in items)
    end
end


function _domain_assign_canvas_frame!(mplot::DomainPlot, lpane, rpane, tpane, bpane)
    width, height = mplot.width, mplot.height
    base_x = mplot.figure_frame.x
    base_y = mplot.figure_frame.y
    xmin = base_x + mplot.outerpad + lpane
    ymin = base_y + mplot.outerpad + tpane
    xmax = base_x + width - rpane - mplot.outerpad
    ymax = base_y + height - bpane - mplot.outerpad
    mplot.canvas.frame = Frame(xmin, ymin, xmax - xmin, ymax - ymin)
    return nothing
end


function _domain_assign_title_frame!(mplot::DomainPlot, title_height)
    if isempty(mplot.title_box.text)
        mplot.title_box.frame = Frame()
        return nothing
    end

    canvas = mplot.canvas.frame
    base_y = mplot.figure_frame.y
    mplot.title_box.frame = Frame(canvas.x, base_y + mplot.outerpad, canvas.width, title_height)
    mplot.title_box.angle = 0.0
    return nothing
end


function _domain_assign_canvas_limits!(mplot::DomainPlot)
    xmin_data, xmax_data = extrema(node.coord[1] for node in mplot.nodes)
    ymin_data, ymax_data = extrema(node.coord[2] for node in mplot.nodes)

    frame = mplot.canvas.frame
    ratio = min(frame.width/(xmax_data-xmin_data), frame.height/(ymax_data-ymin_data))
    dX = 0.5 * (frame.width - ratio*(xmax_data-xmin_data))
    dY = 0.5 * (frame.height - ratio*(ymax_data-ymin_data))

    xmin_data -= dX/ratio
    xmax_data += dX/ratio
    ymin_data -= dY/ratio
    ymax_data += dY/ratio

    mplot.canvas.limits = [xmin_data, ymin_data, xmax_data, ymax_data]
    return nothing
end


function _domain_configure_axes!(mplot::DomainPlot)
    if mplot.axes_loc == :none
        mplot.axes = AxesWidget(location=:none)
        return nothing
    end

    ndim = mplot.ndim
    canvas = mplot.canvas.frame
    mplot.axes = AxesWidget(
        location=mplot.axes_loc,
        font_size=mplot.font_size,
        font=mplot.font,
        azimuth=mplot.azimuth,
        elevation=mplot.elevation,
        distance=mplot.distance,
        up=mplot.up,
        arrow_length=max(28.0, 4.0*mplot.font_size),
        labels=mplot.axis_labels[1:ndim],
    )
    configure!(mplot.axes)

    if mplot.axes_loc == :bottom_left
        mplot.axes.frame = Frame(canvas.x, canvas.y + canvas.height - mplot.axes.height, mplot.axes.width, mplot.axes.height)
    elseif mplot.axes_loc == :top_left
        mplot.axes.frame = Frame(canvas.x, canvas.y, mplot.axes.width, mplot.axes.height)
    elseif mplot.axes_loc == :bottom_right
        mplot.axes.frame = Frame(canvas.x + canvas.width - mplot.axes.width, canvas.y + canvas.height - mplot.axes.height, mplot.axes.width, mplot.axes.height)
    elseif mplot.axes_loc == :top_right
        mplot.axes.frame = Frame(canvas.x + canvas.width - mplot.axes.width, canvas.y, mplot.axes.width, mplot.axes.height)
    else
        error("DomainPlot: axes location $(mplot.axes_loc) not implemented")
    end

    return nothing
end


function _domain_assign_side_frames!(mplot::DomainPlot, items::Vector{FigureComponent}, side::Symbol, pane_size::Float64)
    isempty(items) && return nothing

    canvas = mplot.canvas.frame
    if side in (:left, :right)
        slot_length = canvas.height / length(items)
        pane_x = side == :left ? canvas.x - pane_size : canvas.x + canvas.width
        for (i, item) in enumerate(items)
            cb = item::Colorbar
            slot_y = canvas.y + (i - 1) * slot_length
            cb.height = cb.length_factor * slot_length
            cb.axis.height = cb.height
            cb.frame = Frame(pane_x, slot_y, pane_size, slot_length)
        end
    else
        slot_length = canvas.width / length(items)
        pane_y = side == :top ? canvas.y - pane_size : canvas.y + canvas.height
        for (i, item) in enumerate(items)
            cb = item::Colorbar
            slot_x = canvas.x + (i - 1) * slot_length
            cb.width = cb.length_factor * slot_length
            cb.axis.width = cb.width
            cb.frame = Frame(slot_x, pane_y, slot_length, pane_size)
        end
    end

    return nothing
end


function configure!(mplot::DomainPlot)
    length(mplot.layers) > 0 || throw(SerendipException("DomainPlot: no plot layers added"))

    mplot.ndim = mplot.layers[1].mesh.ctx.ndim
    _domain_reset_render_state!(mplot)
    field_layers = _domain_prepare_layers!(mplot)

    @check !isempty(mplot.nodes) "DomainPlot: no drawable entities found"

    if mplot.ndim == 3
        project_to_2d!(mplot.nodes, mplot.azimuth, mplot.elevation, mplot.distance, mplot.up)
    end

    _domain_build_render_elems!(mplot)
    active_layer = isempty(field_layers) ? last(mplot.layers) : field_layers[1]
    _domain_sync_legacy_layer_fields!(mplot, active_layer)

    width, height = mplot.width, mplot.height
    mplot.figure_frame = Frame(mplot.figure_frame.x, mplot.figure_frame.y, width, height)
    mplot.outerpad = max(0.01*min(width, height), mplot.outerpad)
    title_height, title_gap = _domain_title_metrics(mplot)
    _domain_configure_colorbars!(mplot)

    lpane = _domain_side_pane_size(mplot.left_items, :left)
    rpane = _domain_side_pane_size(mplot.right_items, :right)
    top_colorbar_pane = _domain_side_pane_size(mplot.top_items, :top)
    bottom_colorbar_pane = _domain_side_pane_size(mplot.bottom_items, :bottom)
    tpane = title_height + title_gap + top_colorbar_pane

    _domain_assign_canvas_frame!(mplot, lpane, rpane, tpane, bottom_colorbar_pane)
    _domain_assign_title_frame!(mplot, title_height)
    _domain_assign_canvas_limits!(mplot)
    _domain_configure_axes!(mplot)
    _domain_assign_side_frames!(mplot, mplot.left_items, :left, lpane)
    _domain_assign_side_frames!(mplot, mplot.right_items, :right, rpane)
    _domain_assign_side_frames!(mplot, mplot.top_items, :top, top_colorbar_pane)
    _domain_assign_side_frames!(mplot, mplot.bottom_items, :bottom, bottom_colorbar_pane)
end


function iscounterclockwise(points::Vector{Node})
    val = 0.0
    n = length(points)
    for i in 1:n
        x1, y1 = points[i].coord
        x2, y2 = points[mod1(i+1, n)].coord
        val += (x1*y2) - (x2*y1)
    end
    return val > 0
end


function draw_contents!(mplot::DomainPlot, ctx::RenderContext)
    cairo_ctx = ctx.cairo_ctx

    set_line_join(cairo_ctx, Cairo.CAIRO_LINE_JOIN_ROUND)
    set_line_cap(cairo_ctx, Cairo.CAIRO_LINE_CAP_ROUND)
    font = get_font(mplot.font)
    select_font_face(cairo_ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(cairo_ctx, mplot.font_size)

    Xmin = mplot.canvas.frame.x
    Ymin = mplot.canvas.frame.y
    Xmax = mplot.canvas.frame.x + mplot.canvas.frame.width
    Ymax = mplot.canvas.frame.y + mplot.canvas.frame.height
    xmin, ymin, xmax, ymax = mplot.canvas.limits

    ratio = (Xmax-Xmin)/(xmax-xmin)
    set_local_matrix!(ctx, CairoMatrix([ratio, 0, 0, -ratio, Xmin-xmin*ratio, Ymax+ymin*ratio]...))

    for render in mplot.render_elems
        layer = render.layer
        elem = render.elem

        elem.role in (:line, :solid, :surface) || continue

        if elem.owner !== nothing && layer.view_mode in (:surface, :surface_with_edges)
            if elem.owner.shape.ndim == 3
                nnodes_basic = elem.shape.base_shape.npoints
                !iscounterclockwise(elem.nodes[1:nnodes_basic]) && continue
            end
        end

        if elem.role == :line
            x, y = elem.nodes[1].coord
            move_to(cairo_ctx, x, y)
            color = Vec3(0.8, 0.2, 0.1)
            if _domain_layer_field_kind(layer) == :element
                color = layer.colormap(layer.values[elem.id])
            end

            pts = bezier_points(elem)
            curve_to(cairo_ctx, pts[2]..., pts[3]..., pts[4]...)
            set_line_width(cairo_ctx, layer.line_elem_width * ctx.width_scale)
            set_source_rgb(cairo_ctx, color...)
            stroke(cairo_ctx)
        else
            draw_surface_cell!(ctx, layer, elem, render.shade)
        end

        if mplot.ndim == 3
            _domain_draw_node_markers!(ctx, mplot, layer, elem.nodes)
        end

        Cairo.save(cairo_ctx)
        reset_matrix!(ctx)
        if layer.node_labels
            for node in elem.nodes
                x, y = data2user(mplot.canvas, node.coord[1], node.coord[2])
                y += mplot.font_size
                set_source_rgb(cairo_ctx, _colors_dict[:blue]...)
                draw_text(cairo_ctx, x, y, string(node.id), halign="center", valign="center", angle=0)
            end
        end
        Cairo.restore(cairo_ctx)
        set_local_matrix!(ctx, CairoMatrix([ratio, 0, 0, -ratio, Xmin-xmin*ratio, Ymax+ymin*ratio]...))
    end

    if mplot.ndim != 3
        for layer in mplot.layers
            _domain_draw_node_markers!(ctx, mplot, layer, layer.nodes)
        end
    end

    for cb in mplot.colorbars
        cb.location != :none && draw!(mplot, ctx, cb)
    end

    if mplot.axes_loc != :none
        draw!(mplot.axes, ctx)
    end

    if !isempty(mplot.title_box.text)
        reset_matrix!(ctx)
        set_source_rgb(cairo_ctx, 0.0, 0.0, 0.0)
        set_font_size(cairo_ctx, _domain_plot_title_font_size(mplot))
        _draw_text_box!(ctx, mplot.title_box)
        set_font_size(cairo_ctx, mplot.font_size)
    end
end


function draw_background!(mplot::DomainPlot, ctx::RenderContext)
    _draw_figure_background!(ctx, mplot.figure_frame, ctx.background)
end


function draw_surface_cell!(ctx::RenderContext, layer::DomainPlotLayer, elem::AbstractCell, shade::Float64)
    cairo_ctx = ctx.cairo_ctx
    set_line_cap(cairo_ctx, Cairo.CAIRO_LINE_CAP_ROUND)
    field_kind = _domain_layer_field_kind(layer)

    function get_edge_color(face_color, edge_color)
        if layer.view_mode == :surface
            face_color
        elseif edge_color == :auto
            face_color.*0.55
        else
            _colors_dict[layer.edge_color]
        end
    end

    function get_outline_color(face_color, edge_color)
        if edge_color == :auto
            face_color.*0.7
        else
            _colors_dict[layer.edge_color].*0.7
        end
    end

    constant_color = field_kind in (:element, :none) || layer.interpolation == :constant

    if constant_color
        if field_kind == :node
            val = sum(layer.values[node.id] for node in elem.nodes)/length(elem.nodes)
            color = layer.colormap(val).*shade
        elseif field_kind == :element
            color = layer.colormap(layer.values[elem.id]).*shade
        else
            color = layer.face_color.*shade
        end

        edge_color = get_edge_color(color, layer.edge_color)
        outline_color = get_outline_color(color, layer.edge_color)
    else
        cm = layer.colormap
        npoints = elem.shape.base_shape.npoints
        nodes = elem.nodes[1:npoints]

        values = [layer.values[node.id] for node in nodes]
        X = [elem.nodes[i].coord[j] for i in 1:npoints, j in 1:3]

        X[:,3] .= 1.0
        A = pinv(X)*values
        V = normalize(A[1:2])

        d = [dot(V, node.coord[1:2]) for node in nodes]
        idxs = sortperm(d)
        nodes = nodes[idxs]
        values = values[idxs]

        Xa = nodes[1].coord[1:2]
        Xb = nodes[end].coord[1:2]
        ℓ = norm(Xb - Xa)

        vmin = dot(A, [Xa[1], Xa[2], 1.0])
        vmax = dot(A, [Xb[1], Xb[2], 1.0])

        stops = Float64[]
        for i in 1:npoints
            val = values[i]
            stop = clamp(round((val-vmin)/(vmax-vmin), digits=8), 0.0, 1.0)
            push!(stops, stop)
        end

        pat = pattern_create_linear(0, 0, V[1], V[2])
        mat = CairoMatrix(1/ℓ, 0, 0, 1/ℓ, -Xa[1]/ℓ, -Xa[2]/ℓ)
        Cairo.set_matrix(pat, mat)
        for (val, stop) in zip(values, stops)
            color = cm(val).*shade
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        edge_pat = pattern_create_linear(0, 0, V[1], V[2])
        Cairo.set_matrix(edge_pat, mat)
        for (val, stop) in zip(values, stops)
            color = get_edge_color(cm(val).*shade, layer.edge_color)
            pattern_add_color_stop_rgb(edge_pat, stop, color...)
        end

        outline_pat = pattern_create_linear(0, 0, V[1], V[2])
        Cairo.set_matrix(outline_pat, mat)
        for (val, stop) in zip(values, stops)
            color = get_outline_color(cm(val).*shade, layer.edge_color)
            pattern_add_color_stop_rgb(outline_pat, stop, color...)
        end
    end

    edges = get_edges(elem)
    show_surface = layer.view_mode in (:surface, :surface_with_edges)
    show_edges = layer.view_mode in (:surface, :surface_with_edges, :wireframe)
    if layer.view_mode != :outline
        if constant_color
            x, y = edges[1].nodes[1].coord
            new_path(cairo_ctx)
            move_to(cairo_ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(cairo_ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            close_path(cairo_ctx)

            if show_surface
                set_source_rgb(cairo_ctx, color...)
                fill_preserve(cairo_ctx)
            end

            if show_edges
                set_source_rgb(cairo_ctx, edge_color...)
                set_line_width(cairo_ctx, layer.edge_width * ctx.width_scale)
                stroke(cairo_ctx)
            end
        elseif layer.interpolation == :linear
            x, y = edges[1].nodes[1].coord
            new_path(cairo_ctx)
            move_to(cairo_ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(cairo_ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            close_path(cairo_ctx)

            if show_surface
                set_source(cairo_ctx, pat)
                fill_preserve(cairo_ctx)
            end

            if show_edges
                set_source(cairo_ctx, edge_pat)
                set_line_width(cairo_ctx, layer.edge_width * ctx.width_scale)
                stroke(cairo_ctx)
            end
        else
            pattern = CairoPatternMesh()
            mesh_pattern_begin_patch(pattern)

            x, y = edges[1].nodes[1].coord
            mesh_pattern_move_to(pattern, x, y)
            nedges = length(edges)

            if length(edges) == 3
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

            color = layer.face_color
            if field_kind == :element
                color = layer.colormap(layer.values[elem.id])
            end

            for (i, node) in enumerate(elem.nodes[1:nedges])
                if field_kind == :node
                    color = layer.colormap(layer.values[node.id])
                end
                scolor = color.*shade
                mesh_pattern_set_corner_color_rgb(pattern, i-1, scolor...)
            end
            mesh_pattern_end_patch(pattern)

            x, y = edges[1].nodes[1].coord
            new_path(cairo_ctx)
            move_to(cairo_ctx, x, y)
            for edge in edges
                pts = bezier_points(edge)
                curve_to(cairo_ctx, pts[2]..., pts[3]..., pts[4]...)
            end

            if show_surface
                set_source(cairo_ctx, pattern)
                fill_preserve(cairo_ctx)
            end

            if show_edges
                if layer.view_mode == :surface
                    set_source(cairo_ctx, pattern)
                else
                    set_source(cairo_ctx, edge_pat)
                end
                set_line_width(cairo_ctx, layer.edge_width * ctx.width_scale)
                stroke(cairo_ctx)
            end
        end
    end

    set_line_cap(cairo_ctx, Cairo.CAIRO_LINE_CAP_ROUND)

    if layer.show_feature_edges
        new_path(cairo_ctx)
        set_line_width(cairo_ctx, layer.outline_width * ctx.width_scale)
        if constant_color
            set_source_rgb(cairo_ctx, outline_color...)
        else
            set_source(cairo_ctx, outline_pat)
        end

        for edge in edges
            haskey(layer.feature_edges_d, _domain_edge_key(edge)) || continue

            x, y = edge.nodes[1].coord
            move_to(cairo_ctx, x, y)
            pts = bezier_points(edge)
            curve_to(cairo_ctx, pts[2]..., pts[3]..., pts[4]...)
            set_line_width(cairo_ctx, layer.outline_width * ctx.width_scale)
            stroke(cairo_ctx)
        end
    end
end
