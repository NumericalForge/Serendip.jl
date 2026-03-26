# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

"""
    Chart(; 
        size=(220,150), font="NewComputerModern", font_size=7.0,
        xlimits, ylimits, aspect_ratio=:auto,
        xmult=1.0, ymult=1.0, xbins=7, ybins=6,
        title="", background=nothing,
        xlabel="\$x\$", ylabel="\$y\$",
        xticks=Float64[], yticks=Float64[],
        xtick_labels=String[], ytick_labels=String[],
        legend=:top_right, legend_font_size=0, legend_background=nothing,
        quiet=false)

Construct a 2D chart figure with axes, legend, and optional tick customization.

# Arguments
- `size::Tuple{<:Real,<:Real}`: width × height in points. Use `cm` as a convenience helper, e.g. `size=(8cm, 6cm)`.
- `font::AbstractString`: font family for axes and legend.
- `font_size::Real`: base font size.
- `xlimits::Vector{<:Real}`, `ylimits::Vector{<:Real}`: axis limits `[min,max]`; use empty vectors for auto scaling.
- `aspect_ratio::Symbol`: `:auto` or `:equal`.
- `xmult::Real`, `ymult::Real`: multiplicative factors applied to tick values.
- `xbins::Int`, `ybins::Int`: target number of major ticks.
- `title::AbstractString`: chart title, centered above the plot area.
- `background`: full-figure background fill; `nothing` leaves the figure unfilled.
- `xlabel::AbstractString`, `ylabel::AbstractString`: axis labels.
- `xticks::Vector{<:Real}`, `yticks::Vector{<:Real}`: explicit tick positions; empty vectors enable auto ticks.
- `xtick_labels::Vector{<:AbstractString}`, `ytick_labels::Vector{<:AbstractString}`: custom tick labels; if provided, lengths must match the corresponding tick arrays.
- `legend::Symbol`: legend location (e.g., `:top_right`, `:top_left`, `:bottom_left`, `:outer_right`).
- `legend_font_size::Real`: legend font size; `0` uses `font_size`.
- `legend_background`: legend box fill color; if unset it defaults to white standalone and follows the grid background when drawn inside a `ChartGrid`.
- `quiet::Bool`: suppress constructor log.

# Notes
- Use `add_series` to append data series to the chart.
- Use `add_annotation` to add plot-relative overlay annotations.
- Inline math in plot text may be delimited with either `\$...\$` or `` `...` ``.
- The legend is drawn after annotations.
- `background=nothing` leaves the chart background unfilled in vector outputs and uses a white page background in PNG.
- Use `save` to export the chart to a file.

# Returns
- A `Chart` object.

# Example
```julia
using Serendip: Chart, cm

ch = Chart(size=(8cm, 6cm),
           title="Response History",
           xlabel="Time [s]",
           ylabel="Displacement [mm]",
           xlimits=[0.0,10.0],
           ylimits=[-5.0,5.0],
           legend=:bottom_right)
```
"""
mutable struct Chart <: Figure
    width::Float64
    height::Float64
    figure_frame::Frame
    background::Union{Nothing,Color}
    title_box::TextBox
    canvas::Canvas
    xaxis::Axis
    yaxis::Axis
    dataseries::Vector{DataSeries}
    legend::Legend
    annotations::AbstractArray

    aspect_ratio::Symbol
    outerpad::Float64
    left_items::Vector{FigureComponent}
    right_items::Vector{FigureComponent}
    top_items::Vector{FigureComponent}
    bottom_items::Vector{FigureComponent}
    overlay_items::Vector{FigureComponent}
    icolor::Int
    iorder::Int

    quiet::Bool

    function Chart(;
        size=(220, 150),
        font="NewComputerModern",
        font_size::Real=7.0,
        xlimits=Float64[],
        ylimits=Float64[],
        aspect_ratio=:auto,
        xmult::Real=1.0,
        ymult::Real=1.0,
        xbins::Int=7,
        ybins::Int=6,
        title::AbstractString="",
        background=nothing,
        xlabel::AbstractString="\$x\$",
        ylabel::AbstractString="\$y\$",
        xticks::Vector{<:Real}=Float64[],
        yticks::Vector{<:Real}=Float64[],
        xtick_labels::Vector{<:AbstractString}=String[],
        ytick_labels::Vector{<:AbstractString}=String[],
        legend::Symbol=:top_right,
        legend_font_size::Real=0,
        legend_background=nothing,
        # colorbar=:right,
        # colorbar_scale=0.9,
        # colorbar_label="",
        # colorbar_limits=[0.0, 0.0],
        # colorbar_font_size=7.0,
        quiet=false
    )
        if legend_font_size == 0
            legend_font_size = font_size
        end

        @check font_size > 0
        @check legend_font_size > 0
        @check aspect_ratio in (:auto, :equal) "Chart: Invalid aspect_ratio: $aspect_ratio. Use :auto or :equal. Got $(repr(aspect_ratio))."

        width, height = size
        outerpad = 0.01 * min(width, height)
        the_legend = Legend(; location=legend, font=font, font_size=legend_font_size, background=legend_background, ncols=1)
        xaxis = Axis(direction=:horizontal, limits=xlimits, label=xlabel, font=font, font_size=font_size, ticks=xticks, tick_labels=xtick_labels, mult=xmult, bins=xbins)
        yaxis = Axis(direction=:vertical, limits=ylimits, label=ylabel, font=font, font_size=font_size, ticks=yticks, tick_labels=ytick_labels, mult=ymult, bins=ybins)
        background = resolve_color(background)
        title_box = TextBox(title)

        this = new(width, height, Frame(0.0, 0.0, width, height), background, title_box, Canvas(), xaxis, yaxis, [], the_legend, [],
            aspect_ratio, outerpad, FigureComponent[], FigureComponent[], FigureComponent[], FigureComponent[], FigureComponent[], 1, 1, quiet)

        if !quiet
            printstyled("Chart figure\n", bold=true, color=:cyan)
            println("  size: $(this.width) x $(this.height) pt")
            title != "" && println("  title: $(title)")
            background !== nothing && println("  background: $(repr(background))")
            println("  axes: $(xlabel) vs $(ylabel)")
        end

        return this
    end
end


_figure_background(c::Chart) = c.background

function _chart_gridline_color(c::Chart, ctx::RenderContext)
    background = ctx.background === nothing ? c.background : ctx.background
    return background === nothing ? Color(0.9, 0.9, 0.9) : darken(grayscale(background), 0.15)
end

"""
    add_series(chart::Chart, kind::Symbol, X::AbstractArray, Y::AbstractArray; kwargs...)
    add_series(chart::Chart, X::AbstractArray, Y::AbstractArray; kwargs...)

Append a `DataSeries` to `chart`.
The second version uses `kind = :line`.

# Arguments
- `chart::Chart` : Target chart (mutated).
- `kind::Symbol` : Plot type: `:line`, `:scatter`, `:bar`.
- `X, Y::AbstractArray` : Data vectors.

# Keyword options
- `line_style::Symbol = :solid` : Line style (e.g. `:solid`, `:dash`, ...).
- `dash::Vector{Float64} = Float64[]` : Custom dash pattern. If nonempty, overrides `line_style`.
- `color = :default` : Line/marker color. `:default` selects from the chart palette cyclically.
- `line_width::Float64 = 0.5` : Line width (> 0).
- `mark::Symbol = :none` : Mark shape.
- `mark_size::Float64 = 2.5` : Mark size (> 0).
- `mark_color = :white` : Mark fill color.
- `mark_stroke_color = :default` : Mark edge color (`:default` follows `color`).
- `label::AbstractString = ""` : Legend label.
- `tag::AbstractString = ""` : On-curve annotation text.
- `tag_location::Symbol = :top` : Relative location of tag (`:top`, `:bottom`, `:left`, `:right`).
- `tag_position::Float64 = 0.5` : Position along the curve in [0,1].
- `tag_alignment::Symbol = :horizontal` : Tag orientation (`:horizontal`, `:vertical`, `:parallel`).
- `bar_width::Float64 = 0.0` : Bar width in x-data units (`0` enables auto width).
- `bar_base::Float64 = 0.0` : Bar baseline in y-data units.
- `order::Int = 0` : Z-order. If `0`, an incremental order is assigned.

# Returns
- The series object.

# Examples
```julia
ch = Chart(size=(300,200), xlabel="Time [s]", ylabel="Displacement [mm]",
           xlimits=[0.0,10.0], ylimits=[-5.0,5.0], legend=:bottom_right)

add_line(ch, 0:0.1:10, sin.(0:0.1:10); label="sin")
```
"""
function add_series(chart::Chart, kind::Symbol, X::AbstractArray, Y::AbstractArray;
    line_style=:solid, dash=Float64[],
    color=:default,
    line_width=0.5,
    mark=:none, mark_size=2.5,
    mark_color=:white, mark_stroke_color=:default,
    label="", tag="", tag_location=:top, tag_position=0.5,
    tag_alignment=:horizontal,
    bar_width=0.0,
    bar_base=0.0,
    order=0
)

    @check line_width > 0 "Line width must be positive"
    @check mark_size > 0 "Mark size must be positive"
    @check tag_position >= 0 && tag_position <= 1 "Tag position must be in [0,1]"
    @check order >= 0 "Order must be non-negative"
    @check bar_width >= 0 "Bar width must be non-negative"
    @check kind in (:line, :scatter, :bar) "Invalid series kind: $kind. Use :line, :scatter, or :bar"
    @check length(X) == length(Y) "X and Y must have the same length"
    @check mark in _mark_list "Invalid mark: $mark. Use one of $_mark_list"
    @check tag_location in (:top, :bottom, :left, :right) "Invalid tag location: $tag_location. Use :top, :bottom, :left, or :right"
    @check tag_alignment in (:horizontal, :vertical, :parallel) "Invalid tag alignment: $tag_alignment. Use :horizontal, :vertical, or :parallel"

    series = DataSeries(kind, X, Y;
        line_style=line_style, dash=dash,
        color=color,
        line_width=line_width,
        mark=mark, mark_size=mark_size,
        mark_color=mark_color, mark_stroke_color=mark_stroke_color,
        label=label, tag=tag, tag_location=tag_location, tag_position=tag_position,
        tag_alignment=tag_alignment,
        bar_width=bar_width, bar_base=bar_base,
        order=order
    )

    if series.color === :default # update colors
        series.color = _default_colors[chart.icolor]
        chart.icolor = mod(chart.icolor, length(_default_colors)) + 1
    end

    if series.order === 0
        series.order = chart.iorder
        chart.iorder += 1
    end

    push!(chart.dataseries, series)

    return series
end


function add_series(chart::Chart, X::AbstractArray, Y::AbstractArray; kwargs...)
    return add_series(chart, :line, X, Y; kwargs...)
end

# Aliases
const add_series = add_series

function add_line(chart::Chart, X::AbstractArray, Y::AbstractArray; kwargs...)
    return add_series(chart, :line, X, Y; kwargs...)
end

function add_scatter(chart::Chart, X::AbstractArray, Y::AbstractArray; kwargs...)
    defaults = (line_style=:none, mark=:circle)
    return add_series(chart, :scatter, X, Y; merge(defaults, kwargs)...)
end

function add_bar(chart::Chart, X::AbstractArray, Y::AbstractArray; kwargs...)
    defaults = (line_style=:none, mark=:none)
    return add_series(chart, :bar, X, Y; merge(defaults, kwargs)...)
end


function configure!(c::Chart)

    length(c.dataseries) > 0 || throw(SerendipException("No dataseries added to the chart"))

    c.outerpad = 0.01 * min(c.width, c.height)
    c.figure_frame = Frame(c.figure_frame.x, c.figure_frame.y, c.width, c.height)

    configure!(c, c.xaxis, c.yaxis)

    if c.aspect_ratio == :equal
        plot_frame = _chart_plot_frame(c)
        xmin, xmax = c.xaxis.limits
        ymin, ymax = c.yaxis.limits
        r = min(plot_frame.width / (xmax - xmin), plot_frame.height / (ymax - ymin))
        dx = 0.5 * (plot_frame.width / r - (xmax - xmin))
        dy = 0.5 * (plot_frame.height / r - (ymax - ymin))

        c.xaxis.limits = [xmin - dx, xmax + dx]
        c.yaxis.limits = [ymin - dy, ymax + dy]
        if !c.xaxis.manual_ticks
            c.xaxis.ticks = Float64[]
        end
        if !c.yaxis.manual_ticks
            c.yaxis.ticks = Float64[]
        end
        if !c.xaxis.manual_tick_labels
            c.xaxis.tick_labels = String[]
        end
        if !c.yaxis.manual_tick_labels
            c.yaxis.tick_labels = String[]
        end

        configure!(c.xaxis)
        configure!(c.yaxis)
    end

    _assign_chart_frames!(c)

end


function configure!(c::Chart, canvas::Canvas)
    canvas.limits = [c.xaxis.limits[1], c.yaxis.limits[1], c.xaxis.limits[2], c.yaxis.limits[2]]
end


function configure!(chart::Chart, xax::Axis, yax::Axis)

    for ax in (xax, yax)
        if ax.auto_limits
            extent = ax.manual_ticks ? collect(extrema(ax.ticks)) : _chart_axis_data_extent(chart, ax)
            ax.limits = compute_auto_limits(extent; mult=ax.mult)
        end
    end

    configure!(xax)
    configure!(yax)

end


function _chart_bar_width(series::DataSeries)
    w = series.bar_width
    if w == 0
        Xu = unique(sort(collect(series.X)))
        if length(Xu) > 1
            w = 0.56 * abs(minimum(diff(Xu)))
        else
            xspan = abs(maximum(series.X) - minimum(series.X))
            w = xspan > 0 ? 0.035 * xspan : 1.0
        end
    end
    return w
end


function _chart_axis_data_extent(chart::Chart, ax::Axis)
    lower = Inf
    upper = -Inf

    for series in chart.dataseries
        series isa DataSeries || continue
        if ax.direction == :horizontal
            isempty(series.X) && continue
            if series.kind == :bar
                w = _chart_bar_width(series)
                lower = min(lower, minimum(series.X) - 0.5 * w)
                upper = max(upper, maximum(series.X) + 0.5 * w)
            else
                lower = min(lower, minimum(series.X))
                upper = max(upper, maximum(series.X))
            end
        else
            isempty(series.Y) && continue
            if series.kind == :bar
                base = series.bar_base
                lower = min(lower, minimum(min.(series.Y, base)))
                upper = max(upper, maximum(max.(series.Y, base)))
            else
                lower = min(lower, minimum(series.Y))
                upper = max(upper, maximum(series.Y))
            end
        end
    end

    return isfinite(lower) && isfinite(upper) ? [lower, upper] : [0.0, 1.0]
end


function configure!(c::Chart, legend::Legend)
    legend.handle_length = 1.9 * legend.font_size
    legend.row_sep = 0.3 * legend.font_size
    legend.col_sep = 1.5 * legend.font_size
    legend.inner_pad = 1.5 * legend.row_sep
    legend.outer_pad = legend.inner_pad

    plots = [p for p in c.dataseries if p.label != ""]

    nlabels     = length(plots)
    label_width = maximum(getsize(plot.label, legend.font_size)[1] for plot in plots)
    label_heigh = maximum(getsize(plot.label, legend.font_size)[2] for plot in plots)

    handle_length = legend.handle_length
    row_sep       = legend.row_sep
    inner_pad     = legend.inner_pad
    col_sep       = legend.col_sep
    ncols         = legend.ncols

    nrows = ceil(Int, nlabels / legend.ncols)

    col_witdhs = zeros(ncols)
    for k in 1:length(plots)
        j = k % ncols == 0 ? ncols : k % ncols # column
        item_width = handle_length + 2 * inner_pad + label_width
        col_witdhs[j] = max(col_witdhs[j], item_width)
    end

    legend.height = nrows * label_heigh + (nrows - 1) * row_sep + 2 * inner_pad
    legend.width = sum(col_witdhs) + (ncols - 1) * col_sep + 2 * inner_pad

end


function _chart_title_font_size(c::Chart)
    return 1.2 * c.xaxis.font_size
end


function _chart_title_height(c::Chart)
    isempty(c.title_box.text) && return 0.0
    surf = CairoImageSurface(4, 4, Cairo.FORMAT_ARGB32)
    cc = CairoContext(surf)
    select_font_face(cc, get_font(c.xaxis.font), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(cc, _chart_title_font_size(c))
    return getsize(cc, c.title_box.text, _chart_title_font_size(c))[2]
end


function _chart_has_legend(c::Chart)
    return any(ds.label != "" for ds in c.dataseries)
end


function _chart_plot_frame(c::Chart)
    left_margin = c.outerpad + c.yaxis.width
    right_margin = c.outerpad
    top_margin = c.outerpad + axis_top_overhang(c.yaxis)
    bottom_margin = c.outerpad + c.xaxis.height
    title_height = _chart_title_height(c)
    title_gap = isempty(c.title_box.text) ? 0.0 : 0.6 * c.xaxis.font_size

    top_margin += title_height + title_gap

    if _chart_has_legend(c)
        legend = c.legend
        if legend.location in (:outer_top_right, :outer_right, :outer_bottom_right)
            right_margin += c.outerpad + legend.width
        elseif legend.location in (:outer_top_left, :outer_left, :outer_bottom_left)
            left_margin += c.outerpad + legend.width
        elseif legend.location == :outer_top
            top_margin += legend.height + c.outerpad
        elseif legend.location == :outer_bottom
            bottom_margin += legend.height + c.outerpad
        end
    end

    width = c.width - left_margin - right_margin
    height = c.height - top_margin - bottom_margin
    frame = Frame(c.figure_frame.x + left_margin, c.figure_frame.y + top_margin, width, height)

    # Ensure the last horizontal tick label fits inside the figure frame.
    right_gap = minimum((frame.x + frame.width) - (frame.x + frame.width / (c.xaxis.limits[2] - c.xaxis.limits[1]) * (tick - c.xaxis.limits[1])) for tick in c.xaxis.ticks)
    extra_right = max(0.0, getsize(c.xaxis.tick_labels[end], c.xaxis.font_size)[1] / 2 - right_gap)

    # Ensure the top-most vertical tick label fits inside the figure frame.
    top_gap = minimum(frame.y + frame.height / (c.yaxis.limits[2] - c.yaxis.limits[1]) * (c.yaxis.limits[2] - tick) - frame.y for tick in c.yaxis.ticks)
    extra_top = max(0.0, getsize(c.yaxis.tick_labels[end], c.yaxis.font_size)[2] / 2 - top_gap)

    return Frame(c.figure_frame.x + left_margin, c.figure_frame.y + top_margin + extra_top, width - extra_right, height - extra_top)
end


function _assign_chart_frames!(c::Chart)
    _chart_has_legend(c) && configure!(c, c.legend)

    plot_frame = _chart_plot_frame(c)
    @check plot_frame.width > 0 && plot_frame.height > 0 "Chart: insufficient space for plot area"

    c.canvas.frame = plot_frame
    configure!(c, c.canvas)

    c.xaxis.width = plot_frame.width
    c.xaxis.frame = Frame(plot_frame.x, plot_frame.y + plot_frame.height, plot_frame.width, c.xaxis.height)

    c.yaxis.height = plot_frame.height
    c.yaxis.frame = Frame(plot_frame.x - c.yaxis.width, plot_frame.y, c.yaxis.width, plot_frame.height)

    c.left_items = FigureComponent[c.yaxis]
    c.right_items = FigureComponent[]
    c.top_items = FigureComponent[]
    c.bottom_items = FigureComponent[c.xaxis]
    c.overlay_items = FigureComponent[a for a in c.annotations]

    if !isempty(c.title_box.text)
        title_height = _chart_title_height(c)
        y = c.figure_frame.y + c.outerpad
        c.title_box.frame = Frame(plot_frame.x, y, plot_frame.width, title_height)
        c.title_box.angle = 0.0
    else
        c.title_box.frame = Frame()
    end

    if _chart_has_legend(c)
        _assign_legend_frame!(c, c.legend)
    else
        c.legend.frame = Frame()
    end
end


function _assign_legend_frame!(c::Chart, legend::Legend)
    plot = c.canvas.frame
    outer_pad = legend.outer_pad

    if legend.location in (:top_right, :right, :bottom_right)
        x1 = plot.x + plot.width - outer_pad - legend.width
    elseif legend.location in (:top, :bottom, :outer_top, :outer_bottom)
        x1 = plot.x + 0.5 * (plot.width - legend.width)
    elseif legend.location in (:top_left, :left, :bottom_left)
        x1 = plot.x + outer_pad
    elseif legend.location in (:outer_top_left, :outer_left, :outer_bottom_left)
        x1 = c.figure_frame.x + c.outerpad
    elseif legend.location in (:outer_top_right, :outer_right, :outer_bottom_right)
        x1 = c.figure_frame.x + c.width - legend.width - c.outerpad
    else
        error("Chart: unsupported legend location $(legend.location)")
    end

    if legend.location in (:top_left, :top, :top_right)
        y1 = plot.y + outer_pad
    elseif legend.location in (:left, :right, :outer_left, :outer_right)
        y1 = plot.y + 0.5 * (plot.height - legend.height)
    elseif legend.location in (:bottom_left, :bottom, :bottom_right)
        y1 = plot.y + plot.height - outer_pad - legend.height
    elseif legend.location == :outer_top
        y1 = c.figure_frame.y + c.outerpad + _chart_title_height(c) + (isempty(c.title_box.text) ? 0.0 : 0.6 * c.xaxis.font_size)
    elseif legend.location == :outer_bottom
        y1 = c.figure_frame.y + c.height - legend.height - c.outerpad
    elseif legend.location in (:outer_top_left, :outer_top_right)
        y1 = plot.y
    elseif legend.location in (:outer_bottom_left, :outer_bottom_right)
        y1 = plot.y + plot.height - legend.height
    else
        error("Chart: unsupported legend location $(legend.location)")
    end

    legend.frame = Frame(x1, y1, legend.width, legend.height)

    if legend.location in (:outer_top_right, :outer_right, :outer_bottom_right)
        push!(c.right_items, legend)
    elseif legend.location in (:outer_top_left, :outer_left, :outer_bottom_left)
        push!(c.left_items, legend)
    elseif legend.location in (:outer_top,)
        push!(c.top_items, legend)
    elseif legend.location in (:outer_bottom,)
        push!(c.bottom_items, legend)
    else
        push!(c.overlay_items, legend)
    end
end


function draw!(c::Chart, ctx::RenderContext, canvas::Canvas)
    # draw grid
    cairo_ctx = ctx.cairo_ctx
    reset_matrix!(ctx)
    set_source_rgba(cairo_ctx, rgba(_chart_gridline_color(c, ctx))...)
    set_line_width(cairo_ctx, 0.2 * ctx.width_scale)
    x0 = canvas.frame.x
    y0 = canvas.frame.y
    x1 = canvas.frame.x + canvas.frame.width
    y1 = canvas.frame.y + canvas.frame.height

    xmin, xmax = c.xaxis.limits
    for x in c.xaxis.ticks
        min(xmax, xmin) <= x <= max(xmax, xmin) || continue
        xc = x0 + canvas.frame.width / (xmax - xmin) * (x - xmin)
        move_to(cairo_ctx, xc, y0)
        line_to(cairo_ctx, xc, y1)
        stroke(cairo_ctx)
    end

    ymin, ymax = c.yaxis.limits
    for y in c.yaxis.ticks
        min(ymax, ymin) <= y <= max(ymax, ymin) || continue
        yc = y0 + canvas.frame.height / (ymax - ymin) * (ymax - y)
        move_to(cairo_ctx, x0, yc)
        line_to(cairo_ctx, x1, yc)
        stroke(cairo_ctx)
    end

    # draw border
    set_source_rgb(cairo_ctx, 0.0, 0.0, 0.0)
    set_line_width(cairo_ctx, 0.5 * ctx.width_scale)
    rectangle(cairo_ctx, x0, y0, canvas.frame.width, canvas.frame.height)
    stroke(cairo_ctx)
end


function draw!(chart::Chart, ctx::RenderContext, p::DataSeries)
    cairo_ctx = ctx.cairo_ctx

    p.mark_color = p.mark_color == :default ? p.color : p.mark_color
    p.mark_stroke_color = p.mark_stroke_color == :default ? p.color : p.mark_stroke_color

    # p.mark_color = get_color(p.mark_color, p.color)
    # p.mark_stroke_color = get_color(p.mark_stroke_color, p.color)

    reset_matrix!(ctx)
    set_source_rgb(cairo_ctx, rgb(p.color)...)
    set_line_width(cairo_ctx, p.line_width * ctx.width_scale)
    set_line_join(cairo_ctx, Cairo.CAIRO_LINE_JOIN_ROUND)

    # Draw lines
    new_path(cairo_ctx)
    n = length(p.X)
    X = p.X * chart.xaxis.mult
    Y = p.Y * chart.yaxis.mult
    if p.kind == :bar
        xmin, xmax = chart.xaxis.limits
        xspan = abs(xmax - xmin)

        w = p.bar_width
        if w == 0
            if n > 1
                Xu = unique(sort(collect(X)))
                if length(Xu) > 1
                    dx = minimum(diff(Xu))
                    w = 0.56 * abs(dx)
                end
            end
            w == 0 && (w = xspan > 0 ? 0.035 * xspan : 1.0)
        end
        w *= abs(chart.xaxis.mult)
        base = p.bar_base * chart.yaxis.mult

        for (x, y) in zip(X, Y)
            y0 = base
            h = y - y0
            xleft = x - 0.5 * w
            ytop = h >= 0 ? y : y0
            rect_x, rect_y = data2user(chart.canvas, xleft, ytop)
            xden = abs(chart.xaxis.limits[2] - chart.xaxis.limits[1])
            yden = abs(chart.yaxis.limits[2] - chart.yaxis.limits[1])
            rect_w = xden > 0 ? chart.canvas.frame.width / xden * w : 0.0
            rect_h = yden > 0 ? chart.canvas.frame.height / yden * abs(h) : 0.0

            rectangle(cairo_ctx, rect_x, rect_y, rect_w, rect_h)
            fill_preserve(cairo_ctx)
            set_source_rgb(cairo_ctx, 0.0, 0.0, 0.0)
            set_line_width(cairo_ctx, p.line_width * ctx.width_scale)
            stroke(cairo_ctx)
            set_source_rgb(cairo_ctx, rgb(p.color)...)
        end
        return
    end

    if p.kind != :scatter && p.line_style !== :none
        x1, y1 = data2user(chart.canvas, X[1], Y[1])

        if p.line_style == :solid
            move_to(cairo_ctx, x1, y1)
            for i in 2:n
                x, y = data2user(chart.canvas, X[i], Y[i])
                line_to(cairo_ctx, x, y)
            end
            stroke(cairo_ctx)
        else # dashed
            len = sum(p.dash)
            offset = 0.0
            set_dash(cairo_ctx, p.dash, offset)
            move_to(cairo_ctx, x1, y1)
            for i in 2:n
                x, y = data2user(chart.canvas, X[i], Y[i])
                line_to(cairo_ctx, x, y)
                offset = mod(offset + norm((x1 - x, y1 - y)), len)
                set_dash(cairo_ctx, p.dash, offset)
                x1, y1 = x, y
            end
            stroke(cairo_ctx)
            set_dash(cairo_ctx, Float64[])
        end
    end

    # Draw marks
    for (x, y) in zip(X, Y)
        x, y = data2user(chart.canvas, x, y)
        draw_mark(cairo_ctx, x, y, p.mark, p.mark_size, p.mark_color, p.mark_stroke_color)
    end

    # Draw tag
    if p.tag != ""
        len = 0.0
        L = [len] # lengths
        for i in 2:length(X)
            len += norm((X[i] - X[i-1], Y[i] - Y[i-1]))
            push!(L, len)
        end
        lpos = p.tag_position * len # length to position

        i = findfirst(z -> z > lpos, L)
        i = min(i, length(L) - 1)

        # location coordinates
        x = X[i] + (lpos - L[i]) / (L[i+1] - L[i]) * (X[i+1] - X[i])
        y = Y[i] + (lpos - L[i]) / (L[i+1] - L[i]) * (Y[i+1] - Y[i])

        # location coordinates in user units
        x, y = data2user(chart.canvas, x, y)
        x1, y1 = data2user(chart.canvas, X[i], Y[i])
        x2, y2 = data2user(chart.canvas, X[i+1], Y[i+1])
        α = -atand(y2 - y1, x2 - x1) # tilt

        # pads
        pad = chart.xaxis.font_size * 0.3

        dx = pad * abs(sind(α))
        dy = pad * abs(cosd(α))

        # Default location "top"
        if p.tag_location == :top
            va = "bottom"
            if 0 < α <= 90 || -180 < α <= -90
                ha = "right"
                dx, dy = -dx, -dy
            else
                ha = "left"
                dy = -dy
            end
        else
            va = "top"
            if 0 < α <= 90 || -180 < α <= -90
                ha = "left"
            else
                ha = "right"
                dx = -dx
            end
        end

        if p.tag_alignment == :parallel
            ha = "center"
            dx = 0.0
            dy = p.tag_location == :top ? -pad : 0.0
        else
            α = 0.0
        end

        set_font_size(cairo_ctx, chart.xaxis.font_size * 0.9)
        font = get_font(chart.xaxis.font)
        select_font_face(cairo_ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
        set_source_rgb(cairo_ctx, 0, 0, 0)
        draw_text(cairo_ctx, x + dx, y + dy, p.tag, halign=ha, valign=va, angle=α)
    end

end


function draw!(c::Chart, ctx::RenderContext, legend::Legend)
    cairo_ctx = ctx.cairo_ctx

    plots = [p for p in c.dataseries if p.label != ""]

    set_font_size(cairo_ctx, legend.font_size)
    font = get_font(legend.font)
    select_font_face(cairo_ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)

    handle_length = legend.handle_length
    row_sep = legend.row_sep
    inner_pad = legend.inner_pad
    outer_pad = legend.outer_pad
    col_sep = legend.col_sep
    ncols = legend.ncols

    # update the width
    col_witdhs = zeros(ncols)
    for (k, plot) in enumerate(plots)
        j = k % ncols == 0 ? ncols : k % ncols # column
        label_width = getsize(cairo_ctx, plot.label, legend.font_size)[1]
        item_width = handle_length + 2 * inner_pad + label_width
        col_witdhs[j] = max(col_witdhs[j], item_width)
    end

    # legend.height = nrows*label_heigh + (nrows-1)*row_sep + 2*inner_pad
    legend.width = sum(col_witdhs) + (ncols - 1) * col_sep + 2 * inner_pad

    x1 = legend.frame.x
    y1 = legend.frame.y
    x2 = legend.frame.x + legend.frame.width
    y2 = legend.frame.y + legend.frame.height

    reset_matrix!(ctx)

    # draw rounded rectangle
    r = 0.02 * min(c.canvas.frame.width, c.canvas.frame.height)
    move_to(cairo_ctx, x1, y1 + r)
    line_to(cairo_ctx, x1, y2 - r)
    curve_to(cairo_ctx, x1, y2, x1, y2, x1 + r, y2)
    line_to(cairo_ctx, x2 - r, y2)
    curve_to(cairo_ctx, x2, y2, x2, y2, x2, y2 - r)
    line_to(cairo_ctx, x2, y1 + r)
    curve_to(cairo_ctx, x2, y1, x2, y1, x2 - r, y1)
    line_to(cairo_ctx, x1 + r, y1)
    curve_to(cairo_ctx, x1, y1, x1, y1, x1, y1 + r)
    close_path(cairo_ctx)
    legend_background = something(legend.background, ctx.background, Color(:white))
    set_source_rgba(cairo_ctx, rgba(legend_background)...)
    fill_preserve(cairo_ctx)
    set_source_rgb(cairo_ctx, 0, 0, 0) # black
    set_line_width(cairo_ctx, 0.4 * ctx.width_scale)
    stroke(cairo_ctx)

    # draw labels
    label_heigh = maximum(getsize(plot.label, legend.font_size)[2] for plot in plots)

    for (k, plot) in enumerate(plots)
        i = ceil(Int, k / ncols)  # line
        j = k % ncols == 0 ? ncols : k % ncols # column
        x2 = x1 + inner_pad + sum(col_witdhs[1:j-1]) + (j - 1) * col_sep

        y2 = y1 + inner_pad + label_heigh / 2 + (i - 1) * (label_heigh + row_sep)

        set_source_rgb(cairo_ctx, rgb(plot.color)...)
        if plot.kind == :bar
            hbar = 0.6 * legend.font_size
            rectangle(cairo_ctx, x2, y2 - 0.5 * hbar, handle_length, hbar)
            fill_preserve(cairo_ctx)
            set_source_rgb(cairo_ctx, 0.0, 0.0, 0.0)
            set_line_width(cairo_ctx, max(plot.line_width, 0.4) * ctx.width_scale)
            stroke(cairo_ctx)
            set_source_rgb(cairo_ctx, rgb(plot.color)...)
        elseif plot.line_style != :none
            move_to(cairo_ctx, x2, y2)
            rel_line_to(cairo_ctx, handle_length, 0)
            set_line_width(cairo_ctx, plot.line_width * ctx.width_scale)
            plot.line_style != :solid && set_dash(cairo_ctx, plot.dash)
            stroke(cairo_ctx)
            set_dash(cairo_ctx, Float64[])
        end

        # draw mark
        if plot.kind != :bar
            x = x2 + handle_length / 2
            draw_mark(cairo_ctx, x, y2, plot.mark, plot.mark_size, plot.mark_color, plot.mark_stroke_color)
        end

        # draw label
        x = x2 + handle_length + 2 * inner_pad
        y = y2

        set_source_rgb(cairo_ctx, 0, 0, 0)
        draw_text(cairo_ctx, x, y, plot.label, halign="left", valign="center", angle=0)
    end

end


function draw_background!(c::Chart, ctx::RenderContext)
    _draw_figure_background!(ctx, c.figure_frame, ctx.background)
end

function draw_contents!(c::Chart, ctx::RenderContext)
    cairo_ctx = ctx.cairo_ctx
    reset_matrix!(ctx)
    # draw canvas grid
    draw!(c, ctx, c.canvas)

    # draw axes
    draw!(c.xaxis, ctx)

    draw!(c.yaxis, ctx)

    # draw plots
    rectangle(cairo_ctx, c.canvas.frame.x, c.canvas.frame.y, c.canvas.frame.width, c.canvas.frame.height)
    Cairo.clip(cairo_ctx)

    # draw dataseries
    sorted = sort(c.dataseries, by=x -> x.order)
    for p in sorted
        draw!(c, ctx, p)
    end
    reset_clip(cairo_ctx)

    if !isempty(c.title_box.text)
        reset_matrix!(ctx)
        select_font_face(cairo_ctx, get_font(c.xaxis.font), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
        set_font_size(cairo_ctx, _chart_title_font_size(c))
        set_source_rgb(cairo_ctx, 0.0, 0.0, 0.0)
        _draw_text_box!(ctx, c.title_box)
    end

    # draw overlay annotations before legend
    for item in c.overlay_items
        item isa Annotation || continue
        draw!(c, ctx, item)
    end

    # draw legend last
    if _chart_has_legend(c)
        draw!(c, ctx, c.legend)
    end
end


function add_annotation(c::Chart, a::Annotation)
    push!(c.annotations, a)
    push!(c.overlay_items, a)
    return a
end
