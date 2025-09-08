# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

"""
    Chart(; kwargs...)

Creates a customizable `Chart` instance for plotting data with labeled axes,
legends, ticks, and optional colorbar support.

# Parameters
- `size::Tuple{Int,Int}=(220,150)` : Chart drawing size in pt. 1 cm = 28.35 pt.
- `font::AbstractString="NewComputerModern"` : Font family for axis and legend text.
- `font_size::Float64=7.0` : Font size for axis labels and ticks.
- `xlimits, xlims::Vector{Float64}` : Limits of the x-axis `[xmin, xmax]`.
- `ylimits, ylims::Vector{Float64}` : Limits of the y-axis `[ymin, ymax]`.
- `aspect_ratio::Symbol` : Ratio of y-unit to x-unit (`:auto` or `:equal`, default = `:auto`).
- `xmult::Float64` : Scaling factor for x-axis values (default = 1.0).
- `ymult::Float64` : Scaling factor for y-axis values (default = 1.0).
- `xbins::Int` : Number of tick bins along the x-axis (default = 7).
- `ybins::Int` : Number of tick bins along the y-axis (default = 6).
- `xlabel::AbstractString` : Label for the x-axis (default = L"\$x\$").
- `ylabel::AbstractString` : Label for the y-axis (default = L"\$y\$").
- `xticks::AbstractArray{Float64}` : Explicit tick positions along the x-axis (default = []).
- `yticks::AbstractArray{Float64}` : Explicit tick positions along the y-axis (default = []).
- `xtick_labels::AbstractArray{String}` : Labels for x-axis ticks (default = []).
- `ytick_labels::AbstractArray{String}` : Labels for y-axis ticks (default = []).
- `legend::Symbol` : Legend location (`:top_right`, `:top_left`, `:bottom_right`, `:bottom_left`, etc.).
- `legend_font_size::Union{Float64,Symbol}` : Font size for legend text (> 0, default = `:font_size`).
- `quiet::Bool` : Suppress console output if `true` (default = `false`).

# Notes
- A `Chart` manages the figure canvas, axes, legend, data series, and annotations.
- Axis properties (limits, bins, ticks, labels) are set at initialization but can be modified later.
- The `legend` is initialized with the chosen font and location.
- `aspect_ratio=:equal` enforces equal unit length on both axes.
- If `quiet=false`, chart information is printed to the console at creation.

# Example
```julia
ch = Chart(size=(300,200),
           xlabel="Time [s]",
           ylabel="Displacement [mm]",
           xlimits=[0.0,10.0],
           ylimits=[-5.0,5.0],
           legend=:bottom_right)
```
"""




"""
    Chart(; 
        size=(220,150), font="NewComputerModern", font_size=7.0,
        xlimits, ylimits, aspect_ratio=:auto,
        xmult=1.0, ymult=1.0, xbins=7, ybins=6,
        xlabel=L"\$x\$", ylabel=L"\$y\$",
        xticks=Float64[], yticks=Float64[],
        xtick_labels=String[], ytick_labels=String[],
        legend=:top_right, legend_font_size=0,
        quiet=false)

Construct a 2D chart figure with axes, legend, and optional tick customization.

# Arguments
- `size::Tuple{Int,Int}`: width × height in points. 1 cm = 28.35 points.
- `font::AbstractString`: font family for axes and legend.
- `font_size::Real`: base font size.
- `xlimits::Vector{<:Real}`, `ylimits::Vector{<:Real}`: axis limits `[min,max]`; `[0,0]` enables auto scaling.
- `aspect_ratio::Symbol`: `:auto` or `:equal`.
- `xmult::Real`, `ymult::Real`: multiplicative factors applied to tick values.
- `xbins::Int`, `ybins::Int`: target number of major ticks.
- `xlabel::AbstractString`, `ylabel::AbstractString`: axis labels.
- `xticks::Vector{<:Real}`, `yticks::Vector{<:Real}`: explicit tick positions; empty vectors enable auto ticks.
- `xtick_labels::Vector{<:AbstractString}`, `ytick_labels::Vector{<:AbstractString}`: custom tick labels; if provided, lengths must match the corresponding tick arrays.
- `legend::Symbol`: legend location (e.g., `:top_right`, `:top_left`, `:bottom_left`, `:none`).
- `legend_font_size::Real`: legend font size; `0` uses `font_size`.
- `quiet::Bool`: suppress constructor log.

# Notes
- Use `add_series` to append data series to the chart.
- Use `save` to export the chart to a file.

# Returns
- A `Chart` object.

# Example
```julia
ch = Chart(size=(300,200),
           xlabel="Time [s]",
           ylabel="Displacement [mm]",
           xlimits=[0.0,10.0],
           ylimits=[-5.0,5.0],
           legend=:bottom_right)
```
"""
mutable struct Chart<:Figure
    width::Float64
    height::Float64
    canvas::Canvas
    xaxis::Axis
    yaxis::Axis
    dataseries::Vector{DataSeries}
    legend::Legend
    annotations::AbstractArray
    # colorbar::Union{FigureComponent, Nothing}

    aspect_ratio::Symbol
    outerpad::Float64
    leftpad::Float64
    rightpad::Float64
    toppad::Float64
    bottompad::Float64
    icolor::Int
    iorder::Int

    quiet::Bool

    function Chart(;
        size=(220,150),
        font="NewComputerModern",
        font_size::Real=7.0,
        xlimits=[0.0,0.0],
        ylimits=[0.0,0.0],
        aspect_ratio=:auto,
        xmult::Real=1.0,
        ymult::Real=1.0,
        xbins::Int=7,
        ybins::Int=6,
        xlabel::AbstractString=L"$x$",
        ylabel::AbstractString=L"$y$",
        xticks::Vector{<:Real}=Float64[],
        yticks::Vector{<:Real}=Float64[],
        xtick_labels::Vector{<:AbstractString}=String[],
        ytick_labels::Vector{<:AbstractString}=String[],
        legend::Symbol=:top_right,
        legend_font_size::Real=0,
        # colorbar=:right,
        # colorbar_scale=0.9,
        # colorbar_label="",
        # colorbar_limits=[0.0, 0.0],
        # colorbar_font_size=7.0,
        quiet=false
    )
        if legend_font_size==0
            legend_font_size = font_size
        end

        @check font_size > 0
        @check legend_font_size > 0
        @check aspect_ratio in (:auto, :equal) "Chart: Invalid aspect_ratio: $aspect_ratio. Use :auto or :equal. Got $(repr(aspect_ratio))."

        width, height = size
        outerpad = 0.01*min(width, height)
        the_legend = Legend(; location=legend, font=font, font_size=legend_font_size, ncols=1)
        xaxis = Axis(direction=:horizontal, limits=xlimits, label=xlabel, font=font, font_size=font_size, ticks=xticks, tick_labels=xtick_labels, mult=xmult, bins=xbins)
        yaxis = Axis(direction= :vertical, limits=ylimits, label=ylabel, font=font, font_size=font_size, ticks=yticks, tick_labels=ytick_labels, mult=ymult, bins=ybins)

        this = new(width, height, Canvas(), xaxis, yaxis, [], the_legend, [],
            aspect_ratio, outerpad, outerpad, outerpad, outerpad, outerpad, 1, 1, quiet)

        if !quiet
            printstyled("Chart figure\n", bold=true, color=:cyan)
            println("  size: $(this.width) x $(this.height) pt")
            println("  axes: $(xlabel) vs $(ylabel)")
        end

        return this
    end
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
- `color = :default` : Series color. `:default` selects from the chart palette cyclically.
- `line_width::Float64 = 0.5` : Line width (> 0).
- `marker::Symbol = :none` : Marker shape.
- `marker_size::Float64 = 2.5` : Marker size (> 0).
- `marker_color = :white` : Marker fill color.
- `marker_stroke_color = :default` : Marker edge color (`:default` follows `color`).
- `label::AbstractString = ""` : Legend label.
- `tag::AbstractString = ""` : On-curve annotation text.
- `tag_location::Symbol = :top` : Relative location of tag (`:top`, `:bottom`, `:left`, `:right`).
- `tag_position::Float64 = 0.5` : Position along the curve in [0,1].
- `tag_alignment::Symbol = :horizontal` : Tag orientation (`:horizontal`, `:vertical`, `:parallel`).
- `order::Int = 0` : Z-order. If `0`, an incremental order is assigned.

# Returns
- The series object.

# Examples
```julia
ch = Chart(size=(300,200), xlabel="Time [s]", ylabel="Displacement [mm]",
           xlimits=[0.0,10.0], ylimits=[-5.0,5.0], legend=:bottom_right)

add_series(ch, 0:0.1:10, sin.(0:0.1:10); label="sin")

add_series(ch, :scatter, x, y;
           marker=:circle, marker_size=3,
           color=:black, label="samples")

add_series(ch, :bar, 1:5, [2,3,1,4,2]; color=:gray, label="counts")
```
"""
function add_series(chart::Chart, kind::Symbol, X::AbstractArray, Y::AbstractArray;
    line_style=:solid, dash=Float64[],
    color=:default,
    line_width=0.5,
    marker=:none, marker_size=2.5,
    marker_color=:white, marker_stroke_color=:default,
    label="", tag="", tag_location=:top, tag_position=0.5,
    tag_alignment=:horizontal,
    order=0
)

    @check line_width > 0 "Line width must be positive"
    @check marker_size > 0 "Marker size must be positive"
    @check tag_position >= 0 && tag_position <= 1 "Tag position must be in [0,1]"
    @check order >= 0 "Order must be non-negative"
    @check kind in (:line, :scatter, :bar) "Invalid series kind: $kind. Use :line, :scatter, or :bar"
    @check length(X) == length(Y) "X and Y must have the same length"
    @check marker in _marker_list "Invalid marker: $marker. Use one of $_default_markers"
    @check tag_location in (:top, :bottom, :left, :right) "Invalid tag location: $tag_location. Use :top, :bottom, :left, or :right"
    @check tag_alignment in (:horizontal, :vertical, :parallel) "Invalid tag alignment: $tag_alignment. Use :horizontal, :vertical, or :parallel"

    series = DataSeries(kind, X, Y;
        line_style=line_style, dash=dash,
        color=color,
        line_width=line_width,
        marker=marker, marker_size=marker_size,
        marker_color=marker_color, marker_stroke_color=marker_stroke_color,
        label=label, tag=tag, tag_location=tag_location, tag_position=tag_position,
        tag_alignment=tag_alignment,
        order=order
    )

    if series.color===:default # update colors
        series.color = _default_colors[chart.icolor]
        chart.icolor = mod(chart.icolor, length(_default_colors)) + 1
    end

    if series.order===0
        series.order = chart.iorder
        chart.iorder += 1
    end

    push!(chart.dataseries, series)

    return series
end


function add_series(chart::Chart, X::AbstractArray, Y::AbstractArray; kwargs...)
    return add_series(chart, :line, X, Y; kwargs...)
end


function configure!(c::Chart)

    length(c.dataseries) > 0 || throw(SerendipException("No dataseries added to the chart"))

    c.outerpad  = 0.01*min(c.width, c.height)
    c.leftpad   = c.outerpad
    c.rightpad  = c.outerpad
    c.toppad    = c.outerpad
    c.bottompad = c.outerpad

    # configure legend
    if any( ds.label!="" for ds in c.dataseries )
        configure!(c, c.legend)
    end

    # configure axes, may change chart pads
    configure!(c, c.xaxis, c.yaxis)

    # configure the canvas
    configure!(c, c.canvas)

    # for p in c.dataseries
    #     configure!(c, p)
    # end

end


function configure!(c::Chart, canvas::Canvas)
    xmin, xmax = c.xaxis.limits
    ymin, ymax = c.yaxis.limits
    # canvas.limits = [ xmin, ymin, xmax, ymax ]

    if c.aspect_ratio==:equal
        # compute extra limits
        width = c.width - c.yaxis.width - c.leftpad - c.rightpad
        height = c.height - c.xaxis.height - c.toppad - c.rightpad
        r = min(width/(xmax-xmin), height/(ymax-ymin))
        dx = 0.5*(width/r - (xmax-xmin))
        dy = 0.5*(height/r - (ymax-ymin))

        # update limits
        c.xaxis.limits = [ xmin-dx, xmax+dx ]
        c.yaxis.limits = [ ymin-dy, ymax+dy ]

        # force recompute ticks
        c.xaxis.ticks = []
        c.yaxis.ticks = []

        # reconfigure axes
        configure!(c, c.xaxis, c.yaxis)
        xmin, xmax = c.xaxis.limits
        ymin, ymax = c.yaxis.limits

        # udpa
    end

    canvas.width = c.width - c.yaxis.width - c.leftpad - c.rightpad
    canvas.height = c.height - c.xaxis.height - c.toppad - c.bottompad
    canvas.box = [ c.leftpad + c.yaxis.width, c.toppad, c.width-c.rightpad, c.height - c.xaxis.height-c.bottompad ]
    canvas.limits = [ xmin, ymin, xmax, ymax ]
end


function configure!(chart::Chart, xax::Axis, yax::Axis)

    # check limits
    for ax in (xax, yax)
        if ax.limits==[0.0,0.0]
            if length(ax.ticks)==0
                limits = [Inf, -Inf]
                for p in chart.dataseries
                    p isa DataSeries || continue
                    if ax.direction==:horizontal
                        # if p.x!==nothing
                            # limits[1] = min(limits[1], p.x)
                            # limits[2] = max(limits[2], p.x)
                        # end
                        if length(p.X)>0
                            limits[1] = min(limits[1], minimum(p.X))
                            limits[2] = max(limits[2], maximum(p.X))
                        end
                    else
                        # if p.y!==nothing
                            # limits[1] = min(limits[1], p.y)
                            # limits[2] = max(limits[2], p.y)
                        # end
                        if length(p.Y)>0
                            limits[1] = min(limits[1], minimum(p.Y))
                            limits[2] = max(limits[2], maximum(p.Y))
                        end

                    end
                end
                if limits[1]==limits[2]
                    limits = [limits[1]-0.1, limits[1]+0.1]
                end
            else
                limits = collect(extrema(chart.args.xticks))
            end

            # extend limits
            f = ax.direction == :horizontal ? 0.03 : 0.03*chart.width/chart.height
            dx = f*(limits[2]-limits[1])
            limits = [ limits[1]-dx, limits[2]+dx ]
            limits = [ limits[1]*ax.mult, limits[2]*ax.mult ]
            ax.limits = limits
        end
    end

    configure!(xax)
    configure!(yax)

    # set width and height of axes
    width, height = chart.width, chart.height
    xax.width  = width - yax.width - chart.leftpad - chart.rightpad
    yax.height = height - xax.height - chart.toppad - chart.bottompad

    # update chart.rightpad if required
    label_width = getsize(xax.tick_labels[end], xax.font_size)[1]
    xdist = xax.width*(xax.limits[2]-xax.ticks[end])/(xax.limits[2]-xax.limits[1]) # distance of the right most tick to the right side of axis
    if xdist-label_width/2 < 0
        chart.rightpad += label_width/2 - xdist
        xax.width = width - yax.width - chart.leftpad - chart.rightpad
    end

    # update chart.toppad if required
    label_height = getsize(yax.tick_labels[end], yax.font_size)[2]
    ydist = yax.height * (yax.limits[2] - yax.ticks[end])/(yax.limits[2]-yax.limits[1]) # distance of the upper most tick to the top side of axis
    if ydist-label_height/2 < 0
        chart.toppad += label_height/2 - ydist
        yax.height = height - xax.height - chart.toppad - chart.bottompad
    end

end


# function configure!(chart::Chart, p::DataSeries)
#     xmin, ymin, xmax, ymax = chart.canvas.limits
#     # if p.x !== nothing
#         # p.X = [ p.x, p.x ]
#         # p.Y = [ ymin, ymax ]
#     # elseif p.y !== nothing
#         # p.X = [ xmin, xmax ]
#         # p.Y = [ p.y, p.y ]
#     # end
# end


function configure!(c::Chart, legend::Legend)
    legend.handle_length = 1.9*legend.font_size
    legend.row_sep = 0.3*legend.font_size
    legend.col_sep = 1.5*legend.font_size
    legend.inner_pad = 1.5*legend.row_sep
    legend.outer_pad = legend.inner_pad

    plots = [ p for p in c.dataseries if p.label != ""]

    nlabels = length(plots)
    label_width = maximum( getsize(plot.label, legend.font_size)[1] for plot in plots )
    label_heigh = maximum( getsize(plot.label, legend.font_size)[2] for plot in plots )

    handle_length = legend.handle_length
    row_sep       = legend.row_sep
    inner_pad     = legend.inner_pad
    col_sep       = legend.col_sep
    ncols         = legend.ncols

    nrows = ceil(Int, nlabels/legend.ncols)

    col_witdhs = zeros(ncols)
    for k in 1:length(plots)
        j = k%ncols==0 ? ncols : k%ncols # column
        item_width = handle_length + 2*inner_pad + label_width
        col_witdhs[j] = max(col_witdhs[j], item_width)
    end

    legend.height = nrows*label_heigh + (nrows-1)*row_sep + 2*inner_pad
    legend.width = sum(col_witdhs) + (ncols-1)*col_sep + 2*inner_pad

    if legend.location in (:outer_top_right, :outer_right, :outer_bottom_right)
        c.rightpad += legend.width + c.outerpad
    elseif legend.location in ( :outer_top_left, :outer_left, :outer_bottom_left)
        c.leftpad += legend.width + c.outerpad
    elseif legend.location == :outer_top
        c.toppad += legend.height + 2*c.outerpad
    elseif legend.location == :outer_bottom
        c.bottompad += legend.height + 2*c.outerpad
    end

end


function draw!(c::Chart, ctx::CairoContext, canvas::Canvas)
    # draw grid
    set_source_rgb(ctx, 0.9, 0.9, 0.9) # gray
    set_line_width(ctx, 0.2)

    xmin, xmax = c.xaxis.limits
    for x in c.xaxis.ticks
        min(xmax, xmin) <= x<= max(xmax,xmin) || continue
        x1 = canvas.box[1] + c.xaxis.width/(xmax-xmin)*(x-xmin)
        y1 = canvas.box[2]
        y2 = canvas.box[4]
        move_to(ctx, x1, y1); line_to(ctx, x1, y2); stroke(ctx)
    end

    ymin, ymax = c.yaxis.limits
    for y in c.yaxis.ticks
        min(ymax, ymin) <= y<= max(ymax,ymin) || continue
        y1 = canvas.box[2] + c.yaxis.height/(ymax-ymin)*(ymax-y)
        x1 = canvas.box[1]
        x2 = canvas.box[3]
        move_to(ctx, x1, y1); line_to(ctx, x2, y1); stroke(ctx)
    end

    # draw border
    set_source_rgb(ctx, 0.0, 0.0, 0.0)
    set_line_width(ctx, 0.5)
    x, y = canvas.box[1:2]
    w, h = canvas.box[3:4] - canvas.box[1:2]
    rectangle(ctx, x, y, w, h)
    stroke(ctx)
end


function draw!(chart::Chart, ctx::CairoContext, p::DataSeries)

    p.marker_color = p.marker_color==:default ? p.color : p.marker_color
    p.marker_stroke_color = p.marker_stroke_color==:default ? p.color : p.marker_stroke_color

    # p.marker_color = get_color(p.marker_color, p.color)
    # p.marker_stroke_color = get_color(p.marker_stroke_color, p.color)

    set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
    set_source_rgb(ctx, rgb(p.color)...)
    set_line_width(ctx, p.line_width)
    set_line_join(ctx, Cairo.CAIRO_LINE_JOIN_ROUND)

    # Draw lines
    new_path(ctx)
    n = length(p.X)
    X = p.X*chart.xaxis.mult
    Y = p.Y*chart.yaxis.mult

    if p.line_style !== :none
        x1, y1 = data2user(chart.canvas, X[1], Y[1])

        if p.line_style==:solid
            move_to(ctx, x1, y1)
            for i in 2:n
                x, y = data2user(chart.canvas, X[i], Y[i])
                line_to(ctx, x, y);
            end
            stroke(ctx)
        else # dashed
            len = sum(p.dash)
            offset = 0.0
            set_dash(ctx, p.dash, offset)
            move_to(ctx, x1, y1)
            for i in 2:n
                x, y = data2user(chart.canvas, X[i], Y[i])
                line_to(ctx, x, y);
                offset = mod(offset + norm((x1-x,y1-y)), len)
                set_dash(ctx, p.dash, offset)
                x1, y1 = x, y
            end
            stroke(ctx)
            set_dash(ctx, Float64[])
        end
    end

    # Draw markers
    for (x,y) in zip(X, Y)
        x, y = data2user(chart.canvas, x, y)
        draw_marker(ctx, x, y, p.marker, p.marker_size, p.marker_color, p.marker_stroke_color)
    end

    # Draw tag
    if p.tag!=""
        len = 0.0
        L = [ len ] # lengths
        for i in 2:length(X)
            len += norm((X[i]-X[i-1], Y[i]-Y[i-1]))
            push!(L, len)
        end
        lpos = p.tagpos*len # length to position

        i = findfirst(z->z>lpos, L)
        i = min(i, length(L)-1)

        # location coordinates
        x  = X[i] + (lpos-L[i])/(L[i+1]-L[i])*(X[i+1]-X[i])
        y  = Y[i] + (lpos-L[i])/(L[i+1]-L[i])*(Y[i+1]-Y[i])

        # location coordinates in user units
        x, y = data2user(chart.canvas, x, y)
        x1, y1 = data2user(chart.canvas, X[i], Y[i])
        x2, y2 = data2user(chart.canvas, X[i+1], Y[i+1])
        α = -atand(y2-y1, x2-x1) # tilt

        # pads
        pad = chart.args.font_size*0.3

        dx = pad*abs(sind(α))
        dy = pad*abs(cosd(α))

        # Default location "top"
        if p.tagloc==:top
            va = "bottom"
            if 0<α<= 90 || -180 <α<= -90
                ha = "right"
                dx, dy = -dx, -dy
            else
                ha = "left"
                dy = -dy
            end
        else
            va = "top"
            if 0<α<=90 || -180<α<=-90
                ha = "left"
            else
                ha = "right"
                dx = -dx
            end
        end

        if p.tagalong
            ha = "center"
            dx = 0.0
            dy = p.tagloc==:top ? -pad : 0.0
        else
            α = 0.0
        end

        set_font_size(ctx, chart.args.font_size*0.9)
        font = get_font(chart.args.font)
        select_font_face(ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
        set_source_rgb(ctx, 0, 0, 0)
        draw_text(ctx, x+dx, y+dy, p.tag, halign=ha, valign=va, angle=α)
    end

end


function draw!(c::Chart, ctx::CairoContext, legend::Legend)

    plots = [ p for p in c.dataseries if p.label != ""]

    set_font_size(ctx, legend.font_size)
    font = get_font(legend.font)
    select_font_face(ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )

    handle_length = legend.handle_length
    row_sep       = legend.row_sep
    inner_pad     = legend.inner_pad
    outer_pad     = legend.outer_pad
    col_sep       = legend.col_sep
    ncols         = legend.ncols

    # update the width
    col_witdhs = zeros(ncols)
    for (k, plot) in enumerate(plots)
        j = k%ncols==0 ? ncols : k%ncols # column
        label_width = getsize(ctx, plot.label, legend.font_size)[1]
        item_width = handle_length + 2*inner_pad + label_width
        col_witdhs[j] = max(col_witdhs[j], item_width)
    end

    # legend.height = nrows*label_heigh + (nrows-1)*row_sep + 2*inner_pad
    legend.width = sum(col_witdhs) + (ncols-1)*col_sep + 2*inner_pad

    # set legend location
    if legend.location in (:top_right, :right, :bottom_right)
        x1 = c.canvas.box[3] - outer_pad - legend.width
    elseif legend.location in (:top, :bottom, :outer_top, :outer_bottom)
        x1 = 0.5*(c.canvas.box[1]+c.canvas.box[3]) - legend.width/2
    elseif legend.location in (:top_left, :left, :bottom_left)
        x1 = c.canvas.box[1] + outer_pad
    elseif legend.location in (:outer_top_left, :outer_left, :outer_bottom_left)
        x1 = c.outerpad
    elseif legend.location in (:outer_top_right, :outer_right, :outer_bottom_right)
        x1 = c.width - legend.width - c.outerpad
    end

    if legend.location in (:top_left, :top, :top_right)
        y1 = c.canvas.box[2] + outer_pad
    elseif legend.location in (:left, :right, :outer_left, :outer_right)
        y1 = 0.5*(c.canvas.box[2]+c.canvas.box[4]) - legend.height/2
    elseif legend.location in (:bottom_left, :bottom, :bottom_right)
        y1 = c.canvas.box[4] - outer_pad - legend.height
    elseif legend.location == :outer_top
        y1 = c.outerpad
    elseif legend.location == :outer_bottom
        y1 = c.height - legend.height - c.outerpad
    elseif legend.location in (:outer_top_left, :outer_top_right)
        y1 = c.canvas.box[2]
    elseif legend.location in (:outer_bottom_left, :outer_bottom_right)
        y1 = c.canvas.box[4] - legend.height
    end
    x2, y2 = x1 + legend.width, y1 + legend.height

    set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))

    # draw rounded rectangle
    r = 0.02*min(c.canvas.width, c.canvas.height)
    move_to(ctx, x1, y1+r)
    line_to(ctx, x1, y2-r)
    curve_to(ctx, x1, y2, x1, y2, x1+r, y2)
    line_to(ctx, x2-r, y2)
    curve_to(ctx, x2, y2, x2, y2, x2, y2-r)
    line_to(ctx, x2, y1+r)
    curve_to(ctx, x2, y1, x2, y1, x2-r, y1)
    line_to(ctx, x1+r, y1)
    curve_to(ctx, x1, y1, x1, y1, x1, y1+r)
    close_path(ctx)
    set_source_rgb(ctx, 1, 1, 1) # white
    fill_preserve(ctx)
    set_source_rgb(ctx, 0, 0, 0) # black
    set_line_width(ctx, 0.4)
    stroke(ctx)

    # draw labels
    label_heigh = maximum( getsize(plot.label, legend.font_size)[2] for plot in plots )

    for (k, plot) in enumerate(plots)
        i = ceil(Int, k/ncols)  # line
        j = k%ncols==0 ? ncols : k%ncols # column
        x2 = x1 + inner_pad + sum(col_witdhs[1:j-1]) + (j-1)*col_sep

        y2 = y1 + inner_pad + label_heigh/2 + (i-1)*(label_heigh+row_sep)

        set_source_rgb(ctx, rgb(plot.color)...)
        if plot.line_style!=:none
            move_to(ctx, x2, y2)
            rel_line_to(ctx, handle_length, 0)
            set_line_width(ctx, plot.line_width)
            plot.line_style!=:solid && set_dash(ctx, plot.dash)
            stroke(ctx)
            set_dash(ctx, Float64[])
        end

        # draw marker
        x = x2 + handle_length/2
        draw_marker(ctx, x, y2, plot.marker, plot.marker_size, plot.marker_color, plot.marker_stroke_color)

        # draw label
        x = x2 + handle_length + 2*inner_pad
        y = y2 + 0.25*legend.font_size

        set_source_rgb(ctx, 0, 0, 0)
        draw_text(ctx, x, y, plot.label, halign="left", valign="bottom", angle=0)
    end

end


function draw!(c::Chart, ctx::CairoContext)
    # draw canvas grid
    draw!(c, ctx, c.canvas)

    # draw axes
    x = c.leftpad+c.yaxis.width
    y = c.toppad+c.yaxis.height
    move_to(ctx, x, y)
    draw!(ctx, c.xaxis)

    x = c.leftpad
    y = c.toppad
    move_to(ctx, x, y)
    draw!(ctx, c.yaxis)

    # draw plots
    x, y = c.canvas.box[1:2]
    w, h = c.canvas.box[3:4] - c.canvas.box[1:2]
    rectangle(ctx, x, y, w, h)
    Cairo.clip(ctx)

    # draw dataseries
    sorted = sort(c.dataseries, by=x->x.order)
    for p in sorted
        draw!(c, ctx, p)
    end
    reset_clip(ctx)

    # draw annotations
    for a in c.annotations
        draw!(c, ctx, a)
    end

    # draw legend
    if any( map(s->s.label!="", c.dataseries) )
        draw!(c, ctx, c.legend)
    end
end
