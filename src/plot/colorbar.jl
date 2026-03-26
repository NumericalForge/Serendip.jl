# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct Colorbar<:FigureComponent
    location::Symbol
    colormap::Colormap
    axis::Union{Axis,Nothing}
    length_factor::Float64
    inner_sep::Float64
    width::Float64
    height::Float64
    thickness::Float64  # only bar thickness
    frame::Frame

    function Colorbar(;
        location::Symbol = :right,
        colormap = :coolwarm,
        limits = [0.0, 0.0],
        label = "",
        font_size = 9.0,
        font = "NewComputerModern",
        ticks = Float64[],
        tick_labels = String[],
        tick_length = 3,
        bins = 6,
        inner_sep = 3,
        length_factor = 1.0,
    )
        @check location in (:none, :left, :right, :top, :bottom) "Colorbar location must be :none, :left, :right, :top or :bottom"
        @check length_factor>0 "Colorbar length_factor must be >0"
        @check font_size>0 "Colorbar font_size must be >0"
        @check length(limits)==2 "Colorbar limits must be a vector of length 2"
    
        colormap = colormap isa Symbol ? Colormap(colormap) : colormap

        axis = nothing
        if location != :none
            direction = location in (:left, :right) ? :vertical : :horizontal
            axis =  Axis(;
                direction   = direction,
                location    = location,
                limits      = limits,
                label       = label,
                font_size   = font_size,
                font        = font,
                ticks       = ticks,
                tick_labels = tick_labels,
                tick_length = tick_length,
                bins        = bins,
            )
        end

        return new(location, colormap, axis, length_factor, inner_sep, 0.0, 0.0, 0.0, Frame())

    end
end


function configure!(fig::Figure, cb::Colorbar)
    if cb.location!==:none
        configure!(cb.axis)
        cb.thickness = 1.33*cb.axis.font_size
        cb.inner_sep = cb.thickness
        if cb.location in (:left, :right)
            cb.height = cb.length_factor*(fig.height - 2*fig.outerpad)
            cb.axis.height = cb.height
            cb.width = cb.inner_sep + cb.thickness + cb.axis.tick_length + cb.axis.width
        elseif cb.location in (:top, :bottom)
            cb.width = cb.length_factor*(fig.width - 2*fig.outerpad)
            cb.axis.width = cb.width
            cb.height = cb.inner_sep + cb.thickness + cb.axis.tick_length + cb.axis.height
        end
    end
end


function draw!(fig::Figure, ctx::RenderContext, cb::Colorbar)
    cairo_ctx = ctx.cairo_ctx
    reset_matrix!(ctx)
    cb.location == :none && return

    x0 = cb.frame.x
    y0 = cb.frame.y
    x1 = cb.frame.x + cb.frame.width
    y1 = cb.frame.y + cb.frame.height
    fmin, fmax = cb.axis.limits

    if cb.location == :right
        bar_x = x0 + cb.inner_sep
        bar_y = y0 + (cb.frame.height - cb.height) / 2
        cb.axis.frame = Frame(bar_x + cb.thickness + cb.axis.tick_length, bar_y, cb.axis.width, cb.height)

        draw!(cb.axis, ctx)

        pat = pattern_create_linear(0.0, bar_y + cb.height, 0.0, bar_y)
        nstops = length(cb.colormap.stops)
        for i in 1:nstops
            stop = cb.colormap.stops[i]
            color = cb.colormap.colors[i]
            stop = round((stop-fmin)/(fmax-fmin), digits=8)
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        set_source(cairo_ctx, pat)
        rectangle(cairo_ctx, bar_x, bar_y, cb.thickness, cb.height)
        fill(cairo_ctx)
    elseif cb.location == :left
        bar_x = x1 - cb.inner_sep - cb.thickness
        bar_y = y0 + (cb.frame.height - cb.height) / 2
        cb.axis.frame = Frame(x0, bar_y, cb.axis.width, cb.height)

        draw!(cb.axis, ctx)

        pat = pattern_create_linear(0.0, bar_y + cb.height, 0.0, bar_y)
        nstops = length(cb.colormap.stops)
        for i in 1:nstops
            stop = cb.colormap.stops[i]
            color = cb.colormap.colors[i]
            stop = round((stop-fmin)/(fmax-fmin), digits=8)
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        set_source(cairo_ctx, pat)
        rectangle(cairo_ctx, bar_x, bar_y, cb.thickness, cb.height)
        fill(cairo_ctx)
    elseif cb.location == :bottom
        bar_x = x0 + (cb.frame.width - cb.width) / 2
        bar_y = y0 + cb.inner_sep
        cb.axis.frame = Frame(bar_x, bar_y + cb.thickness + cb.axis.tick_length, cb.width, cb.axis.height)

        draw!(cb.axis, ctx)

        pat = pattern_create_linear(bar_x, 0.0, bar_x + cb.width, 0.0)
        nstops = length(cb.colormap.stops)
        for i in 1:nstops
            stop = cb.colormap.stops[i]
            color = cb.colormap.colors[i]
            stop = round((stop-fmin)/(fmax-fmin), digits=8)
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        set_source(cairo_ctx, pat)
        rectangle(cairo_ctx, bar_x, bar_y, cb.width, cb.thickness)
        fill(cairo_ctx)
    elseif cb.location == :top
        bar_x = x0 + (cb.frame.width - cb.width) / 2
        bar_y = y1 - cb.inner_sep - cb.thickness
        cb.axis.frame = Frame(bar_x, y0, cb.width, cb.axis.height)

        draw!(cb.axis, ctx)

        pat = pattern_create_linear(bar_x, 0.0, bar_x + cb.width, 0.0)
        nstops = length(cb.colormap.stops)
        for i in 1:nstops
            stop = cb.colormap.stops[i]
            color = cb.colormap.colors[i]
            stop = round((stop-fmin)/(fmax-fmin), digits=8)
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        set_source(cairo_ctx, pat)
        rectangle(cairo_ctx, bar_x, bar_y, cb.width, cb.thickness)
        fill(cairo_ctx)
    end
end
