# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct Colorbar<:FigureComponent
    location::Symbol
    colormap::Colormap
    axis::Union{Axis,Nothing}
    length_factor::Float64
    inner_sep::Float64
    box::Vector{Float64}
    width::Float64
    height::Float64
    thickness::Float64  # only bar thickness

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
        @check location in (:none, :right, :bottom) "Colorbar location must be :none, :right or :bottom"
        @check length_factor>0 "Colorbar length_factor must be >0"
        @check font_size>0 "Colorbar font_size must be >0"
        @check length(limits)==2 "Colorbar limits must be a vector of length 2"
    
        colormap = colormap isa Symbol ? Colormap(colormap) : colormap

        axis = nothing
        if location != :none
            direction = location == :right ? :vertical : :horizontal
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

        return new(location, colormap, axis, length_factor, inner_sep, Float64[0,0,0,0])

    end
end


function configure!(fig::Figure, cb::Colorbar)
    if cb.location!==:none
        configure!(cb.axis)
        # cb.thickness = 0.035*max(fig.width, fig.height)
        cb.thickness = 1.33*cb.axis.font_size
        cb.inner_sep = cb.thickness
        if cb.location==:right
            cb.height = cb.length_factor*(fig.height - 2*fig.outerpad)
            cb.axis.height = cb.height
            cb.width = cb.inner_sep + cb.thickness + cb.axis.tick_length + cb.axis.width
        elseif cb.location==:bottom
            cb.width = cb.length_factor*(fig.width - 2*fig.outerpad)
            cb.axis.width = cb.width
            cb.height = cb.inner_sep + cb.thickness + cb.axis.tick_length + cb.axis.height
        end
    end
end


function draw!(fig::Figure, cc::CairoContext, cb::Colorbar)
    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))

    if cb.location==:right
        # Axis
        fmin, fmax = cb.axis.limits
        x = fig.canvas.box[3]+cb.inner_sep+cb.thickness+cb.axis.tick_length
        h = cb.height

        y = fig.height/2 - h/2
        move_to(cc, x, y)
        draw!(cc, cb.axis)

        # Colorbar
        x = fig.canvas.box[3] + cb.inner_sep
        w = cb.thickness
        y = fig.height/2 + h/2

        pat = pattern_create_linear(0.0, y,  0.0, y-h)
        nstops = length(cb.colormap.stops)
        for i in 1:nstops
            stop = cb.colormap.stops[i]
            color = cb.colormap.colors[i]
            stop = round((stop-fmin)/(fmax-fmin), digits=8)
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        set_source(cc, pat)
        rectangle(cc, x, y, w, -h)
        fill(cc)
    elseif cb.location==:bottom
        # Axis
        fmin, fmax = cb.axis.limits
        w = cb.width
        x = fig.width/2 - w/2
        y = fig.canvas.box[4]+cb.inner_sep+cb.thickness+cb.axis.tick_length
        move_to(cc, x, y)
        draw!(cc, cb.axis)

        # Colorbar
        x = fig.width/2 - w/2
        y = fig.canvas.box[4] + cb.inner_sep
        h = cb.thickness

        pat = pattern_create_linear(x, 0.0,  x+w, 0.0)
        nstops = length(cb.colormap.stops)
        for i in 1:nstops
            stop = cb.colormap.stops[i]
            color = cb.colormap.colors[i]
            stop = round((stop-fmin)/(fmax-fmin), digits=8)
            pattern_add_color_stop_rgb(pat, stop, color...)
        end

        set_source(cc, pat)
        rectangle(cc, x, y, w, h)
        fill(cc)
    end
end
