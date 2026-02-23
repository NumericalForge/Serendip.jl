# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

const _line_style_list = [:none, :solid, :dot, :dash, :dashdot]
const _mark_list = [:none, :circle, :square, :triangle, :utriangle, :cross, :xcross, :diamond, :pentagon, :hexagon, :star]


mutable struct DataSeries
    kind             ::Symbol  # type of data series, e.g., :line, :scatter, :bar, etc.
    X                ::AbstractArray
    Y                ::AbstractArray
    line_style       ::Symbol
    line_width       ::Float64
    color            ::Union{Symbol,Color}
    dash             ::Vector{Float64}
    mark             ::Symbol
    mark_size        ::Float64
    mark_color       ::Union{Symbol,Color}
    mark_stroke_color::Union{Symbol,Color}
    label            ::AbstractString
    tag              ::AbstractString
    tag_location     ::Symbol
    tag_position     ::Float64
    tag_alignment    ::Symbol
    bar_width        ::Float64
    bar_base         ::Float64
    order            ::Int

    function DataSeries(kind::Symbol, X::AbstractArray, Y::AbstractArray;
                            line_style=:solid,
                            line_width=0.5,
                            color=:default,
                            dash=Float64[],
                            mark=:none, mark_size=2.5,
                            mark_color=:white, mark_stroke_color=:default,
                            label="", tag="", tag_location=:top, tag_position=0.5,
                            tag_alignment=:horizontal,
                            bar_width=0.0,
                            bar_base=0.0,
                            order=0
    )
        @check kind in [:line, :scatter, :bar] "Invalid data series kind: $(repr(kind))"
        @check length(X)==length(Y) "Length of X and Y arrays must be equal"
        @check line_style in _line_style_list "Invalid line style: $(repr(line_style))"
        @check mark in _mark_list "Invalid mark shape: $(repr(mark))"
        @check tag_location in [:bottom, :top, :left, :right] "Invalid tag location: $(repr(tag_location))"
        @check tag_alignment in [:horizontal, :vertical, :parallel] "Invalid tag alignment: $(repr(tag_alignment))"
        @check line_width > 0 "Line width must be greater than zero"
        @check mark_size > 0 "Mark size must be greater than zero"
        @check bar_width >= 0 "Bar width must be non-negative"


        if color != :default && color isa Symbol
            color = Color(color)
        elseif color isa Tuple
            color = Color(color)
        end

        if mark_color != :default && mark_color isa Symbol
            mark_color = Color(mark_color)
        elseif mark_color isa Tuple
            mark_color = Color(mark_color)
        end

        if mark_stroke_color != :default && mark_stroke_color isa Symbol
            mark_stroke_color = Color(mark_stroke_color)
        elseif mark_stroke_color isa Tuple
            mark_stroke_color = Color(mark_stroke_color)
        end

        n = min(length(X), length(Y))

        if length(dash)==0
            if line_style==:dash
                dash = [4.0, 2.4]*line_width
            elseif line_style==:dashdot
                dash = [2.0, 1.0, 2.0, 1.0]*line_width
            elseif line_style==:dot
                dash = [1.0, 1.0]*line_width
            end
        else
            line_style = :dash
        end

        return new(kind, X[1:n], Y[1:n],
                    line_style, line_width, color,
                    dash,
                    mark, mark_size, mark_color, mark_stroke_color,
                    label, tag, tag_location, tag_position, tag_alignment,
                    bar_width, bar_base,
                    order)
    end
end


function DataSeries(X::AbstractArray, Y::AbstractArray; args...)
    return DataSeries(:line, X, Y; args...)
end


function draw_polygon(cc::CairoContext, x, y, n, length, color, strokecolor; angle=0)
    Δθ = 360/n
    minθ = angle + 90
    maxθ = angle + 360 + 90

    for θ in minθ:Δθ:maxθ
        xi = x + length*cosd(θ)
        yi = y - length*sind(θ)
        if θ==angle
            move_to(cc, xi, yi)
        else
            line_to(cc, xi, yi)
        end
    end

    close_path(cc)
    set_source_rgb(cc, rgb(color)...)
    fill_preserve(cc)
    set_source_rgb(cc, rgb(strokecolor)...)
    stroke(cc)
end


function draw_star(cc::CairoContext, x, y, n, length, color, strokecolor; angle=0)
    Δθ = 360/n/2
    minθ = angle + 90
    maxθ = angle + 360 + 90


    for (i,θ) in enumerate(minθ:Δθ:maxθ)
        if i%2==1
            xi = x + length*cosd(θ)
            yi = y - length*sind(θ)
        else
            xi = x + 0.5*length*cosd(θ)
            yi = y - 0.5*length*sind(θ)
        end
        if θ==angle
            move_to(cc, xi, yi)
        else
            line_to(cc, xi, yi)
        end
    end

    close_path(cc)
    set_source_rgb(cc, rgb(color)...)
    fill_preserve(cc)
    set_source_rgb(cc, rgb(strokecolor)...)
    stroke(cc)
end


function draw_mark(cc::CairoContext, x, y, mark, size, color, strokecolor)
    radius = size/2
    new_path(cc)

    if mark==:circle
        arc(cc, x, y, radius, 0, 2*pi)
        set_source_rgb(cc, rgb(color)...)
        fill_preserve(cc)
        set_source_rgb(cc, rgb(strokecolor)...)
        stroke(cc)
    elseif mark==:square
        draw_polygon(cc, x, y, 4, 1.2*radius, color, strokecolor, angle=45)
    elseif mark==:diamond
        draw_polygon(cc, x, y, 4, 1.2*radius, color, strokecolor, angle=0)
    elseif mark==:triangle
        draw_polygon(cc, x, y, 3, 1.3*radius, color, strokecolor, angle=0)
    elseif mark==:utriangle
        draw_polygon(cc, x, y, 3, 1.3*radius, color, strokecolor, angle=180)
    elseif mark==:pentagon
        draw_polygon(cc, x, y, 5, 1.1*radius, color, strokecolor, angle=0)
    elseif mark==:hexagon
        draw_polygon(cc, x, y, 6, 1.1*radius, color, strokecolor, angle=0)
    elseif mark==:star
        draw_star(cc, x, y, 5, 1.25*radius, color, strokecolor, angle=0)
    elseif mark==:cross
        radius = 1.35*radius
        set_line_width(cc, radius/3)
        move_to(cc, x, y-radius)
        line_to(cc, x, y+radius)
        stroke(cc)
        move_to(cc, x-radius, y)
        line_to(cc, x+radius, y)
        stroke(cc)
    elseif mark==:xcross
        radius = 1.35*radius
        set_line_width(cc, radius/3)
        move_to(cc, x-radius, y-radius)
        line_to(cc, x+radius, y+radius)
        stroke(cc)
        move_to(cc, x+radius, y-radius)
        line_to(cc, x-radius, y+radius)
        stroke(cc)
    end
end
