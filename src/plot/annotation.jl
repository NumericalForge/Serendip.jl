
# Annotation_params = [
#     FunInfo(:Annotation, "Creates an `Annotation` instance."),
#     ArgInfo(:text, "Text to be displayed", type=AbstractString),
#     ArgInfo(:x, "x-coordinate of the annotation", cond=:(0<=x<=1), type=Real),
#     ArgInfo(:y, "y-coordinate of the annotation", cond=:(0<=y<=1), type=Real),
#     KwArgInfo(:text_alignment, "Alignment of the text", :auto, values=(:auto, :left, :right, :top, :bottom)),
#     KwArgInfo(:target, "Coordinates of the target point of the arrow relative to data", [0.0,0.0], length=2),
#     KwArgInfo(:line_width, "Edge weight", 0.4, cond=:(line_width>0)),
#     KwArgInfo(:font, "Name of the font", "NewComputerModern", type=AbstractString),
#     KwArgInfo(:fontsize, "Size of the font in dpi", 6.0, cond=:(fontsize>0)),
#     KwArgInfo(:color, "Color of the text", :default),
# ]


mutable struct Annotation <: FigureComponent
    text::AbstractString
    x::Float64
    y::Float64
    text_alignment::Symbol
    target::Vector{Float64}
    line_width::Float64
    font::String
    fontsize::Float64
    color::Symbol
    function Annotation(text::AbstractString, x::Real, y::Real;
        text_alignment::Symbol=:auto,
        target::AbstractArray{<:Real,1}=[0.0, 0.0],
        line_width::Real=0.4,
        font::AbstractString="NewComputerModern",
        fontsize::Real=6.0,
        color::Symbol=:black
    )

        @check 0 <= x <= 1 "x must be in the range [0,1]"
        @check 0 <= y <= 1 "y must be in the range [0,1]"
        @check text_alignment in (:auto, :left, :right, :top, :bottom) "Invalid text_alignment: $(repr(text_alignment))"
        @check length(target) == 2 "target must be a 2D point"
        @check line_width > 0 "line_width must be positive"
        @check fontsize > 0 "fontsize must be positive"

        # args = checkargs([text, x, y], kwargs, Annotation_params)
        target = float.(target)

        return new(text, x, y, text_alignment, target, line_width, font, fontsize, color)
    end
end


function add_annotation(c::Figure, a::Annotation)
    push!(c.annotations, a)
end


function draw!(c::Figure, cc::CairoContext, a::Annotation)

    set_font_size(cc, a.fontsize)
    font = get_font(a.font)
    select_font_face(cc, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)

    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))

    # convert from axes to Cairo coordinates
    x = c.canvas.box[1] + a.x * c.canvas.width
    y = c.canvas.box[2] + (1 - a.y) * c.canvas.height
    halign = a.text_alignment == :right ? "right" : "left"
    valign = a.text_alignment == :top ? "top" : "bottom"
    set_source_rgb(cc, 0, 0, 0)
    draw_text(cc, x, y, a.text, halign=halign, valign=valign, angle=0)

    if a.text_alignment == :auto
        a.text_alignment = :left
    end

    # draw arrow
    if a.target !== nothing

        # compute text size
        w, h = getsize(cc, a.text, a.fontsize)
        text_outerpad = 0.1 * min(w, h)

        if halign == "left"
            x += w / 2
        else
            x -= w / 2
        end

        if valign == "top"
            y += 0.5 * h
        else
            y -= 0.5 * h
        end

        w += text_outerpad
        h += text_outerpad

        # target coordinates
        xa, ya = data2user(c.canvas, a.target[1], a.target[2])

        # deltas
        dx = xa - x
        dy = ya - y

        # compute lines
        if abs(dx) > abs(dy)
            if abs(dy) < h / 2
                lines = "-|"
                if dx > 0
                    x += w / 2 # right
                else
                    x -= w / 2 # left
                end
            else # two lines
                lines = "|-"
                if dy > 0 # top
                    y += h / 2
                else # bottom
                    y -= h / 2
                end
            end
        else
            if abs(dx) < w / 2
                lines = "|-"
                if dy > 0
                    y += h / 2 # top
                else
                    y -= h / 2 # bottom
                end
            else # two lines
                lines = "-|"
                if dx > 0 # right
                    x += w / 2
                else # left
                    x -= w / 2
                end
            end
        end

        set_source_rgb(cc, _colors_dict[a.color]...)

        set_line_join(cc, Cairo.CAIRO_LINE_JOIN_ROUND)
        set_line_width(cc, a.line_width)

        # update deltas
        dx = xa - x
        dy = ya - y

        # Draw line 1
        move_to(cc, x, y)
        if lines[1] == '|'
            rel_line_to(cc, 0, dy)
            y += dy
        else
            rel_line_to(cc, dx, 0)
            x += dx
        end

        # Draw line 2
        dx += sign(dx) * a.line_width
        move_to(cc, x, y)
        if lines[2] == '|'
            rel_line_to(cc, 0, dy)
            y += dy
        else
            rel_line_to(cc, dx, 0)
            x += dx
        end

        stroke(cc)
    end
end