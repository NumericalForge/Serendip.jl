
mutable struct Annotation <: FigureComponent
    text       ::AbstractString
    x          ::Float64
    y          ::Float64
    alignment  ::Symbol
    target     ::Vector{Float64}
    has_target ::Bool
    line_width ::Float64
    font       ::String
    font_size  ::Float64
    color      ::Symbol
    function Annotation(text::AbstractString, x::Real, y::Real;
        alignment::Symbol=:auto,
        target::Union{Nothing,AbstractArray{<:Real,1}}=[0.0, 0.0],
        line_width::Real=0.4,
        font::AbstractString="NewComputerModern",
        font_size::Real=6.0,
        color::Symbol=:black
    )

        @check 0 <= x <= 1 "x must be in the range [0,1]"
        @check 0 <= y <= 1 "y must be in the range [0,1]"
        @check alignment in (:auto, :left, :right, :top, :bottom) "Invalid alignment: $(repr(alignment))"
        @check line_width > 0 "line_width must be positive"
        @check font_size > 0 "font_size must be positive"

        has_target = target !== nothing
        if has_target
            @check length(target) == 2 "target must be a 2D point"
            target = float.(target)
        else
            target = [0.0, 0.0]
        end

        return new(text, x, y, alignment, target, has_target, line_width, font, font_size, color)
    end
end


function add_annotation(c::Figure, a::Annotation)
    push!(c.annotations, a)
end


function draw!(c::Figure, cc::CairoContext, a::Annotation)

    set_font_size(cc, a.font_size)
    font = get_font(a.font)
    select_font_face(cc, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)

    set_matrix(cc, CairoMatrix([1, 0, 0, 1, 0, 0]...))

    # convert from axes to Cairo coordinates
    x = c.canvas.box[1] + a.x * c.canvas.width
    y = c.canvas.box[2] + (1 - a.y) * c.canvas.height
    halign = a.alignment == :right ? "right" : "left"
    valign = a.alignment == :top ? "top" : "bottom"
    set_source_rgb(cc, 0, 0, 0)
    draw_text(cc, x, y, a.text, halign=halign, valign=valign, angle=0)

    if a.alignment == :auto
        a.alignment = :left
    end

    # draw arrow
    if a.has_target

        # compute text size
        w, h = getsize(cc, a.text, a.font_size)
        text_outerpad = 0.1 * min(w, h)

        x = halign == "left" ? x + w/2 : x - w/2
        y = valign == "top" ? y + h/2 : y - h/2

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
        x_prev, y_prev = x, y
        move_to(cc, x, y)
        if lines[2] == '|'
            rel_line_to(cc, 0, dy)
            y += dy
        else
            rel_line_to(cc, dx, 0)
            x += dx
        end

        stroke(cc)

        # Draw a concave (chevron/notched) arrowhead at target.
        vx = xa - x_prev
        vy = ya - y_prev
        vlen = hypot(vx, vy)
        if vlen < 1e-8
            vx = xa - x
            vy = ya - y
            vlen = hypot(vx, vy)
        end
        if vlen > 1e-8
            ux = vx / vlen
            uy = vy / vlen
            nx = -uy
            ny = ux

            head_len = max(6.0 * a.line_width, 4.0)
            head_w = 0.9 * head_len
            notch_depth = 0.45 * head_len

            tipx, tipy = xa, ya
            basex = tipx - head_len * ux
            basey = tipy - head_len * uy
            leftx = basex + 0.5 * head_w * nx
            lefty = basey + 0.5 * head_w * ny
            rightx = basex - 0.5 * head_w * nx
            righty = basey - 0.5 * head_w * ny
            notchx = tipx - (head_len - notch_depth) * ux
            notchy = tipy - (head_len - notch_depth) * uy

            move_to(cc, tipx, tipy)
            line_to(cc, leftx, lefty)
            line_to(cc, notchx, notchy)
            line_to(cc, rightx, righty)
            close_path(cc)
            fill(cc)
        end
    end
end
