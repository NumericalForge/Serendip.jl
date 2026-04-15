# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

_axes_widget_locations = [:none, :top_right, :top_left, :bottom_right, :bottom_left]

function axis_widget_label(label::AbstractString)
    (occursin('$', label) || occursin('`', label)) && return label
    return "\$" * label * "\$"
end

axis_widget_head_length(aw) = 7.0

"""
    AxesWidget(;
        location=:bottom_left,
        labels=["x", "y"],
        font="NewComputerModern",
        font_size=9.0,
        azimuth=0.0,
        elevation=0.0,
        distance=0.0,
        up=:z,
        arrow_length=20.0,
    )

Create a small orientation triad overlay for plot figures.

# Keywords
- `location::Symbol`: widget anchor location (`:none`, `:top_right`, `:top_left`, `:bottom_right`, `:bottom_left`).
- `labels::Vector{<:AbstractString}`: axis labels; use length 2 for a planar widget or length 3 for a 3D triad.
- `font::AbstractString`: label font family.
- `font_size::Real`: label font size.
- `azimuth::Real`: 3D azimuth angle in degrees.
- `elevation::Real`: 3D elevation angle in degrees.
- `distance::Real`: camera distance used for 3D projection.
- `up::Symbol`: 3D up direction (`:x | :y | :z`).
- `arrow_length::Real`: target side length of the widget box in screen-space units.

# Notes
- For two labels, the widget draws the original 2D right/up arrows.
- For three labels, the triad is projected with the same 3D camera convention used by `DomainPlot`.
"""
mutable struct AxesWidget<:FigureComponent
    location::Symbol
    labels::Vector{AbstractString}
    font::String
    font_size::Float64
    azimuth::Float64
    elevation::Float64
    distance::Float64
    up::Symbol
    arrow_length::Float64
    projected_origin::Vec3
    projected_axes::Vector{Vec3}
    projected_limits::Vector{Float64}
    width::Float64
    height::Float64
    frame::Frame

    function AxesWidget(;
        location::Symbol = :bottom_left,
        labels::Vector{<:AbstractString} = ["`x`", "`y`"],
        font::String = "NewComputerModern",
        font_size::Real = 9.0,
        azimuth::Real = 0.0,
        elevation::Real = 0.0,
        distance::Real = 0.0,
        up::Symbol = :z,
        arrow_length::Real = 20.0,
    )

        @check up in (:x, :y, :z) "up must be one of :x, :y, :z. Got $up"

        this = new(location, labels, font, font_size, azimuth, elevation, distance, up, arrow_length, Vec3(0.0, 0.0, 0.0), Vec3[], Float64[], 0.0, 0.0, Frame())
        return this
    end
end


function configure!(aw::AxesWidget)
    ndim = length(aw.labels)
    head_width = 0.15*aw.arrow_length
    pad = head_width + 0.09*aw.arrow_length

    if ndim == 2
        aw.projected_origin = Vec3(0.0, 0.0, 0.0)
        aw.projected_axes = Vec3[
            Vec3(aw.arrow_length, 0.0, 0.0),
            Vec3(0.0, aw.arrow_length, 0.0),
        ]
    elseif ndim == 3
        s = aw.arrow_length
        axis_points = Vec3[
            Vec3(0.0, 0.0, 0.0),
            Vec3(s, 0.0, 0.0),
            Vec3(0.0, s, 0.0),
            Vec3(0.0, 0.0, s),
        ]
        cube_points = Vec3[
            Vec3(i*s, j*s, k*s)
            for i in (0.0, 1.0), j in (0.0, 1.0), k in (0.0, 1.0)
        ]

        projected_axes = [
            project_view_point(point, aw.azimuth, aw.elevation, aw.distance, up=aw.up, center=Vec3(0.0, 0.0, 0.0), reflength=s)
            for point in axis_points
        ]
        projected_cube = [
            project_view_point(point, aw.azimuth, aw.elevation, aw.distance, up=aw.up, center=Vec3(0.0, 0.0, 0.0), reflength=s)
            for point in cube_points
        ]

        cube_xmin = minimum(point[1] for point in projected_cube)
        cube_xmax = maximum(point[1] for point in projected_cube)
        cube_ymin = minimum(point[2] for point in projected_cube)
        cube_ymax = maximum(point[2] for point in projected_cube)
        cube_side = max(cube_xmax - cube_xmin, cube_ymax - cube_ymin)
        scale = cube_side > 1e-8 ? aw.arrow_length / cube_side : 1.0

        origin = projected_axes[1]
        axes = projected_axes[2:end]
        aw.projected_origin = Vec3((origin[1] - cube_xmin) * scale, (origin[2] - cube_ymin) * scale, origin[3])
        aw.projected_axes = Vec3[
            Vec3((point[1] - cube_xmin) * scale, (point[2] - cube_ymin) * scale, point[3])
            for point in axes
        ]
    else
        error("AxesWidget: labels must have length 2 or 3")
    end

    xs = Float64[aw.projected_origin[1]]
    ys = Float64[aw.projected_origin[2]]
    label_pad = 0.09*aw.arrow_length
    for point in aw.projected_axes
        push!(xs, point[1])
        push!(ys, point[2])

        dx = point[1] - aw.projected_origin[1]
        dy = point[2] - aw.projected_origin[2]
        ℓ = hypot(dx, dy)
        if ℓ > 1e-8
            push!(xs, point[1] + label_pad*dx/ℓ)
            push!(ys, point[2] + label_pad*dy/ℓ)
        end
    end

    xmin = minimum(xs)
    xmax = maximum(xs)
    ymin = minimum(ys)
    ymax = maximum(ys)

    aw.projected_limits = [xmin, ymin, xmax, ymax]
    aw.width = max(xmax - xmin + 2*pad, pad)
    aw.height = max(ymax - ymin + 2*pad, pad)
end


function draw_arrow(ctx::CairoContext, x1, y1, x2, y2; head_length=7)
    Δx = x2 - x1
    Δy = y2 - y1
    ℓ = hypot(Δx, Δy)

    # Draw the shaft
    move_to(ctx, x1, y1)
    line_to(ctx, x2, y2)
    stroke(ctx)

    ℓ < 1e-8 && return
    ℓ < 0.8*head_length && return

    θ = π/8

    N  = Vec3(0, 0, 1)
    v1 = normalize(Vec3(Δx, Δy, 0))
    v3 = v1*cos(π+θ) + cross(N, v1)*sin(π+θ)
    v4 = v1*cos(π-θ) + cross(N, v1)*sin(π-θ)
    p2 = Vec3(x2, y2, 0)
    p0 = p2 - v1*head_length*0.7
    p3 = p2 + v3*head_length/cos(θ/2)
    p4 = p2 + v4*head_length/cos(θ/2)


    move_to(ctx, x2, y2)
    line_to(ctx, p3[1], p3[2])
    line_to(ctx, p0[1], p0[2])
    line_to(ctx, p4[1], p4[2])

    close_path(ctx)
    fill(ctx)

end


function draw!(aw::AxesWidget, ctx::RenderContext)
    cairo_ctx = ctx.cairo_ctx

    Cairo.save(cairo_ctx)
    x0 = aw.frame.x
    y0 = aw.frame.y

    font = get_font(aw.font)
    select_font_face(cairo_ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
    set_font_size(cairo_ctx, aw.font_size)
    ndim = length(aw.labels)
    line_width = ndim == 3 ? 1.0 * ctx.width_scale : 0.7 * ctx.width_scale
    set_line_width(cairo_ctx, line_width)

    reset_matrix!(ctx)
    pad = 0.15*aw.arrow_length + 0.09*aw.arrow_length
    xmin, ymin = aw.projected_limits[1], aw.projected_limits[2]

    function to_screen(point::Vec3)
        x = x0 + pad + (point[1] - xmin)
        y = y0 + aw.height - pad - (point[2] - ymin)
        return x, y
    end

    colors = [_colors_dict[:indianred], _colors_dict[:green], _colors_dict[:blue]]
    order = length(aw.labels) == 3 ? sortperm([point[3] for point in aw.projected_axes], rev=true) : collect(eachindex(aw.projected_axes))
    origin_x, origin_y = to_screen(aw.projected_origin)
    label_pad = 0.09*aw.arrow_length
    head_length = axis_widget_head_length(aw)

    for idx in order
        endpoint = aw.projected_axes[idx]
        x2, y2 = to_screen(endpoint)
        dx = x2 - origin_x
        dy = y2 - origin_y
        ℓ = hypot(dx, dy)

        set_source_rgb(cairo_ctx, colors[idx]...)
        draw_arrow(cairo_ctx, origin_x, origin_y, x2, y2; head_length=head_length)

        dirx = ℓ > 1e-8 ? dx/ℓ : 1.0
        diry = ℓ > 1e-8 ? dy/ℓ : 0.0
        label_x = x2 + label_pad*dirx
        label_y = y2 + label_pad*diry

        halign = dirx > 0.2 ? "left" : dirx < -0.2 ? "right" : "center"
        valign = diry > 0.2 ? "top" : diry < -0.2 ? "bottom" : "center"
        draw_text(cairo_ctx, label_x, label_y, axis_widget_label(String(aw.labels[idx])), halign=halign, valign=valign, angle=0)
    end

    Cairo.restore(cairo_ctx)
end
