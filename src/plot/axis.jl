# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

"""
    Axis(; kwargs...)

Creates an `Axis` component.

# Parameters
- `direction::Symbol` : Axis direction (`:horizontal` or `:vertical`, default `:horizontal`).
- `location::Symbol` : Axis location (`:none`, `:left`, `:right`, `:top`, `:bottom`, default `:none`).
  If `:none`, it becomes `:bottom` for horizontal and `:left` for vertical.
- `limits::Vector{Float64}` : Axis limits `[min, max]`. Use `Float64[]` for auto limits.
- `label::AbstractString` : Axis label text (default `""`).
- `font::AbstractString` : Font family for label and ticks (default `"NewComputerModern"`).
- `font_size::Float64` : Font size (> 0, default `7.0`).
- `ticks::AbstractArray` : Explicit tick positions (default `Float64[]`).
- `tick_labels::AbstractArray` : Labels for ticks (default `String[]`). Must match `length(ticks)` when provided.
- `tick_length::Float64` : Tick mark length input (default `0.4*font_size`).
- `bins::Int` : Hint for the number of tick bins (default `6`).
- `mult::Float64` : Multiplier applied to axis values (default `1.0`).

# Notes
- `direction` and `location` control orientation and side.
- If `tick_labels` are given, their count must equal `length(ticks)`.
- `bins` is a hint. Explicit `ticks` take precedence.

"""
mutable struct Axis<:FigureComponent
    direction  ::Symbol
    location   ::Symbol
    limits     ::Vector{Float64}
    auto_limits::Bool
    label      ::AbstractString
    font       ::String
    font_size  ::Float64
    ticks      ::AbstractArray
    tick_labels::AbstractArray
    manual_ticks::Bool
    manual_tick_labels::Bool
    tick_length::Float64
    nbins      ::Int
    mult       ::Float64
    tick_exponent::Int
    exponent_box::TextBox
    exponent_width::Float64
    exponent_height::Float64
    inner_sep  ::Float64
    width      ::Float64
    height     ::Float64
    frame      ::Frame

    function Axis(;
        direction::Symbol=:horizontal,
        location::Symbol=:none,
        limits::AbstractVector{<:Real}=Float64[],
        label::AbstractString="",
        font::String="NewComputerModern",
        font_size::Real=7.0,
        ticks::AbstractArray=Real[],
        tick_labels::AbstractArray=String[],
        tick_length::Real=3.0,
        bins::Int=6,
        mult::Real=1.0,
        )

        if location==:none
            if direction==:horizontal
                location = :bottom
            else
                location = :left
            end
        end

        if length(tick_labels)>0
            @check length(ticks)==length(tick_labels) "Axis: the length of labels must match the number of ticks"
        end

        @check length(limits) in (0, 2) "Axis: limits must have length 2 (manual) or be empty (auto)"
        auto_limits = length(limits) == 0
        limits = auto_limits ? [0.0, 0.0] : collect(float.(limits))
        manual_ticks = length(ticks) > 0
        manual_tick_labels = length(tick_labels) > 0

        return new(direction, location, limits, auto_limits, label, font, font_size, ticks, tick_labels,
            manual_ticks, manual_tick_labels, tick_length, bins, mult, 0, TextBox(), 0.0, 0.0, 3.0, 0.0, 0.0, Frame())
    end
end


function get_bin_length(vinf, vsup, n)
    dx = (vsup-vinf)/n*sign(vsup-vinf)
    abs(dx) < eps(Float64) && return dx == 0 ? 1.0e-12 : dx
    ex = floor(log10(abs(dx))) # exponent
    mant = dx/10^ex

    if mant<=1.4
        mant=1
    elseif mant<=1.8 && vinf==0
        mant=1.5
    elseif mant<=2.25
        mant=2.0
    elseif mant<=3.0
        mant=2.5
    elseif mant<=4.5
        mant=4.0
    elseif mant<=7.5
        mant=5.0
    else
        mant=10.0
    end

    return round(mant*10^ex, sigdigits=2)*sign(vsup-vinf)
end


function _tick_round_digits(dv::Real)
    dv == 0 && return 12
    exponent = floor(Int, log10(abs(float(dv))))
    return max(0, -exponent + 2)
end


function compute_auto_limits(extent::AbstractVector{<:Real};
    padding_fraction::Real=0.03,
    degenerate_fraction::Real=0.05,
    degenerate_floor::Real=1.0e-6,
    mult::Real=1.0,
)
    @check length(extent) == 2 "Axis: extent must have length 2"
    @check padding_fraction >= 0 "Axis: padding_fraction must be non-negative"
    @check degenerate_fraction > 0 "Axis: degenerate_fraction must be positive"
    @check degenerate_floor > 0 "Axis: degenerate_floor must be positive"

    vmin = float(extent[1])
    vmax = float(extent[2])

    if !isfinite(vmin) || !isfinite(vmax)
        return [0.0, 1.0]
    end

    lower = min(vmin, vmax)
    upper = max(vmin, vmax)
    span = upper - lower

    if span == 0.0
        center = lower
        scale = max(abs(center), float(degenerate_floor))
        pad = degenerate_fraction * scale
        limits = [center - pad, center + pad]
    else
        pad = padding_fraction * span
        limits = [lower - pad, upper + pad]
    end

    return [limits[1] * mult, limits[2] * mult]
end


function _axis_common_exponent(values; threshold::Int=3)
    isempty(values) && return 0
    maxabs = maximum(abs, values)
    maxabs == 0 && return 0
    exponent = floor(Int, log10(maxabs))
    return abs(exponent) >= threshold ? exponent : 0
end


function _axis_digits(values; max_digits::Int=5)
    isempty(values) && return 0
    cleaned = [isapprox(v, 0.0; atol=1e-12) ? 0.0 : v for v in values]

    if all(v -> isapprox(v, round(v); atol=1e-8, rtol=1e-8), cleaned)
        return 0
    end

    for digits in 1:max_digits
        atol = 10.0^(-digits - 1)
        if all(v -> isapprox(v, round(v; digits=digits); atol=atol, rtol=1e-8), cleaned)
            return digits
        end
    end

    return max_digits
end


function _format_axis_value(value::Real, digits::Int)
    v = isapprox(value, 0.0; atol=1e-12) ? 0.0 : float(value)
    return Printf.format(Printf.Format("%." * string(digits) * "f"), v)
end


function make_ticklabels(ticks)
    values = collect(float.(ticks))
    exponent = _axis_common_exponent(values)
    scale = exponent == 0 ? 1.0 : 10.0^exponent
    mantissas = values ./ scale
    digits = _axis_digits(mantissas)
    labels = [_format_axis_value(m, digits) for m in mantissas]
    return labels, exponent
end


function _configure_exponent_box!(ax::Axis, exponent::Int)
    ax.tick_exponent = exponent
    if exponent == 0
        ax.exponent_box.text = ""
        ax.exponent_box.frame = Frame()
        ax.exponent_width = 0.0
        ax.exponent_height = 0.0
    else
        ax.exponent_box.text = "\$× 10^{$exponent}\$"
        ax.exponent_box.angle = 0.0
        ax.exponent_width, ax.exponent_height = getsize(ax.exponent_box.text, ax.font_size)
    end
end


function _axis_tick_label_extent(ax::Axis)
    isempty(ax.tick_labels) && return 0.0, 0.0
    widths = [getsize(lbl, ax.font_size)[1] for lbl in ax.tick_labels]
    heights = [getsize(lbl, ax.font_size)[2] for lbl in ax.tick_labels]
    return maximum(widths), maximum(heights)
end


function axis_right_overhang(ax::Axis)
    return 0.0
end


function axis_top_overhang(ax::Axis)
    if ax.direction == :vertical && !isempty(ax.exponent_box.text)
        return 0.3 * ax.font_size + ax.exponent_height
    end
    return 0.0
end


function configure!(ax::Axis)
    if ax.limits[1] == ax.limits[2]
        ax.limits = compute_auto_limits(ax.limits; mult=1.0)
    end

    # configure ticks
    if !ax.manual_ticks || length(ax.ticks) == 0
        vinf, vsup = ax.limits
        dv = get_bin_length(vinf, vsup, ax.nbins)
        digits = _tick_round_digits(dv)

        # roundup first tick
        m = mod(vinf,dv)
        vinf = m==0 ? vinf : vinf - m + dv

        ax.ticks = unique(round.(collect(vinf:dv:vsup), digits=digits))
    else # check ticks
        vmin = minimum(ax.limits)
        vmax = maximum(ax.limits)
        idxs = [ i for (i,tick) in enumerate(ax.ticks) if vmin<=tick<=vmax ]
        ax.ticks = ax.ticks[idxs]

        ax.manual_tick_labels && (ax.tick_labels = ax.tick_labels[idxs])
    end

    if ax.manual_tick_labels
        _configure_exponent_box!(ax, 0)
    else
        ax.tick_labels, exponent = make_ticklabels(ax.ticks)
        _configure_exponent_box!(ax, exponent)
    end

    ax.nbins = max(length(ax.ticks) - 1, 0)

    # configure size (axis dimensions do do not include tick lengths)
    ax.tick_length = 0.4*ax.font_size
    ax.inner_sep = 0.4*ax.font_size
    tk_lbs_width, tk_lbs_height = _axis_tick_label_extent(ax)
    label_height = getsize(ax.label, ax.font_size)[2]

    if ax.direction == :horizontal
        text_height = max(tk_lbs_height, ax.exponent_height)
        ax.height = label_height + ax.inner_sep + text_height + ax.tick_length
    else
        text_width = max(tk_lbs_width, ax.exponent_width)
        ax.width = label_height + ax.inner_sep + text_width + ax.inner_sep + ax.tick_length
    end

    return ax
end



function draw!(ax::Axis, ctx::RenderContext)
    cairo_ctx = ctx.cairo_ctx

    Cairo.save(cairo_ctx)
    x0 = ax.frame.x
    y0 = ax.frame.y

    font = get_font(ax.font)
    select_font_face(cairo_ctx, font, Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL )
    set_font_size(cairo_ctx, ax.font_size)
    reset_matrix!(ctx)

    set_source_rgb(cairo_ctx, 0, 0, 0) # black
    set_line_width(cairo_ctx, 0.05 * ax.font_size * ctx.width_scale)
    tk_lbs_width, tk_lbs_height = _axis_tick_label_extent(ax)

    if ax.direction==:horizontal
        label_height = getsize(ax.label, ax.font_size)[2]
        xmin, xmax = ax.limits
        baseline = ax.location == :top ? y0 + ax.frame.height : y0
        text_row_y = ax.location == :top ? baseline - ax.tick_length - tk_lbs_height / 2 : baseline + ax.tick_length + tk_lbs_height / 2

        # draw tick labels
        for (x,label) in zip(ax.ticks, ax.tick_labels)
            min(xmin, xmax) <= x <=max(xmin, xmax) || continue
            x1 = x0 + ax.frame.width/(xmax-xmin)*(x-xmin)

            move_to(cairo_ctx, x1, baseline)
            rel_line_to(cairo_ctx, 0, ax.location == :top ? ax.tick_length : -ax.tick_length)
            stroke(cairo_ctx)
            draw_text(cairo_ctx, x1, text_row_y, label, halign="center", valign="center")
        end

        x = x0 + ax.frame.width / 2
        y = ax.location == :top ? y0 + label_height / 2 : y0 + ax.frame.height - label_height / 2
        draw_text(cairo_ctx, x, y, ax.label, halign="center", valign="center", angle=0)

        if !isempty(ax.exponent_box.text)
            draw_text(cairo_ctx, x0 + ax.frame.width, y, ax.exponent_box.text, halign="right", valign="center")
        end

    else # :vertical ax
        label_height = getsize(ax.label, ax.font_size)[2]
        ymin, ymax = ax.limits

        if ax.location==:left
            halign="right"
            tick_length = ax.tick_length
            x1 = x0 + ax.frame.width
            xtext = x1 - tick_length
            xexp = x1 - tick_length - max(tk_lbs_width, ax.exponent_width) / 2
        else
            halign="left"
            tick_length = -ax.tick_length
            x1 = x0
            xtext = x1 - tick_length
            xexp = x1 - tick_length + max(tk_lbs_width, ax.exponent_width) / 2
        end

        # draw tick labels
        for (y,label) in zip(ax.ticks, ax.tick_labels)
            min(ymin, ymax) <= y <=max(ymin, ymax) || continue
            y1 = y0 + ax.frame.height/(ymax-ymin)*(ymax-y)

            move_to(cairo_ctx, x1, y1); rel_line_to(cairo_ctx, tick_length, 0); stroke(cairo_ctx)

            draw_text(cairo_ctx, xtext, y1+ax.font_size/2, label, halign=halign, valign="bottom")
        end

        if !isempty(ax.exponent_box.text)
            yexp = y0 - 0.3 * ax.font_size - ax.exponent_height / 2
            draw_text(cairo_ctx, xexp, yexp, ax.exponent_box.text, halign="center", valign="center")
        end

        # draw label
        if ax.location==:left
            x = x0 + ax.frame.width - ax.tick_length - tk_lbs_width - ax.inner_sep - label_height/2
        else
            x = x0 + ax.tick_length + tk_lbs_width + ax.inner_sep + label_height/2
        end
        y = y0 + ax.frame.height/2

        draw_text(cairo_ctx, x, y, ax.label, halign="center", valign="center", angle=90)
    end

    Cairo.restore(cairo_ctx)
end
