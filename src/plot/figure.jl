abstract type Figure
    # width::Float64
    # height::Float64
end

abstract type FigureComponent end
# abstract type DataSeries end

const cm = 72.0 / 2.54

mutable struct Frame
    x::Float64
    y::Float64
    width::Float64
    height::Float64

    function Frame(x::Real=0.0, y::Real=0.0, width::Real=0.0, height::Real=0.0)
        return new(float(x), float(y), float(width), float(height))
    end
end

mutable struct TextBox
    text::AbstractString
    frame::Frame
    angle::Float64
    visible::Bool

    function TextBox(text::AbstractString=""; frame::Frame=Frame(), angle::Real=0.0, visible::Bool=!isempty(text))
        return new(text, frame, float(angle), visible)
    end
end

xmin(frame::Frame) = frame.x
ymin(frame::Frame) = frame.y
xmax(frame::Frame) = frame.x + frame.width
ymax(frame::Frame) = frame.y + frame.height

const _available_formats = [
    ".pdf",
    ".png",
    ".svg",
    ".ps",
]

function _draw_figure_background!(ctx::CairoContext, frame::Frame, color)
    color === nothing && return

    set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
    rectangle(ctx, frame.x, frame.y, frame.width, frame.height)
    set_source_rgba(ctx, rgba(color)...)
    fill(ctx)
end

_has_text(box::TextBox) = box.visible && !isempty(box.text)

function _measure_text_box(box::TextBox, cc::CairoContext, fontsize::Float64)
    _has_text(box) || return 0.0, 0.0
    return getsize(cc, box.text, fontsize)
end

function _measure_text_box_extent(box::TextBox, cc::CairoContext, fontsize::Float64)
    width, height = _measure_text_box(box, cc, fontsize)
    angle = mod(box.angle, 180.0)
    if isapprox(angle, 90.0; atol=1e-8)
        return height, width
    end
    return width, height
end

function _draw_text_box!(ctx::CairoContext, box::TextBox)
    _has_text(box) || return
    x = box.frame.x + 0.5 * box.frame.width
    y = box.frame.y + 0.5 * box.frame.height
    draw_text(ctx, x, y, box.text, halign="center", valign="center", angle=box.angle)
end


function save(figure::Figure, files::String...)
    configure!(figure)

    for file in files
        width, height = figure.width, figure.height
        dir = dirname(file)
        ispath(dir) || mkpath(dir)

        fmt = splitext(file)[end]
        if fmt==".pdf"
            surf = CairoPDFSurface(file, width, height)
        elseif fmt==".svg"
            surf = CairoSVGSurface(file, width, height)
        elseif fmt==".ps"
            surf = CairoPSSurface(file, width, height)
        elseif fmt==".png"
            surf = CairoImageSurface(round(Int, width), round(Int, height), Cairo.FORMAT_ARGB32)
        else
            formats = join(_available_formats, ", ", " and ")
            throw(SerendipException("Cannot save image to format $fmt. Available formats are: $formats"))
        end

        ctx = CairoContext(surf)
        set_antialias(ctx, Cairo.ANTIALIAS_NONE) # ANTIALIAS_DEFAULT, ANTIALIAS_NONE, ANTIALIAS_GRAY, ANTIALIAS_SUBPIXEL
        # set_operator(ctx, Cairo.OPERATOR_SOURCE) # makes to look rasterized

        draw!(figure, ctx)

        if fmt==".png"
            write_to_png(surf, file)
        else
            finish(surf)
        end

        figure.quiet || printstyled("  file $file saved\n", color=:cyan)
    end

end
