abstract type Figure
    # width::Float64
    # height::Float64
end

abstract type FigureComponent end
# abstract type DataSeries end

const cm = 72.0 / 2.54
const _png_raster_scale = 5.0
const _png_stroke_boost = 5.0

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

    function TextBox(text::AbstractString=""; frame::Frame=Frame(), angle::Real=0.0)
        return new(text, frame, float(angle))
    end
end

struct RenderContext
    cairo_ctx::CairoContext
    backend::Symbol
    base_matrix::CairoMatrix
    background::Union{Nothing,Color}
    width_scale::Float64
end

_identity_matrix() = CairoMatrix([1, 0, 0, 1, 0, 0]...)

function RenderContext(
    cairo_ctx::CairoContext;
    backend::Symbol=:custom,
    base_matrix::CairoMatrix=_identity_matrix(),
    background=nothing,
    width_scale::Real=1.0,
)
    return RenderContext(cairo_ctx, backend, base_matrix, resolve_color(background), float(width_scale))
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

function _compose_matrix(A::CairoMatrix, B::CairoMatrix)
    return CairoMatrix(
        A.xx * B.xx + A.xy * B.yx,
        A.yx * B.xx + A.yy * B.yx,
        A.xx * B.xy + A.xy * B.yy,
        A.yx * B.xy + A.yy * B.yy,
        A.xx * B.x0 + A.xy * B.y0 + A.x0,
        A.yx * B.x0 + A.yy * B.y0 + A.y0,
    )
end

reset_matrix!(ctx::RenderContext) = set_matrix(ctx.cairo_ctx, ctx.base_matrix)

function set_local_matrix!(ctx::RenderContext, local_matrix::CairoMatrix)
    set_matrix(ctx.cairo_ctx, _compose_matrix(ctx.base_matrix, local_matrix))
end

_figure_background(::Figure) = nothing

function _draw_figure_background!(ctx::RenderContext, frame::Frame, color)
    color === nothing && return

    cairo_ctx = ctx.cairo_ctx
    reset_matrix!(ctx)
    rectangle(cairo_ctx, frame.x, frame.y, frame.width, frame.height)
    set_source_rgba(cairo_ctx, rgba(color)...)
    fill(cairo_ctx)
end

function _measure_text_box(box::TextBox, cc::CairoContext, fontsize::Float64)
    isempty(box.text) && return 0.0, 0.0
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

function _draw_text_box!(ctx::RenderContext, box::TextBox)
    isempty(box.text) && return
    cairo_ctx = ctx.cairo_ctx
    x = box.frame.x + 0.5 * box.frame.width
    y = box.frame.y + 0.5 * box.frame.height
    draw_text(cairo_ctx, x, y, box.text, halign="center", valign="center", angle=box.angle)
end

draw_background!(::Figure, ::RenderContext) = nothing
draw_contents!(::Figure, ::RenderContext) = throw(MethodError(draw_contents!, (Figure, RenderContext)))

function draw!(figure::Figure, cairo_ctx::CairoContext)
    ctx = RenderContext(cairo_ctx; background=_figure_background(figure), width_scale=1.0)
    return draw!(figure, ctx)
end

function draw!(figure::Figure, ctx::RenderContext)
    reset_matrix!(ctx)
    draw_background!(figure, ctx)
    draw_contents!(figure, ctx)
    return nothing
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
            backend = :pdf
            base_matrix = _identity_matrix()
        elseif fmt==".svg"
            surf = CairoSVGSurface(file, width, height)
            backend = :svg
            base_matrix = _identity_matrix()
        elseif fmt==".ps"
            surf = CairoPSSurface(file, width, height)
            backend = :ps
            base_matrix = _identity_matrix()
        elseif fmt==".png"
            surf = CairoImageSurface(round(Int, _png_raster_scale * width), round(Int, _png_raster_scale * height), Cairo.FORMAT_ARGB32)
            backend = :png
            base_matrix = CairoMatrix([_png_raster_scale, 0, 0, _png_raster_scale, 0, 0]...)
        else
            formats = join(_available_formats, ", ", " and ")
            throw(SerendipException("Cannot save image to format $fmt. Available formats are: $formats"))
        end

        cairo_ctx = CairoContext(surf)
        render_background = backend == :png ? something(_figure_background(figure), Color(:white)) : _figure_background(figure)
        width_scale = backend == :png ? _png_stroke_boost : 1.0
        ctx = RenderContext(cairo_ctx; backend=backend, base_matrix=base_matrix, background=render_background, width_scale=width_scale)
        reset_matrix!(ctx)
        antialias = backend == :png ? Cairo.ANTIALIAS_DEFAULT : Cairo.ANTIALIAS_NONE
        set_antialias(cairo_ctx, antialias) # ANTIALIAS_DEFAULT, ANTIALIAS_NONE, ANTIALIAS_GRAY, ANTIALIAS_SUBPIXEL
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
