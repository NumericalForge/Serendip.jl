# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

"""
    ChartGrid(; title="", size=(420,320), font="NewComputerModern", font_size=10.0,
                background=nothing,
                column_headers=String[], row_headers=String[],
                outerpad=0.0, hgap=8.0, vgap=8.0, quiet=false)

Create a composed figure that arranges child figures in a uniform grid.

# Keywords
- `title::AbstractString`: grid title centered above the cells.
- `size::Tuple{<:Real,<:Real}`: figure size in points. Use `cm` as a convenience helper, e.g. `size=(16cm, 10cm)`.
- `font::AbstractString`: font family used for the grid title.
- `font_size::Real`: title font size.
- `background`: full-figure background fill; `nothing` leaves the figure unfilled.
- `column_headers::Vector{<:AbstractString}`: optional headers centered above columns.
- `row_headers::Vector{<:AbstractString}`: optional headers drawn on the left, rotated like y-axis labels.
- `outerpad::Real`: outer figure margin.
- `hgap::Real`, `vgap::Real`: horizontal and vertical gaps between cells.
- `quiet::Bool`: suppress constructor log.

# Notes
- Add child figures with `add_chart`.
- The grid can contain any `Figure`, including `Chart`, `DomainPlot`, or another `ChartGrid`.
- `background=nothing` leaves the grid background transparent in PNG and unfilled in vector outputs.
- Child chart backgrounds are ignored when charts are rendered inside a grid.
- Headers and titles use the same math-aware text rendering as the rest of the plotting API.
"""
mutable struct ChartGrid <: Figure
    width::Float64
    height::Float64
    figure_frame::Frame
    title_box::TextBox
    font::String
    font_size::Float64
    background::Union{Nothing,Color}
    column_header_boxes::Vector{TextBox}
    row_header_boxes::Vector{TextBox}
    outerpad::Float64
    hgap::Float64
    vgap::Float64
    children::Dict{Tuple{Int,Int},Figure}
    cell_frames::Dict{Tuple{Int,Int},Frame}
    nrows::Int
    ncols::Int
    quiet::Bool

    function ChartGrid(;
        title::AbstractString="",
        size::Tuple{<:Real,<:Real}=(420, 320),
        font::AbstractString="NewComputerModern",
        font_size::Real=10.0,
        background=nothing,
        column_headers::Vector{<:AbstractString}=String[],
        row_headers::Vector{<:AbstractString}=String[],
        outerpad::Real=0.0,
        hgap::Real=8.0,
        vgap::Real=8.0,
        quiet::Bool=false,
    )
        @check font_size > 0 "ChartGrid: font_size must be positive"
        @check hgap >= 0 "ChartGrid: hgap must be non-negative"
        @check vgap >= 0 "ChartGrid: vgap must be non-negative"

        width, height = size
        background = resolve_color(background)
        title_box = TextBox(title)
        column_header_boxes = TextBox[TextBox(String(s)) for s in column_headers]
        row_header_boxes = TextBox[TextBox(String(s); angle=90.0) for s in row_headers]
        return new(
            width,
            height,
            Frame(0.0, 0.0, width, height),
            title_box,
            String(font),
            float(font_size),
            background,
            column_header_boxes,
            row_header_boxes,
            float(outerpad),
            float(hgap),
            float(vgap),
            Dict{Tuple{Int,Int},Figure}(),
            Dict{Tuple{Int,Int},Frame}(),
            0,
            0,
            quiet,
        )
    end
end


"""
    add_chart(grid::ChartGrid, child::Figure, pos::Tuple{Int,Int})

Insert `child` into `grid` at position `(row, col)`.

# Arguments
- `grid::ChartGrid`: target grid.
- `child::Figure`: child figure to place in the grid.
- `pos::Tuple{Int,Int}`: one-based `(row, col)` cell coordinates.

# Notes
- Grid dimensions grow automatically to accommodate the largest inserted row/column.
- Inserting into an occupied cell throws an error.
"""
function add_chart(grid::ChartGrid, child::Figure, pos::Tuple{Int,Int})
    i, j = pos
    @check i > 0 && j > 0 "ChartGrid: cell indices must be positive"
    @check !haskey(grid.children, pos) "ChartGrid: position $pos is already occupied"

    grid.children[pos] = child
    grid.nrows = max(grid.nrows, i)
    grid.ncols = max(grid.ncols, j)
    return child
end


function _measure_grid_title(grid::ChartGrid)
    surf = CairoImageSurface(4, 4, Cairo.FORMAT_ARGB32)
    cc = CairoContext(surf)
    select_font_face(cc, get_font(grid.font), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(cc, grid.font_size)
    return _measure_text_box(grid.title_box, cc, grid.font_size)
end


function _measure_grid_headers(grid::ChartGrid)
    surf = CairoImageSurface(4, 4, Cairo.FORMAT_ARGB32)
    cc = CairoContext(surf)
    select_font_face(cc, get_font(grid.font), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(cc, grid.font_size)

    ncol_headers = min(length(grid.column_header_boxes), grid.ncols)
    nrow_headers = min(length(grid.row_header_boxes), grid.nrows)

    col_height = 0.0
    for i in 1:ncol_headers
        _, h = _measure_text_box_extent(grid.column_header_boxes[i], cc, grid.font_size)
        col_height = max(col_height, h)
    end

    row_width = 0.0
    for i in 1:nrow_headers
        w, _ = _measure_text_box_extent(grid.row_header_boxes[i], cc, grid.font_size)
        row_width = max(row_width, w)
    end

    return col_height, row_width
end


function configure!(grid::ChartGrid)
    length(grid.children) > 0 || throw(SerendipException("ChartGrid: no child figures added"))

    grid.figure_frame = Frame(0.0, 0.0, grid.width, grid.height)
    grid.outerpad = max(0.01 * min(grid.width, grid.height), grid.outerpad)

    title_width, title_height = _measure_grid_title(grid)
    column_header_height, row_header_width = _measure_grid_headers(grid)
    title_gap = _has_text(grid.title_box) ? 0.6 * grid.font_size : 0.0
    column_header_gap = column_header_height > 0 ? 0.5 * grid.font_size : 0.0
    row_header_gap = row_header_width > 0 ? 0.5 * grid.font_size : 0.0

    content_x = grid.figure_frame.x + grid.outerpad + row_header_width + row_header_gap
    content_y = grid.figure_frame.y + grid.outerpad + title_height + title_gap + column_header_height + column_header_gap
    content_width = grid.width - 2 * grid.outerpad - row_header_width - row_header_gap
    content_height = grid.height - 2 * grid.outerpad - title_height - title_gap - column_header_height - column_header_gap

    @check grid.nrows > 0 && grid.ncols > 0 "ChartGrid: invalid grid dimensions"
    @check content_width > 0 && content_height > 0 "ChartGrid: insufficient space for grid content"

    cell_width = (content_width - (grid.ncols - 1) * grid.hgap) / grid.ncols
    cell_height = (content_height - (grid.nrows - 1) * grid.vgap) / grid.nrows
    @check cell_width > 0 && cell_height > 0 "ChartGrid: insufficient space for grid cells"

    grid.title_box.frame = Frame(content_x, grid.figure_frame.y + grid.outerpad, content_width, title_height)
    grid.title_box.angle = 0.0
    grid.title_box.visible = !isempty(grid.title_box.text)
    empty!(grid.cell_frames)

    for i in 1:grid.nrows, j in 1:grid.ncols
        x = content_x + (j - 1) * (cell_width + grid.hgap)
        y = content_y + (i - 1) * (cell_height + grid.vgap)
        grid.cell_frames[(i, j)] = Frame(x, y, cell_width, cell_height)
    end

    if column_header_height > 0
        y = grid.figure_frame.y + grid.outerpad + title_height + title_gap
        for j in 1:grid.ncols
            j <= length(grid.column_header_boxes) || continue
            x = content_x + (j - 1) * (cell_width + grid.hgap)
            box = grid.column_header_boxes[j]
            box.frame = Frame(x, y, cell_width, column_header_height)
            box.angle = 0.0
            box.visible = !isempty(box.text)
        end
    end

    if row_header_width > 0
        x = grid.figure_frame.x + grid.outerpad
        for i in 1:grid.nrows
            i <= length(grid.row_header_boxes) || continue
            y = content_y + (i - 1) * (cell_height + grid.vgap)
            box = grid.row_header_boxes[i]
            box.frame = Frame(x, y, row_header_width, cell_height)
            box.angle = 90.0
            box.visible = !isempty(box.text)
        end
    end
end


function _set_grid_frame!(child::Chart, frame::Frame)
    child.width = frame.width
    child.height = frame.height
    child.figure_frame = Frame(frame.x, frame.y, frame.width, frame.height)
end


function _set_grid_frame!(child::DomainPlot, frame::Frame)
    child.width = frame.width
    child.height = frame.height
    child.figure_frame = Frame(frame.x, frame.y, frame.width, frame.height)
end


function _set_grid_frame!(child::ChartGrid, frame::Frame)
    child.width = frame.width
    child.height = frame.height
    child.figure_frame = Frame(frame.x, frame.y, frame.width, frame.height)
end


function _draw_grid_child!(ctx::CairoContext, child::Figure, frame::Frame, parent_background=nothing)
    old_width = child.width
    old_height = child.height
    old_frame = Frame(child.figure_frame.x, child.figure_frame.y, child.figure_frame.width, child.figure_frame.height)
    old_background = child isa Union{Chart,ChartGrid} ? child.background : nothing
    old_legend_background = child isa Chart ? child.legend.background : nothing

    try
        _set_grid_frame!(child, frame)
        if child isa Union{Chart,ChartGrid}
            child.background = nothing
        end
        if child isa Chart && child.legend.background === nothing && parent_background !== nothing
            child.legend.background = parent_background
        end
        configure!(child)
        draw!(child, ctx)
    finally
        child.width = old_width
        child.height = old_height
        child.figure_frame = old_frame
        if child isa Union{Chart,ChartGrid}
            child.background = old_background
        end
        if child isa Chart
            child.legend.background = old_legend_background
        end
    end
end


function draw!(grid::ChartGrid, ctx::CairoContext)
    _draw_figure_background!(ctx, grid.figure_frame, grid.background)

    set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
    select_font_face(ctx, get_font(grid.font), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(ctx, grid.font_size)
    set_source_rgb(ctx, 0.0, 0.0, 0.0)

    _draw_text_box!(ctx, grid.title_box)

    for j in 1:min(length(grid.column_header_boxes), grid.ncols)
        _draw_text_box!(ctx, grid.column_header_boxes[j])
    end

    for i in 1:min(length(grid.row_header_boxes), grid.nrows)
        _draw_text_box!(ctx, grid.row_header_boxes[i])
    end

    for pos in sort(collect(keys(grid.children)))
        _draw_grid_child!(ctx, grid.children[pos], grid.cell_frames[pos], grid.background)
    end
end
