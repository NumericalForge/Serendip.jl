# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

"""
    ChartGrid(; title="", size=(420,320), font="NewComputerModern", font_size=10.0,
                outerpad=0.0, hgap=8.0, vgap=8.0, quiet=false)

Create a composed figure that arranges child figures in a uniform grid.

# Keywords
- `title::AbstractString`: grid title centered above the cells.
- `size::Tuple{Int,Int}`: figure size in points.
- `font::AbstractString`: font family used for the grid title.
- `font_size::Real`: title font size.
- `outerpad::Real`: outer figure margin.
- `hgap::Real`, `vgap::Real`: horizontal and vertical gaps between cells.
- `quiet::Bool`: suppress constructor log.

# Notes
- Add child figures with `add_chart`.
- The grid can contain any `Figure`, including `Chart`, `DomainPlot`, or another `ChartGrid`.
"""
mutable struct ChartGrid <: Figure
    width::Float64
    height::Float64
    figure_frame::Frame
    title::AbstractString
    title_frame::Frame
    font::String
    font_size::Float64
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
        size::Tuple{Int,Int}=(420, 320),
        font::AbstractString="NewComputerModern",
        font_size::Real=10.0,
        outerpad::Real=0.0,
        hgap::Real=8.0,
        vgap::Real=8.0,
        quiet::Bool=false,
    )
        @check font_size > 0 "ChartGrid: font_size must be positive"
        @check hgap >= 0 "ChartGrid: hgap must be non-negative"
        @check vgap >= 0 "ChartGrid: vgap must be non-negative"

        width, height = size
        return new(
            width,
            height,
            Frame(0.0, 0.0, width, height),
            title,
            Frame(),
            String(font),
            float(font_size),
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
    if isempty(grid.title)
        return 0.0, 0.0
    end

    surf = CairoImageSurface(4, 4, Cairo.FORMAT_ARGB32)
    cc = CairoContext(surf)
    select_font_face(cc, get_font(grid.font), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(cc, grid.font_size)
    return getsize(cc, grid.title, grid.font_size)
end


function configure!(grid::ChartGrid)
    length(grid.children) > 0 || throw(SerendipException("ChartGrid: no child figures added"))

    grid.figure_frame = Frame(0.0, 0.0, grid.width, grid.height)
    grid.outerpad = max(0.01 * min(grid.width, grid.height), grid.outerpad)

    title_width, title_height = _measure_grid_title(grid)
    title_gap = isempty(grid.title) ? 0.0 : 0.6 * grid.font_size

    content_x = grid.figure_frame.x + grid.outerpad
    content_y = grid.figure_frame.y + grid.outerpad + title_height + title_gap
    content_width = grid.width - 2 * grid.outerpad
    content_height = grid.height - 2 * grid.outerpad - title_height - title_gap

    @check grid.nrows > 0 && grid.ncols > 0 "ChartGrid: invalid grid dimensions"
    @check content_width > 0 && content_height > 0 "ChartGrid: insufficient space for grid content"

    cell_width = (content_width - (grid.ncols - 1) * grid.hgap) / grid.ncols
    cell_height = (content_height - (grid.nrows - 1) * grid.vgap) / grid.nrows
    @check cell_width > 0 && cell_height > 0 "ChartGrid: insufficient space for grid cells"

    grid.title_frame = Frame(content_x, grid.figure_frame.y + grid.outerpad, content_width, title_height)
    empty!(grid.cell_frames)

    for i in 1:grid.nrows, j in 1:grid.ncols
        x = content_x + (j - 1) * (cell_width + grid.hgap)
        y = content_y + (i - 1) * (cell_height + grid.vgap)
        grid.cell_frames[(i, j)] = Frame(x, y, cell_width, cell_height)
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


function _draw_grid_child!(ctx::CairoContext, child::Figure, frame::Frame)
    old_width = child.width
    old_height = child.height
    old_frame = Frame(child.figure_frame.x, child.figure_frame.y, child.figure_frame.width, child.figure_frame.height)

    _set_grid_frame!(child, frame)
    configure!(child)
    draw!(child, ctx)

    child.width = old_width
    child.height = old_height
    child.figure_frame = old_frame
end


function draw!(grid::ChartGrid, ctx::CairoContext)
    if !isempty(grid.title)
        set_matrix(ctx, CairoMatrix([1, 0, 0, 1, 0, 0]...))
        select_font_face(ctx, get_font(grid.font), Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
        set_font_size(ctx, grid.font_size)
        set_source_rgb(ctx, 0.0, 0.0, 0.0)
        x = grid.title_frame.x + 0.5 * grid.title_frame.width
        y = grid.title_frame.y + 0.5 * grid.title_frame.height
        draw_text(ctx, x, y, grid.title, halign="center", valign="center", angle=0)
    end

    for pos in sort(collect(keys(grid.children)))
        _draw_grid_child!(ctx, grid.children[pos], grid.cell_frames[pos])
    end
end
