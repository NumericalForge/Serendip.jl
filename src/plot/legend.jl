# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

const _legend_locations=[
    :right,
    :left,
    :top,
    :bottom,
    :top_right,
    :top_left,
    :bottom_right,
    :bottom_left,
    :outer_right,
    :outer_left,
    :outer_top,
    :outer_top_left,
    :outer_top_right,
    :outer_bottom,
    :outer_bottom_right,
    :outer_bottom_left
]


mutable struct Legend<:FigureComponent
    location::Symbol
    font::String
    font_size::Float64
    ncols::Int
    handle_length::Float64  # length of the line
    row_sep::Float64      # separation between labels
    col_sep::Float64      # separation between labels
    inner_pad::Float64      # padding inside the legend box
    outer_pad::Float64      # padding outside the legend box
    width::Float64          # width of the legend box
    height::Float64         # height of the legend box

    function Legend(;
        location::Symbol = :top_right,
        font::String = "NewComputerModern",
        font_size::Float64 = 7.0,
        ncols::Int = 1,
        # title::String = ""
    )
        @check location in _legend_locations "location must be one of $(_legend_locations)"
        @check font_size>0 "font_size must be positive"
        @check ncols>0 "ncols must be positive"

        return new(
            location,
            font,
            font_size,
            ncols,
        )
    end
end

