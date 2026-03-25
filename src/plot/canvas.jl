# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

mutable struct Canvas<:FigureComponent
    frame::Frame
    limits::Vector{Float64}  # [ xmin, ymin, xmax, ymax ] in data coordinates
    function Canvas()
        return new(Frame(), Float64[0.0, 0.0, 0.0, 0.0])
    end
end

function data2user(canvas::Canvas, x::Float64, y::Float64)
    Xmin = canvas.frame.x
    Ymin = canvas.frame.y
    Xmax = canvas.frame.x + canvas.frame.width
    Ymax = canvas.frame.y + canvas.frame.height
    xmin, ymin, xmax, ymax = canvas.limits
    X = Xmin + (Xmax-Xmin)/(xmax-xmin)*(x-xmin)
    Y = Ymin + (Ymax-Ymin)/(ymax-ymin)*(ymax-y)
    return X, Y
end

function user2data(canvas::Canvas, X::Float64, Y::Float64)
    Xmin = canvas.frame.x
    Ymin = canvas.frame.y
    Xmax = canvas.frame.x + canvas.frame.width
    Ymax = canvas.frame.y + canvas.frame.height
    xmin, ymin, xmax, ymax = canvas.limits
    x = xmin + (xmax-xmin)/(Xmax-Xmin)*(X-Xmin)
    y = ymax + (ymax-ymin)/(Ymax-YMin)*(Ymin-Y)
    return x, y
end
