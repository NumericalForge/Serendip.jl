# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl

__precompile__()

"""
**Serendip.jl**

Serendip module implements functions and types to perform finite element analyses.

**Important data types**

Node, Element, Model, Dof, Ip, NodeBC, SurfaceBC

"""
module Serendip

using StatsBase, Statistics, LinearAlgebra
using StaticArrays, FixedSizeArrays, SparseArrays, Arpack, Pardiso
using Printf, DelimitedFiles, Glob, Dates, JSON
using Gmsh
using Cairo
using QuickCharts: Chart, ChartGrid, DataSeries, Legend, Annotation, Color, Colormap, Colorbar, VideoBuilder, cm
using QuickCharts: render, lighten, darken, gray
import QuickCharts: save, Figure, FigureComponent, Frame, TextBox, RenderContext, Canvas, Axis
import QuickCharts: configure!, draw!, draw_background!, draw_contents!, reset_matrix!, set_local_matrix!
import QuickCharts: draw_text, getsize, get_font, draw_mark, resolve_color, rgb, rgba
import QuickCharts: data2user, user2data, _draw_figure_background!, _draw_text_box!
import QuickCharts: compute_auto_limits, _figure_renderable, _png_raster_scale
import QuickCharts: _capture_scaling_state, _apply_scaling_state!
import QuickCharts: resize
import QuickCharts: add_series, add_line, add_scatter, add_bar, add_annotation, add_chart, add_frame
import FreeTypeAbstraction
import DataStructures: OrderedDict, OrderedSet

macro t_str(s)
    return s
end

macro L_str(s)
    return esc(:(@t_str($s)))
end

export @t_str, @L_str


# Tools module
include("tools/include.jl")

# generic exports
export max, min, sort, getindex, sort, copy!, show

abstract type AbstractCell end
abstract type AbstractDomain end

# FEM
include("dof.jl")
export Dof

include("node.jl")
export Node, add_dof, get_dof, get_values

# Shapes
include("shape/include.jl")

# Geometry
include("geo/include.jl")
export GeoModel, Block, Path, GPath, Point, Edge, Surface, Volume

# Mesh
include("mesh/include.jl")
export get_coords, remove_elements

include("context.jl")
export Context

include("constitutive.jl")
export Constitutive

include("ip.jl")
export Ip, maximum, minimum, sort, get_ips

include("element.jl")
export Element
export get_nodes, change_quadrature, get_ips

include("select.jl")

# Plotting
include("plot/include.jl")
export cm
export Color
export gray, lighten, darken
export render
export Chart, ChartGrid, DataSeries, Legend, Colormap, DomainPlot, Annotation, VideoBuilder
export add_line, add_scatter, add_bar, add_annotation
export add_chart, add_plot, add_frame
# Deprecated/Aliases
export add_series

include("region-mapper.jl")
export RegionMapper, RegionModel, add_mapping, add_map

abstract type Analysis end

# Boundary conditions and constraints
include("bc.jl")
include("constraint.jl")

export add_bc, add_constraint

include("stage.jl")
export Stage, add_stage

include("logger.jl")
export add_logger

include("monitor.jl")
export add_monitor

include("analysis.jl")

include("fe-model.jl")
export FEModel

include("solver.jl")
export SolverSettings, run

include("io.jl")

# Mechanical module
export reset_displacements
include("mech/include.jl")

# Hydromech module
# include("hydromech/include.jl")

# ThermoMech module
include("thermomech/include.jl")

# Acoustic module
# include("acousticmech/include.jl")

include("deprecated.jl")

# show functions for goemetry related types
Base.show(io::IO, obj::GeoModel) = _show(io, obj, 2, "")
Base.show(io::IO, obj::GPath) = _show(io, obj, 3, "")
Base.show(io::IO, obj::Path) = _show(io, obj, 3, "")
Base.show(io::IO, obj::PathCmd) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Point) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Edge) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Surface) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Volume) = _show(io, obj, 2, "")

# show function for mesh related types
Base.show(io::IO, obj::CellShape) = _show(io, obj, 2)
Base.show(io::IO, obj::Mesh) = _show(io, obj, 2)
Base.show(io::IO, obj::Block) = _show(io, obj, 2)

# show function for FE related types
Base.show(io::IO, obj::Dof) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Node) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Ip) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Facet) = _show(io, obj, 2, "")
Base.show(io::IO, obj::ConstState) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Element) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Constitutive) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Logger) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Monitor) = _show(io, obj, 2, "")
Base.show(io::IO, obj::RegionMapper) = _show(io, obj, 4, "")
Base.show(io::IO, obj::Context) = _show(io, obj, 2, "")
Base.show(io::IO, obj::FEModel) = _show(io, obj, 2, "")

Base.show(io::IO, obj::BoundaryCondition) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Analysis) = _show(io, obj, 2, "")
Base.show(io::IO, obj::Stage) = _show(io, obj, 2, "")



end#module
