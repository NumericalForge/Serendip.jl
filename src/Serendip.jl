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
    using StaticArrays, FixedSizeArrays, SparseArrays, Arpack, Gmsh
    using Printf, DelimitedFiles, DataStructures, Glob, Dates, JSON
    using Cairo, LaTeXStrings, MathTeXEngine
    import FreeTypeAbstraction

    export @L_str # reexport LaTeXStrings

    import DataStructures: OrderedDict, OrderedSet

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
    export Node, add_dof, get_dof, get_values, setvalue!

    # Shapes
    include("shape/include.jl")

    # Geometry
    include("geo/include.jl")
    export GeoModel, Block, Path, GPath, Point, Edge, Surface, Volume

    # Mesh
    include("mesh/include.jl")
    export getcoords

    include("context.jl")
    export Context

    abstract type Constitutive end
    export Constitutive

    include("ip.jl")
    export Ip, ip_vals, maximum, minimum, sort

    include("element.jl")
    export Element
    export get_nodes, changequadrature!, get_ips, elems_ip_vals, update_material!, set_state

    # Plotting
    include("plot/include.jl")
    export Chart, LineChart, DataSeries, Legend, Colormap, DomainPlot, GeometryPlot, Annotation
    export add_series, addlegend!, addannotation!

    include("region-mapper.jl")
    export RegionMapper, RegionModel, add_mapping, add_map

    abstract type Analysis end

    include("stage.jl")
    export Stage, add_stage

    # Boundary conditions
    include("bc.jl")
    export add_bc, add_body_load

    include("logger.jl")
    export add_logger

    include("monitor.jl")
    export add_monitor

    include("analysis.jl")
    
    include("fe-model.jl")
    export Model
    export FEModel
    export Domain
    export addlogger!, addmonitor!, addloggers!, addmonitors!
    export setloggers!, setmonitors!
    
    include("solver.jl")
    export SolverSettings, run

    include("io.jl")

    # Mechanical module
    export reset_displacements
    include("mech/include.jl")

    # Hydromech module
    # include("hydromech/include.jl")

    # ThermoMech module
    # include("thermomech/include.jl")

    # Acoustic module
    # include("acousticmech/include.jl")

    include("deprecated.jl")

    # show functions for goemetry related types
    Base.show(io::IO, obj::GeoModel) = _show(io, obj, 2, "")
    Base.show(io::IO, obj::GPath)    = _show(io, obj, 3, "")
    Base.show(io::IO, obj::Path)     = _show(io, obj, 3, "")
    Base.show(io::IO, obj::PathCmd)  = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Point)    = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Edge)     = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Surface)  = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Volume)   = _show(io, obj, 2, "")

    # show function for mesh related types
    Base.show(io::IO, obj::CellShape) = _show(io, obj, 2)
    Base.show(io::IO, obj::Mesh)      = _show(io, obj, 2)
    Base.show(io::IO, obj::Block)     = _show(io, obj, 2)

    # show function for FE related types
    Base.show(io::IO, obj::Dof)          = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Node)         = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Ip)           = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Facet)        = _show(io, obj, 2, "")
    Base.show(io::IO, obj::IpState)      = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Element)      = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Constitutive) = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Logger)       = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Monitor)      = _show(io, obj, 2, "")
    Base.show(io::IO, obj::RegionMapper) = _show(io, obj, 4, "")
    Base.show(io::IO, obj::Context)      = _show(io, obj, 2, "")
    Base.show(io::IO, obj::FEModel)      = _show(io, obj, 2, "")

    Base.show(io::IO, obj::BoundaryCondition) = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Analysis)          = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Stage)             = _show(io, obj, 2, "")



    # testing
    export @run_files
    macro run_files(files)
        return esc(quote
            for file in $files
                printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
                include(file)
                println()
            end
        end)
    end

end#module
