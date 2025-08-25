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
    using Printf, DelimitedFiles, DataStructures, Glob, Dates
    using Cairo, LaTeXStrings, MathTeXEngine
    import FreeTypeAbstraction

    export @L_str # reexport LaTeXStrings

    import DataStructures: OrderedDict, OrderedSet

    # Tools module
    include("tools/include.jl")

    # generic exports
    export max, min, sort, reset, getindex, sort, copy!, show

    abstract type AbstractPoint end
    abstract type AbstractCell end
    # abstract type AbstractBlock<:AbstractCell end
    abstract type AbstractBlock end
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
    export GeoModel

    # Mesh
    include("mesh/include.jl")
    export getcoords

    include("context.jl")
    export Context

    include("material.jl")
    export Material, read_prms

    include("ip.jl")
    export Ip, ip_vals, maximum, minimum, sort

    include("element.jl")
    export Element
    export get_nodes, changequadrature!, get_ips, elems_ip_vals, update_material!, setstate!

    # Plotting
    include("plot/include.jl")
    export Chart, LineChart, DataSeries, Legend, Colormap, DomainPlot, GeometryPlot, Annotation
    export add_series, addlegend!, addannotation!

    include("region-mapper.jl")
    export RegionMapper, RegionModel, add_mapping

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
    export SolverSettings

    include("fe-model.jl")
    export Model
    export FEModel
    export Domain
    export addlogger!, addmonitor!, addloggers!, addmonitors!
    export setloggers!, setmonitors!

    include("solver.jl")

    include("io.jl")

    # Mechanical module
    include("mech/include.jl")

    # Hydromech module
    # include("hydromech/include.jl")

    # ThermoMech module
    # include("thermomech/include.jl")

    # Acoustic module
    # include("acousticmech/include.jl")

    export run

    include("deprecated.jl")

    # show function for FE related types
    Base.show(io::IO, obj::Dof)               = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Node)              = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Ip)                = _show(io, obj, 2, "")
    Base.show(io::IO, obj::IpState)           = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Element)           = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Material)          = _show(io, obj, 2, "")
    Base.show(io::IO, obj::BoundaryCondition) = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Facet)             = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Logger)            = _show(io, obj, 2, "")
    Base.show(io::IO, obj::Monitor)           = _show(io, obj, 2, "")
    Base.show(io::IO, obj::FEModel)           = _show(io, obj, 2, "")


    # testing
    export @runfiles
    macro runfiles(files)
        return esc(quote
            for file in $files
                printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
                include(file)
                println()
            end
        end)
    end

end#module
