using Serendip

const API_PAGES = [
    (
        file = "core.md",
        title = "Core Structures",
        sections = [
            ("Module Helpers", ["@t_str", "@L_str", "@run_files"]),
            ("Basic FEM Objects", ["Context", "Constitutive", "Dof", "Node", "Ip", "Element"]),
            ("Node and Element Utilities", [
                "add_dof",
                "get_dof",
                "get_values",
                "change_quadrature",
                "get_ips",
                "nearest",
                "get_coords",
            ]),
            ("Model Assembly", ["RegionMapper", "RegionModel", "Stage", "FEModel", "SolverSettings"]),
            ("Workflow Helpers", [
                "add_mapping",
                "add_map",
                "add_stage",
                "add_bc",
                "add_logger",
                "add_monitor",
                "run",
                "reset_displacements",
                "save",
            ]),
        ],
    ),
    (
        file = "geometry.md",
        title = "Geometry",
        sections = [
            ("Geometry Objects", ["GeoModel", "Block", "Path", "GPath", "Point", "Edge", "Surface", "Volume"]),
            ("Geometry Creation", [
                "add_path",
                "add_block",
                "add_array",
                "add_point",
                "add_line",
                "add_circle_arc",
                "add_circle",
                "add_bezier",
                "add_loop",
                "add_wire",
                "add_plane_surface",
                "add_polygon",
                "add_surface_loop",
                "add_surface_filling",
                "add_volume",
                "add_box",
                "add_cylinder",
                "add_sphere",
                "add_rectangle",
                "add_disk",
            ]),
            ("Geometry Queries and Editing", [
                "get_boundary",
                "get_entities",
                "get_points",
                "get_curves",
                "get_surfaces",
                "get_volumes",
                "get_point",
                "get_curve",
                "get_surface",
                "get_volume",
                "translate",
                "rotate",
                "extrude",
                "revolve",
                "mirror",
                "cut",
                "cut!",
                "fuse",
                "intersect",
                "fragment",
                "fillet",
                "set_size",
                "set_refinement",
                "set_transfinite_curve",
                "set_transfinite_surface",
                "set_recombine",
                "set_transfinite_volume",
            ]),
        ],
    ),
    (
        file = "shapes.md",
        title = "Shape Functions",
        sections = [
            ("Shape and Interpolation Helpers", [
                "VTKCellType",
                "get_ip_coords",
                "get_shape",
                "inverse_map",
                "extrapolator",
            ]),
        ],
    ),
    (
        file = "mesh.md",
        title = "Mesh",
        sections = [
            ("Mesh Objects", ["Cell", "Mesh"]),
            ("Mesh Queries", [
                "select",
                "get_coords",
                "get_nodes",
                "get_facets",
                "get_patches",
                "cell_extent",
                "cell_quality",
                "cell_aspect_ratio",
                "get_outer_facets",
                "threshold",
                "get_feature_edges",
                "get_feature_mesh",
                "stats",
                "get_segment_data",
            ]),
            ("Mesh Editing and Generation", [
                "hrefine",
                "prefine",
                "move",
                "array",
                "copy",
                "mirror",
                "rotate!",
                "polar",
                "scale",
                "permute_coordinates",
                "extrude",
                "revolve",
                "slice",
                "smooth!",
                "laplacian_smooth!",
                "fast_smooth!",
                "fast_smooth2!",
                "add_cohesive_elements",
                "add_contact_elements",
                "add_boundary_contact_elements",
                "add_boundary_shell_elements",
                "load_cracked_mesh",
                "rand_mesh",
                "save",
            ]),
        ],
    ),
    (
        file = "analyses.md",
        title = "Analyses",
        sections = [
            ("Mechanical Analyses", ["MechAnalysis", "MechModalAnalysis", "DynamicAnalysis"]),
            ("Mechanical Elements", [
                "MechBar",
                "MechBeam",
                "MechBondSlip",
                "MechBondTip",
                "MechBulk",
                "MechCohesive",
                "MechContact",
                "MechFrame",
                "MechShell",
                "MechSolid",
            ]),
            ("Mechanical Constitutive Models", [
                "MechIntegrator",
                "LinearElastic",
                "LinearElasticFluid",
                "LinearInterface",
                "LinearContact",
                "LinearTip",
                "ElasticSpring",
                "LinearBondSlip",
                "ElasticBondSlip",
                "ElasticRSJoint",
                "CyclicBondSlip",
                "PowerExpBondSlip",
                "CebBondSlip",
                "LinearCohesive",
                "MohrCoulombCohesive",
                "PowerYieldCohesive",
                "AsinhYieldCohesive",
                "CoulombContact",
                "MohrCoulombContact",
                "DruckerPrager",
                "VonMises",
                "WillamWarnke",
                "ECP",
                "UCP",
                "LumpedMass",
                "update_state",
            ]),
            ("Thermal and Thermo-Mechanical", [
                "ThermoContext",
                "ThermoAnalysis",
                "ThermoMechAnalysis",
                "ThermoSolid",
                "TMSolid",
                "ConstConductivity",
                "LinearElasticThermo",
            ]),
            ("Acoustic and Coupled Acoustic", [
                "AcousticContext",
                "AcousticAnalysis",
                "AcousticMechAnalysis",
                "AcousticModalAnalysis",
                "AcousticFluid",
                "AcousticMechInterface",
                "AcousticMechInterfaceCoupling",
                "LinearAcousticFluid",
            ]),
        ],
    ),
    (
        file = "plotting-data.md",
        title = "Plotting and Data",
        sections = [
            ("Plotting", [
                "Color",
                "lighten",
                "darken",
                "Chart",
                "DataSeries",
                "Legend",
                "Colormap",
                "DomainPlot",
                "Annotation",
                "add_series",
                "add_line",
                "add_scatter",
                "add_bar",
                "add_annotation",
            ]),
            ("Data Containers", [
                "DataTable",
                "DataBook",
                "split_by",
                "LinearSpline",
                "CubicSpline",
                "resize",
                "filter",
                "clamp!",
                "denoise!",
                "getheader",
                "get_frequency",
                "randtable",
                "randbook",
            ]),
            ("Symbolic and XML Utilities", [
                "Symbolic",
                "evaluate",
                "XmlDocument",
                "XmlElement",
                "to_xml_node",
                "@define",
                "@check",
                "x",
                "y",
                "z",
                "t",
                "and",
                "or",
            ]),
            ("General Utilities", [
                "_show",
                "copy!",
                "getindex",
                "max",
                "min",
                "maximum",
                "minimum",
                "sort",
                "show",
                "push!",
                "stop",
                "sound_alert",
            ]),
        ],
    ),
]

function has_runtime_binding(name::String)
    try
        getfield(Serendip, Symbol(name))
        return true
    catch
        return false
    end
end

function get_runtime_binding(name::String)
    try
        return getfield(Serendip, Symbol(name))
    catch
        return nothing
    end
end

function render_item(name::String; canonical::Bool)
    obj = get_runtime_binding(name)
    if obj !== nothing && !startswith(name, "@")
        directive = canonical ? "@autodocs" : "@autodocs; canonical=false"
        return """
        ```@raw html
        <details>
        <summary><code>$name</code></summary>
        ```

        ```$directive
        Modules = [Serendip]
        Private = false
        Public = true
        Filter = x -> x === Serendip.$name
        ```

        ```@raw html
        </details>
        ```
        """
    end

    return "- `$name`"
end

function render_page(title::String, sections, seen)
    parts = String[
        "```@meta",
        "CurrentModule = Serendip",
        "```",
        "",
        "# $title",
        "",
    ]

    for (heading, items) in sections
        push!(parts, "## $heading", "")
        for item in items
            obj = get_runtime_binding(item)
            key = obj === nothing ? item : obj
            canonical = !any(existing -> existing === key, seen)
            push!(parts, render_item(item; canonical), "")
            push!(seen, key)
        end
    end

    return join(parts, "\n")
end

function generate_api_pages(root::String)
    api_dir = joinpath(root, "src", "api")
    seen = Any[]
    for page in API_PAGES
        content = render_page(page.title, page.sections, seen)
        write(joinpath(api_dir, page.file), content)
    end
end
