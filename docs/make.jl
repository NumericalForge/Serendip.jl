import Pkg

const root = @__DIR__
const repo_root = normpath(joinpath(root, ".."))

# Keep the docs environment in sync with the working tree so local builds
# do not depend on a previously-resolved manifest or a machine-specific path.
Pkg.activate(root)
Pkg.develop(Pkg.PackageSpec(path = repo_root))
Pkg.resolve()
Pkg.instantiate()

using Documenter, Literate, Serendip

include("generate_api_pages.jl")
generate_api_pages(root)

Literate.markdown(
    joinpath(repo_root, "examples", "docs", "simple-truss.jl"),
    joinpath(root, "src", "examples");
    name = "simple-truss",
    documenter = true,
)

Literate.markdown(
    joinpath(repo_root, "examples", "docs", "static-2d.jl"),
    joinpath(root, "src", "examples");
    name = "static-2d",
    documenter = true,
)

makedocs(
    # remotes = nothing,
    root = root,
    modules  = [Serendip],
    sitename = "Serendip",
    pagesonly = true,    # only listed pages are included
    checkdocs = :exports,   # :missing, :exports, :all
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        # assets = ["assets/luxor-docs.css"],
        collapselevel=1,
    ),
    pages = [
        "Introduction" => "index.md",
        "Manual" => [
            "Getting Started" => "manual/getting-started.md",
            "Workflow" => "manual/workflow.md",
            "Tutorials" => "manual/tutorials.md",
        ],
        "Examples" =>  [
            "Overview" => "examples/overview.md",
            "Simple Truss" => "examples/simple-truss.md",
            "Static 2D" => "examples/static-2d.md",
        ],
        "API Reference" => [
            "Overview" => "api/reference.md",
            "Core Structures" => "api/core.md",
            "Geometry" => "api/geometry.md",
            "Shape Functions" => "api/shapes.md",
            "Mesh" => "api/mesh.md",
            "Analyses" => "api/analyses.md",
            "Plotting and Data" => "api/plotting-data.md",
        ]
    ],
    doctest = true,
    repo = Documenter.Remotes.GitHub("NumericalForge", "Serendip.jl"),
    # edit_branch = "main",

    # repo = "github.com/JuliaGraphics/LuxorManual.git"
    # repo = "github.com/NumericalForge/Serendip.jl.git"
)

# ENV["GITHUB_REPOSITORY"] = "NumericalForge/Serendip.jl"
# ENV["GITHUB_EVENT_NAME"] = "push"
# ENV["GITHUB_REF"]        = "main"
# ENV["GITHUB_ACTOR"]      = "usr"
# ENV["GITHUB_TOKEN"]      = raw"${{ secrets.GITHUB_TOKEN }}"

deploydocs(;
    devbranch = "main",
    target    = "build",
    branch    = "gh-pages",
    repo      = "github.com/NumericalForge/Serendip.jl.git",
)
