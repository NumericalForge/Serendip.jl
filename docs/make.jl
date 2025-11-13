# Inside make.jl
push!(LOAD_PATH,"../src/")
using Serendip, Documenter

root = joinpath(dirname(pathof(Serendip)), "..", "docs")

makedocs(
    # remotes = nothing,
    root = root,
    modules  = [Serendip],
    sitename = "Serendip",
    pagesonly = true,    # only listed pages are included
    checkdocs = :all,   # :missing, :all
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        # assets = ["assets/luxor-docs.css"],
        collapselevel=1,
    ),
    pages = [
        "Introduction" => "index.md",
        # "Tutorial" => "tutorial/tutorial.md",
        # "Mesh" =>  "meshing/meshing.md",
        # "Model model" =>  "modelling/modelling.md",
        "Mechanical analysis" =>  "modelling/mech/mech.md",
        "Examples" =>  [
            "2D truss" => "examples/truss.md",
            "Solid beam" => "examples/solid-beam.md",
        ]
        # "Plotting" =>  "plotting/plotting.md",
    ],
    doctest = true,
    repo = "https://github.com/NumericalForge/Serendip.jl",
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