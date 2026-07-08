using Serendip
include("../test-helpers.jl")

for file in [
    "identity.jl",
    "select.jl",
    "shape/shape_deriv.jl",
    "shape/extrapolation.jl",
    "structured.jl",
    "path-tips.jl",
    "path-modes.jl",
    "cohesive-insertion.jl",
    "join-meshes.jl",
    "geo-add-mesh.jl",
    "unstructured/runtests.jl",
    "io.jl",
    # "operations.jl",
    "extrude.jl",
    "smoothing.jl",
    "revolve.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
