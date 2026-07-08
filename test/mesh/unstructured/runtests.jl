using Serendip
include("../../test-helpers.jl")

for file in [
    "select.jl",
    "loops.jl",
    "hole.jl",
    "square.jl",
    "pull.jl",
    "cube.jl",
    "step-import.jl",
    "revolve.jl",
    "circular-wedge.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
