using Serendip
using Test
include("test-helpers.jl")

@testset verbose=true "Serendip" begin
    for file in [
        "mesh/runtests.jl",
        "model/runtests.jl",
        "plot/runtests.jl",
        "tools/runtests.jl",
        "mech/runtests.jl",
        "dynamic/runtests.jl",
        # "thermomech/runtests.jl",
        # "acoustic-mech/runtests.jl",
        # "hydromech/runtests.jl",
    ]
        printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
        @testset verbose=true "$file" begin
            include(file)
        end
    end
end

include("clean.jl")
