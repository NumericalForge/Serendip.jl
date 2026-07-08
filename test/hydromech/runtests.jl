using Serendip
include("../test-helpers.jl")

for file in [
    "seep.jl",
    "cutoff.jl",
    "terzaghi.jl",
    "terzaghi-joint.jl",
    "drain.jl",
    "drain-solid.jl",
    "hm-drain.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
