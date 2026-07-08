using Serendip
include("../test-helpers.jl")

for file in [
    "cantilever-beam.jl",
    "cantilever-solid.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
