using Serendip
include("../test-helpers.jl")
for file in [
    "expr.jl",
    "numerical.jl",
    "table-frequency.jl",
    "xml.jl",
    "show.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
