
include("../test-helpers.jl")
for file in [
    "xml-io.jl",
    "logger.jl",
    "monitor.jl",
    "quadrature-mapping.jl",
    "select-points.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
