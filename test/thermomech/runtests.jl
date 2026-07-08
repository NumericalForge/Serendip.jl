using Serendip
include("../test-helpers.jl")

for file in [
    "thermo.jl",
    "thermomech.jl",
    "tm-shell.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
