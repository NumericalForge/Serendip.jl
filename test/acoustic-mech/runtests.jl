using Serendip
include("../test-helpers.jl")

for file in [
    "acoustic-modal.jl",
    "acoustic-modal-free-surface.jl",
    "acoustic-modal-interface.jl",
    "acoustic-transient-interface.jl",
    "acoustic-transient-mech-compare.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
