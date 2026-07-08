using Serendip
include("../test-helpers.jl")
for file in [
    "axes-widget.jl",
    "domain-colorbar.jl",
    "domain-layers.jl",
    "domain-markers.jl",
    "domain-selectors.jl",
    "domain-view.jl",
    "video.jl",
    "units.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
