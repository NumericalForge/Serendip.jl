using Serendip
include("../test-helpers.jl")

for file in [
    # "elements/elastic-solid.jl",
    "elements/elastic-elems.jl",
    "elements/variable-thickness.jl",
    # "elements/cohesive-thickness.jl",
    # "elements/contact-thickness.jl",
    "elements/elastic-quad4.jl",
    "elements/elastic-hex8.jl",
    "elements/axisymmetric.jl",
    "elements/constraints.jl",
    "elements/rod/frame.jl",
    "elements/rod/truss.jl",
    "elements/beam/cantilever.jl",
    "elements/beam/simple-beam.jl",
    "elements/beam/curved-beam.jl",
    "elements/shell/shell.jl",
    "elements/line_interface/dowel.jl",
    "elements/line_interface/joint1d.jl",
    "elements/bar/embedded.jl",
    "elements/tip/tip.jl",
    "constitutive/vm-bar.jl",
    "constitutive/vm-2d.jl",
    "constitutive/vm-3d.jl",
    "constitutive/vm-beam-shell.jl",
    "constitutive/vm-beam-2d.jl",
    "constitutive/vm-beam-3d.jl",
    "constitutive/ceb.jl",
    "constitutive/dp.jl",
    "constitutive/escp-uniaxial-biaxial.jl",
    "constitutive/bondslip-monotonic.jl",
    "constitutive/cohesive.jl",
    "constitutive/contact.jl",
    # "constitutive/coulomb-contact-update.jl",
    "constitutive/coulomb-contact.jl",
]
    printstyled("\nRunning file ", file, "...\n", color=:yellow, bold=true)
    @testset "$file" begin
        include(file)
    end
end
