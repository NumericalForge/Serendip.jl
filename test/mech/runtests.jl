using Serendip

@run_files [
    # Static analysis: bulk elements
    "elem/elastic-solid.jl",
    "elem/elastic-elems.jl",
    "elem/elastic-quad4.jl",
    "elem/elastic-hex8.jl",
    "elem/axisymmetric.jl",

    # Static frame element
    "elem/rod/frame.jl",

    # Static analysis: rod elems
    "elem/rod/truss.jl",
    "elem/rod/pprod.jl",

    # Static analysis: beam elems
    "elem/beam/beam.jl",

    # Static analysis: shell elems
    "elem/shell/shell.jl",

    # Static analysis: joint models
    "elem/joint/joint1d.jl",

    # Static analysis: embeddeed and semi-embedded elements
    "elem/inset/embedded.jl",
    "elem/inset/ceb.jl",
    "elem/inset/tip.jl",

    # Constitutive models
    "mat/dp.jl",
    "mat/vm-2d.jl",
    "mat/vm-3d.jl",
    "mat/vm-beam-shell.jl",
    "mat/vm-beam-2d.jl",
    "mat/vm-beam-3d.jl",

    "elem/joint/joint2d.jl",
    "mat/power-yield-crack.jl",
]

