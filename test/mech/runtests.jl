using Serendip

@run_files [
    # Elements

    # bulk
    "elements/elastic-solid.jl",
    "elements/elastic-elems.jl",
    "elements/elastic-quad4.jl",
    "elements/elastic-hex8.jl",
    "elements/axisymmetric.jl",

    # frame
    "elements/rod/frame.jl",

    # rod
    "elements/rod/truss.jl",
    
    # beam
    "elements/beam/cantilever.jl",
    
    # shell
    "elements/shell/shell.jl",
    
    # joint
    "elements/line_interface/dowel.jl",
    "elements/line_interface/joint1d.jl",
    
    # embeddeed and semi-embedded
    "elements/bar/embedded.jl",
    "elements/tip/tip.jl",
    
    # Constitutive models
    "constitutive/vm-bar.jl",
    "constitutive/vm-2d.jl",
    "constitutive/vm-3d.jl",
    "constitutive/vm-beam-shell.jl",
    "constitutive/vm-beam-2d.jl",
    "constitutive/vm-beam-3d.jl",
    "constitutive/ceb.jl",
    "constitutive/dp.jl",

    # "elements/joint/joint2d.jl",
    "constitutive/cohesive.jl",
    "constitutive/contact.jl",
]