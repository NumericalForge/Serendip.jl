using Glob

files = [
    "truss.jl",
    "2d-static.jl",
    "3d-static.jl",
    "embedded-rebar.jl",
    "discrete-rebar.jl",
    "dynamic.jl",
]

for file in files
    printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
    include(file)
    println()
    rm.(glob("*.pdf"), force=true)
    rm.(glob("*.vtk"), force=true)
    rm.(glob("*.log"), force=true)
    rm("output", recursive=true, force=true)
end