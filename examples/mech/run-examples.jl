using Glob

files = [
    "truss-2d.jl",
    "static-2d.jl",
    "static-3d.jl",
    "embedded-rebar.jl",
    "reinforced-beam.jl",
    "composite-beam.jl",
    "beam-crack-2d.jl",
    "beam-crack-3d.jl",
    # "dynamic.jl",
]

for file in files
    printstyled("\nRunning file ", file,"...\n", color=:yellow, bold=true)
    include(file)
    println()
    rm.(glob("*.step"), force=true)
    rm.(glob("*.vt?"), force=true)
    rm.(glob("*.log"), force=true)
    rm.(glob("*.pdf"), force=true)
    rm("output", recursive=true, force=true)
end