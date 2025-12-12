# Remove output files from tests

path  = dirname(@__FILE__)

println("  deleting temporary files...")
for (root, dirs, files) in walkdir(path)
    for file in files
        ext = split(file*".", ".")[2]
        if contains(ext, r"step|vtk|vtu|pdf|dat|table|book|log")
            fullname = joinpath(root, file)
            rm(fullname)
        end
    end
end
