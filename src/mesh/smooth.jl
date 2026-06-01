# This file is part of Serendip package. See copyright license in https://github.com/NumericalForge/Serendip.jl


function faces_normal(faces::Vector{Cell}, facetol)
    ndim = 1 + faces[1].shape.ndim
    normals = Vector{Float64}[]

    for face in faces
        C = get_coords(face, ndim)

        # Shift coordinates to avoid singular regressions when the best-fit
        # line/plane crosses the origin.
        if ndim == 2
            C .+= [pi pi^1.1]
        else
            C .+= [pi pi^1.1 pi^1.2]
        end

        I = ones(size(C, 1))
        N = pinv(C) * I
        normalize!(N)

        if all(norm(N - NN) > facetol for NN in normals)
            push!(normals, N)
        end
    end

    return normals
end


mutable struct sNode
    node::Node
    faces::Array{Cell}
    normals
end


function str_histogram(hist::Vector{Int64})
    m = maximum(hist)
    H = m == 0 ? zeros(Int, length(hist)) : round.(Int, hist ./ m * 7)
    chars = [" ", "_", "▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"]
    return "[" * join(hist[i] == 0 ? " " : chars[H[i] + 2] for i in 1:length(H)) * "]"
end


include("smooth-lap.jl")
include("smooth-def.jl")


"""
    smooth(mesh; algorithm=:laplacian,
           maxit = algorithm == :deformation ? 10 : 20,
           quiet=true,
           fixed_boundary=false,
           mintol = algorithm == :deformation ? 2e-2 : 1e-2,
           tol = algorithm == :deformation ? 1e-3 : 1e-4,
           facetol = algorithm == :deformation ? 1e-4 : 1e-5,
           binsize=0.05,
           smart=false,
           weighted=false,
           alpha=1.0,
           extended=false,
           conds=nothing)

Smooth a finite element mesh in place and return the mutated `mesh`.

The `algorithm` keyword selects the smoothing method. It must be a `Symbol`.
Supported values are:

- `:laplacian`: Laplacian smoothing, with optional smart and weighted updates.
- `:deformation`: Physics-based deformation smoothing with Lagrange multiplier
  boundary constraints.

Keyword arguments used by both algorithms are:

- `maxit`: maximum number of iterations. Defaults to `20` for `:laplacian`
  and `10` for `:deformation`.
- `quiet`: suppress progress output when `true`. Defaults to `true`.
- `fixed_boundary`: keep boundary nodes fixed when `true`. Defaults to `false`.
- `mintol`: minimum-quality convergence tolerance. Defaults to `1e-2` for
  `:laplacian` and `2e-2` for `:deformation`.
- `tol`: average-quality convergence tolerance. Defaults to `1e-4` for
  `:laplacian` and `1e-3` for `:deformation`.
- `facetol`: tolerance used to identify boundary normals. Defaults to `1e-5`
  for `:laplacian` and `1e-4` for `:deformation`.
- `binsize`: quality histogram bin size used in progress output. Defaults to
  `0.05`.
- `smart`: reject local moves that reduce element quality when possible.
  Defaults to `false`.

Algorithm-specific keyword arguments are:

- `weighted`: only for `:laplacian`; use centroidal weighted Laplacian updates.
  Defaults to `false`.
- `alpha`: only for `:deformation`; force scaling parameter. Defaults to `1.0`.
- `extended`: only for `:deformation`; use quality-weighted force scaling.
  Defaults to `false`.
- `conds`: only for `:deformation`; additional boundary-condition expressions.
  Defaults to `nothing`.

# Examples

```julia
smooth(mesh)
smooth(mesh; algorithm = :laplacian, smart = true)
smooth(mesh; algorithm = :deformation, fixed_boundary = true, maxit = 5)
```
"""
function smooth(
    mesh::Mesh;
    algorithm::Symbol = :laplacian,
    maxit::Int64 = algorithm == :deformation ? 10 : 20,
    quiet = true,
    fixed_boundary = false,
    mintol::Float64 = algorithm == :deformation ? 2e-2 : 1e-2,
    tol::Float64 = algorithm == :deformation ? 1e-3 : 1e-4,
    facetol::Float64 = algorithm == :deformation ? 1e-4 : 1e-5,
    binsize::Float64 = 0.05,
    smart = false,
    weighted = false,
    alpha::Float64 = 1.0,
    extended = false,
    conds = nothing,
)
    if algorithm == :laplacian
        alpha == 1.0 || error("smooth: alpha only applies to algorithm=:deformation")
        extended == false || error("smooth: extended only applies to algorithm=:deformation")
        conds === nothing || error("smooth: conds only applies to algorithm=:deformation")
        return laplacian_smooth(
            mesh;
            maxit,
            quiet,
            fixed_boundary,
            mintol,
            tol,
            facetol,
            binsize,
            smart,
            weighted,
        )
    elseif algorithm == :deformation
        weighted == false || error("smooth: weighted only applies to algorithm=:laplacian")
        return deformation_smooth(
            mesh;
            maxit,
            quiet,
            fixed_boundary,
            mintol,
            tol,
            facetol,
            binsize,
            smart,
            alpha,
            extended,
            conds,
        )
    else
        error("smooth: unknown algorithm $(repr(algorithm)). Expected :laplacian or :deformation")
    end
end
