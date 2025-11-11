# From the book:
# Finite Elements for Structural Analysis. pg 168
# William Weaver & Paul Johnston

using Serendip
using Test

# Mesh generation

geo = GeoModel()
add_block(geo, [0, 0], 1, 1, 0, nx=1, ny=1, shape=QUAD4)
mesh = Mesh(geo, quiet=true)

mapper = RegionModel(MechBulk, LinearElastic, E=1.0, nu=0.25)

model = FEModel(mesh, mapper, stress_state=:plane_strain)
ana = MechAnalysis(model)

stage = add_stage(ana, nincs=1)
add_bc(stage, :face, x==0, ux=0)
add_bc(stage, :face, y==0, uy=0)
add_bc(stage, :face, y==1, ty=-1)

run(ana)


dis =
    [
        0.0       0.0000e+00
        0.3125    0.0000e+00
        0.0      -9.3750e-01
        0.3125   -9.3750e-01
    ]


println("Displacements:")

D = get_values(model.nodes)[[:ux, :uy]]
println(D)

@test dis ≈ Array(D) atol=1e-5

println("Stress:")
S = get_values(model.elems[1])[[:σxx, :σyy, :σxy]]
println(S)

println("Support reactions:")
F = get_values(model.nodes)[[:fx, :fy]]
println(F)