# From the book:
# Finite Elements for Structural Analysis. pg 168
# William Weaver & Paul Johnston

using Serendip
using Test

# Mesh generation

geo = GeoModel()
add_block(geo, [0, 0, 0], 1, 1, 1, nx=1, ny=1, nz=1, shape=HEX8)
mesh = Mesh(geo, quiet=true)

# Model definition
mapper = RegionModel(MechBulk, LinearElastic, E=1.0, nu=0.3)

# Load cases

bcs1 = [
    [(x==0, y==0, z==0), :node, (ux=0, uy=0)],
    [(x==1, y==0, z==0), :node, (uy=0,)],
    [(x==0, y==1, z==0), :node, (ux=0,)],
    [(z==0)            , :node, (uz=0,)],
    [(z==1)            , :node, (fz=1,)],
]

bcs2 = [
    [(x==0, y==0, z==0), :node, (ux=0, uy=0)],
    [(x==1, y==0, z==0), :node, (uy=0,)],
    [(x==0, y==1, z==0), :node, (ux=0,)],
    [(z==0)            , :node, (uz=0,)],
    [(y==1, z==1)      , :edge, (qy=2,)],
]

bcs3 = [
    [(x==0, y==0, z==0), :node, (ux=0, uy=0)],
    [(x==1, y==0, z==0), :node, (uy=0,)],
    [(x==0, y==1, z==0), :node, (ux=0,)],
    [(z==0)            , :node, (uz=0,)],
    [(x==1)            , :face, (tx=3*z,)],
]

bcs4 = [
    [(x==0, y==0, z==0), :node, (ux=0, uy=0)],
    [(x==1, y==0, z==0), :node, (uy=0,)],
    [(x==0, y==1, z==0), :node, (ux=0,)],
    [(z==0)            , :node, (uz=0,)],
    [(x>=0)            , :body, (wz=-1,)],
]

lc_list = ["Nodal load", "Edge load", "Triangular face load", "Volume load"]
bcs_list = [bcs1, bcs2, bcs3, bcs4]
dis_list = [
            [ 0.0, 0.0, 0.0, 0.0, 4.0    ,    4.0,     4.0,      4.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 3.32088, 3.1998, -4.6002, -4.32047 ],
            [ 0.0, 0.0, 0.0, 0.0, 1.51044,-2.4501,  1.4499, -2.31023 ],
            [ 0.0, 0.0, 0.0, 0.0, -0.5   ,   -0.5,    -0.5,     -0.5 ] ]

for (lc, bcs, dis) in zip(lc_list, bcs_list, dis_list)

    println("\nLoad case: $lc \n")

    global model = FEModel(mesh, mapper)
    ana = MechAnalysis(model)
    stage = add_stage(ana, nouts=1)
    for bc in bcs
        selector, kind, conds = bc
        add_bc(stage, kind, selector; conds...)
    end

    run(ana)

    println("Displacements:")
    D = get_values(model.nodes)[[:ux, :uy, :uz]]
    println(D)

    println(@test dis ≈ D[:uz] atol=1e-5)

    println("Stress:")
    S = get_values(model.elems[1])[[:σxx, :σyy, :σzz, :σyz, :σxz, :σxy]]
    println(S)

    println("Support reactions:")
    # F = get_values(model.nodes[:(z==0)])[[:fx, :fy, :fz]]
    F = get_values(select(model, :node, z==0))[[:fx, :fy, :fz]]
    println(F)
end


