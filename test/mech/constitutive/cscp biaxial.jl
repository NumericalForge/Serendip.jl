using Serendip
using Test

geo = GeoModel()

add_block(geo, [0,0,0], [0.1,0.1,0.1], shape=HEX8, nx=1, ny=1, nz=1, tag="solids")

mesh= Mesh(geo)

fc = -30.87e3
ft = 2.95e3


mapper = RegionMapper()
add_mapping(mapper, "solids", MechBulk, CSCP, E=30e6, nu=0.2, alpha=0.666, beta=1.15, fc=fc, ft=ft, epsc=0.002, GF=0.1, wc=0.00005)

model = FEModel(mesh, mapper)
ana = MechAnalysis(model)

add_logger(ana, :ip, [0.05, 0.05, 0.05], "cscp.dat")
# add_monitor(ana, :ip, [0.05, 0.05, 0.05], :(sxx, syy), stop=:( rho<0.3*rho_max))


θ = 225
α = θ

U      = 0.0014
ΔT     = 0.01
T      = 0.0
factor = 1.1
dT0    = 0.001

while T<1.0
    ux = cosd(α)*U*ΔT
    uy = sind(α)*U*ΔT
    # boundary conditions
    bcs = [
        x==0 => NodeBC(ux=0),
        y==0 => NodeBC(uy=0),
        z==0 => NodeBC(uz=0),
    ]

    if cosd(θ)!=0
        bcs = [ bcs
            x==0.1 => NodeBC(ux=ux)
        ]
    end
    if sind(θ)!=0
        bcs = [ bcs
            y==0.1 => NodeBC(uy=uy)
        ]
    end

    addstage!(ana, bcs)
    res = solve!(ana, autoinc=true, tol=0.01, rspan=0.02, dT0=dT0)
    global dT0 = 0.05

    sx = ana.loggers[1].table[:σxx][end]
    sy = ana.loggers[1].table[:σyy][end]
    β  = (atand(sy, sx) + 360) % 360

    global T = T+ΔT
    global ΔT *= factor

    global α = α + 2*(θ-β)

    @show θ
    @show β

    contains(res.message, "Stop") && break
    contains(res.message, "not converge") && break

end


include("plot_cscp.jl")