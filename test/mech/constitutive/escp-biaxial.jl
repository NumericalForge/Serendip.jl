using Serendip
using Test


function run_target_angle(mesh::Mesh, mapper::RegionMapper, θ::Float64, ℓ::Float64, h::Float64)
    model = FEModel(mesh, mapper, quiet=true)
    ana   = MechAnalysis(model, outkey="escp-biaxial", outdir=@__DIR__)
    ip    = nearest(select(model, :ip, :all), [ℓ / 2, ℓ / 2, h / 2])
    log   = add_logger(ana, :ip, collect(ip.coord))
    add_monitor(ana, :ip, ip.coord, stop=:(ρ < 0.95 * ρ_max))

    α = θ
    U = 0.006
    ΔT = 0.01
    T = 0.0
    factor = 1.05
    dT0 = 0.001
    Tmax = 1.5

    while T < Tmax
        ux = cosd(α) * U * ΔT
        uy = sind(α) * U * ΔT

        stage = add_stage(ana)
        add_bc(stage, :node, x == 0, ux=0)
        add_bc(stage, :node, y == 0, uy=0)
        add_bc(stage, :node, z == 0, uz=0)

        abs(cosd(θ)) > 1e-12 && add_bc(stage, :node, x == ℓ, ux=ux)
        abs(sind(θ)) > 1e-12 && add_bc(stage, :node, y == ℓ, uy=uy)

        dTmax = θ<150 ? 0.005 : 0.1

        status = run(ana, autoinc=true, tol=0.01, rspan=0.01, dT0=dT0, dTmax=dTmax, quiet=true)
        if !status.successful
            @test occursin("monitor condition", lowercase(status.message))
            break
        end

        table = log.table
        σxx = table["σxx"][end]
        σyy = table["σyy"][end]
        β = mod(atand(σyy, σxx), 360.0)

        T += ΔT
        ΔT *= factor
        α += 2 * (θ - β)
        dT0 = 0.01
    end

    table = log.table
    ipeak = argmax(table["ρ"])
    return (
        table = table,
        ipeak = ipeak,
    )
end


function to_table(cols::Pair...)
    header = [string(k) for (k, _) in cols]
    data = AbstractVector[v for (_, v) in cols]
    return DataTable(header, data)
end


function main()
    ℓ = 0.1
    h = 0.1
    fc = -30.0e3
    ft = 3.0e3
    f0 = abs(fc)

    # Mesh
    geo = GeoModel()
    add_block(geo, [0, 0, 0], ℓ, ℓ, h, shape=:hex8, nx=1, ny=1, nz=1, tag="solids")
    mesh = Mesh(geo, quadratic=true)

    # Mapper
    mapper = RegionMapper()
    add_mapping( mapper, "solids", MechSolid, ESCP, E=30.0e6, nu=0.25, beta=1.15, fc=fc, epsc=-0.002, ft=ft, wc=0.0001 )

    # One symmetric branch is enough; the plot script mirrors it across σxx = σyy.
    angles = collect(45.0:15.0:225.0)
    angles = [ 45.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 175.0, 180.0, 188.0, 196.0, 205.0, 215.0, 225.0 ]

    θ_env = Float64[]
    σxx_env = Float64[]
    σyy_env = Float64[]
    f0_env = Float64[]

    for θ in angles
        printstyled("Running target angle $(θ)°...\n", color=:yellow)
        result = run_target_angle(mesh, mapper, θ, ℓ, h)
        table = result.table
        ipeak = result.ipeak

        n = length(table["σxx"])
        @test n > 5

        push!(θ_env, θ)
        σxx_peak = table["σxx"][ipeak]
        σyy_peak = table["σyy"][ipeak]
        push!(σxx_env, σxx_peak)
        push!(σyy_env, σyy_peak)
        @show σxx_peak, σyy_peak
        push!(f0_env, f0)
    end

    envelope = to_table(
        :θ => θ_env,
        :σxx => σxx_env,
        :σyy => σyy_env,
        :f0 => f0_env,
    )

    save(envelope, joinpath(@__DIR__, "escp-biaxial-envelope.dat"))

    include("plot_escp_biaxial.jl")
end


main()
