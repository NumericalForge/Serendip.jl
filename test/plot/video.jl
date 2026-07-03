using Serendip
using Test

function limits_changed(a, b; atol=1.0e-8)
    return !all(isapprox.(a, b; atol=atol))
end

function warped_mesh_2d()
    geo = GeoModel(quiet=true)
    add_block(geo, [0, 0, 0], 1, 1, 0, nx=1, ny=1, shape=:quad4, tag="solids")
    mesh = Mesh(geo, quiet=true)
    mesh.node_fields["U"] = hcat(
        [0.4 * node.coord[1] + 0.1 for node in mesh.nodes],
        [0.2 * node.coord[2] for node in mesh.nodes],
        zeros(length(mesh.nodes)),
    )
    return mesh
end

function warped_mesh_3d()
    geo = GeoModel(quiet=true)
    add_block(geo, [0, 0, 0], 1, 1, 1, nx=1, ny=1, nz=1, tag="solids")
    mesh = Mesh(geo, quiet=true)
    mesh.node_fields["U"] = hcat(
        [0.2 * node.coord[1] + 0.05 for node in mesh.nodes],
        [0.15 * node.coord[2] for node in mesh.nodes],
        [0.25 * node.coord[3] + 0.05 for node in mesh.nodes],
    )
    return mesh
end

frame_root = joinpath("output", "domain-video-frames")
ispath(frame_root) && rm(frame_root; recursive=true, force=true)

mesh2d = warped_mesh_2d()
plot1_2d = DomainPlot(mesh2d; face_color=:aliceblue, title="2D frame 1", warp=0.5)
plot2_2d = DomainPlot(mesh2d; face_color=:lightblue, title="2D frame 2", warp=3.0)

video2d = VideoBuilder(framerate=2, cleanup=true, tempdir=frame_root, freeze_scale=true)
add_frame(video2d, plot1_2d)
add_frame(video2d, plot2_2d)

outfile2d = joinpath("output", "domain-video-2d.mp4")
save(video2d, outfile2d)

@test isfile(outfile2d)
@test filesize(outfile2d) > 0
@test plot2_2d.frozen_scaling_state isa Serendip.DomainPlotScalingState2D
@test all(isapprox.(plot1_2d.canvas.limits, plot2_2d.canvas.limits; atol=1.0e-8))

plot3_2d = DomainPlot(mesh2d; face_color=:old_paper, title="2D frame 3", warp=0.5)
plot4_2d = DomainPlot(mesh2d; face_color=:linen, title="2D frame 4", warp=3.0)
video2d_pad = VideoBuilder(framerate=2, cleanup=true, tempdir=frame_root, freeze_scale=true, bounds_factor=1.1)
add_frame(video2d_pad, plot3_2d)
add_frame(video2d_pad, plot4_2d)
outfile2d_pad = joinpath("output", "domain-video-2d-pad.mp4")
save(video2d_pad, outfile2d_pad)

@test isfile(outfile2d_pad)
@test filesize(outfile2d_pad) > 0
@test plot3_2d.frozen_scaling_state isa Serendip.DomainPlotScalingState2D
@test plot4_2d.frozen_scaling_state isa Serendip.DomainPlotScalingState2D
@test all(isapprox.(plot3_2d.canvas.limits, plot3_2d.frozen_scaling_state.canvas_limits; atol=1.0e-8))
@test all(isapprox.(plot3_2d.canvas.limits, plot4_2d.canvas.limits; atol=1.0e-8))
@test all(isapprox.(plot4_2d.canvas.limits, plot4_2d.frozen_scaling_state.canvas_limits; atol=1.0e-8))
@test (plot3_2d.canvas.limits[3] - plot3_2d.canvas.limits[1]) > (plot1_2d.canvas.limits[3] - plot1_2d.canvas.limits[1])

plot5_2d = DomainPlot(mesh2d; title="2D free 1", warp=0.5)
plot6_2d = DomainPlot(mesh2d; title="2D free 2", warp=3.0)
video2d_free = VideoBuilder(framerate=2, cleanup=true, tempdir=frame_root, freeze_scale=false)
add_frame(video2d_free, plot5_2d)
add_frame(video2d_free, plot6_2d)
outfile2d_free = joinpath("output", "domain-video-2d-free.mp4")
save(video2d_free, outfile2d_free)

@test isfile(outfile2d_free)
@test filesize(outfile2d_free) > 0
@test plot6_2d.frozen_scaling_state === nothing
@test limits_changed(plot5_2d.canvas.limits, plot6_2d.canvas.limits)

mesh3d = warped_mesh_3d()
plot1_3d = DomainPlot(mesh3d; face_color=:aliceblue, title="3D frame 1", warp=0.5, azimuth=35, elevation=20)
plot2_3d = DomainPlot(mesh3d; face_color=:lightblue, title="3D frame 2", warp=3.0, azimuth=35, elevation=20)

video3d = VideoBuilder(framerate=2, cleanup=true, tempdir=frame_root, freeze_scale=true)
add_frame(video3d, plot1_3d)
add_frame(video3d, plot2_3d)

outfile3d = joinpath("output", "domain-video-3d.mp4")
save(video3d, outfile3d)

@test isfile(outfile3d)
@test filesize(outfile3d) > 0
@test plot2_3d.frozen_scaling_state isa Serendip.DomainPlotScalingState3D
@test all(isapprox.(plot1_3d.canvas.limits, plot2_3d.canvas.limits; atol=1.0e-8))

Serendip._apply_scaling_state!(plot2_3d, plot2_2d.frozen_scaling_state)
@test_throws Serendip.SerendipException Serendip.configure!(plot2_3d)
