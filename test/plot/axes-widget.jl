using Test
using Serendip


aw2d = Serendip.AxesWidget(labels=["x", "y"])
Serendip.configure!(aw2d)
@test aw2d.width > 0
@test aw2d.height > 0

aw3d = Serendip.AxesWidget(
    labels=["x", "y", "z"],
    azimuth=30.0,
    elevation=20.0,
    distance=10.0,
    up=:z,
)
Serendip.configure!(aw3d)
@test aw3d.width > 0
@test aw3d.height > 0

aw3d_x = Serendip.AxesWidget(
    labels=["x", "y", "z"],
    azimuth=30.0,
    elevation=20.0,
    distance=10.0,
    up=:x,
)
Serendip.configure!(aw3d_x)
@test any(
    i -> !isapprox(aw3d.projected_axes[i][1], aw3d_x.projected_axes[i][1]) ||
         !isapprox(aw3d.projected_axes[i][2], aw3d_x.projected_axes[i][2]),
    eachindex(aw3d.projected_axes),
)

aw3d_rot = Serendip.AxesWidget(
    labels=["x", "y", "z"],
    azimuth=55.0,
    elevation=20.0,
    distance=10.0,
    up=:z,
)
Serendip.configure!(aw3d_rot)
@test any(
    i -> !isapprox(aw3d.projected_axes[i][1], aw3d_rot.projected_axes[i][1]) ||
         !isapprox(aw3d.projected_axes[i][2], aw3d_rot.projected_axes[i][2]),
    eachindex(aw3d.projected_axes),
)

@test Serendip.effective_arrow_head_length(20.0, 7.0) == 7.0
@test isapprox(Serendip.effective_arrow_head_length(6.0, 7.0), 3.6)
