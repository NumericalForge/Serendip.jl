using Test
using Serendip

limits = Serendip.compute_auto_limits([2.0, 2.0])
@test limits[1] < 2.0 < limits[2]
@test isapprox(limits[2] - 2.0, 2.0 - limits[1])

ax = Serendip.Axis(direction=:horizontal, limits=[1.0e4, 5.0e4], label="x")
Serendip.configure!(ax)
@test ax.tick_exponent == 4
@test !isempty(ax.exponent_box.text)
@test all(!occursin("10", lbl) for lbl in ax.tick_labels)

ay = Serendip.Axis(direction=:vertical, limits=[1.0e-5, 5.0e-5], label="y")
Serendip.configure!(ay)
@test ay.tick_exponent == -5
@test !isempty(ay.exponent_box.text)

ax_small = Serendip.Axis(direction=:horizontal, limits=[0.0, 3.1e-3], label="x", ticks=[0.0, 1.0e-3, 2.0e-3, 3.1e-3])
Serendip.configure!(ax_small)
@test ax_small.tick_exponent == -3
@test "3.1" in ax_small.tick_labels
@test !isempty(ax_small.exponent_box.text)

aint = Serendip.Axis(direction=:horizontal, limits=[1.0e4, 4.0e4], ticks=[1.0e4, 2.0e4, 3.0e4, 4.0e4])
Serendip.configure!(aint)
@test all(!occursin(".", lbl) for lbl in aint.tick_labels)

amanual = Serendip.Axis(direction=:horizontal, limits=[1.0e4, 3.0e4], ticks=[1.0e4, 2.0e4, 3.0e4], tick_labels=["A", "B", "C"])
Serendip.configure!(amanual)
@test amanual.tick_labels == ["A", "B", "C"]
@test isempty(amanual.exponent_box.text)

line_chart = Chart(quiet=true)
add_line(line_chart, [1.0e4, 2.0e4, 3.0e4], [1.0e-5, 2.0e-5, 3.0e-5]; label="scaled")
Serendip.configure!(line_chart)
@test line_chart.xaxis.tick_exponent == 4
@test line_chart.yaxis.tick_exponent == -5

bar_chart = Chart(quiet=true)
add_bar(bar_chart, [1.0, 2.0], [3.0, 4.0]; bar_base=-2.0, label="bars")
Serendip.configure!(bar_chart)
@test bar_chart.yaxis.limits[1] < -2.0
@test bar_chart.yaxis.limits[2] > 4.0

small_span = Serendip.Axis(direction=:horizontal, limits=[1.0e-10, 1.5e-10])
Serendip.configure!(small_span)
@test length(small_span.ticks) > 1

aw2d = Serendip.AxesWidget(labels=["x", "y"])
Serendip.configure!(aw2d)
@test aw2d.width > 0
@test aw2d.height > 0

aw3d = Serendip.AxesWidget(labels=["x", "y", "z"], azimuth=30.0, elevation=20.0, distance=10.0, up=:z)
Serendip.configure!(aw3d)
@test aw3d.width > 0
@test aw3d.height > 0

aw3d_x = Serendip.AxesWidget(labels=["x", "y", "z"], azimuth=30.0, elevation=20.0, distance=10.0, up=:x)
Serendip.configure!(aw3d_x)
@test any(i -> !isapprox(aw3d.projected_axes[i][1], aw3d_x.projected_axes[i][1]) || !isapprox(aw3d.projected_axes[i][2], aw3d_x.projected_axes[i][2]), eachindex(aw3d.projected_axes))

aw3d_rot = Serendip.AxesWidget(labels=["x", "y", "z"], azimuth=55.0, elevation=20.0, distance=10.0, up=:z)
Serendip.configure!(aw3d_rot)
@test any(i -> !isapprox(aw3d.projected_axes[i][1], aw3d_rot.projected_axes[i][1]) || !isapprox(aw3d.projected_axes[i][2], aw3d_rot.projected_axes[i][2]), eachindex(aw3d.projected_axes))
