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
