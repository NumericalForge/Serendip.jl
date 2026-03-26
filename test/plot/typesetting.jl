using Test
using Serendip

w1, h1 = Serendip.getsize(raw"Axial $sigma_n$", 8.0)
w2, h2 = Serendip.getsize("Axial `sigma_n`", 8.0)
@test isapprox(w1, w2; atol=1.0e-6)
@test isapprox(h1, h2; atol=1.0e-6)

wm, hm = Serendip.getsize(raw"Load `P` at $x$", 8.0)
@test wm > 0
@test hm > 0


chart = Chart(
    title="Axial `sigma_n`",
    xlabel="`x`",
    ylabel="`u_x`",
    quiet=true,
)
add_line(chart, [0.0, 1.0], [0.0, 1.0]; label="`u_x`")
add_annotation(chart, Annotation(raw"Load `P` at $x$", 0.2, 0.2))
save(chart, "output/typesetting-backtick.pdf")
@test isfile("output/typesetting-backtick.pdf")
