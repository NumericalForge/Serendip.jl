using Serendip


X = collect(0:0.4:8)
Y1 = sin.(X)
Y2 = cos.(X)
Y3 = sin.(X) .* cos.(X)

chart = Chart(
    title="Trigonometric Responses",
    background=:white,
    xlabel="`x` coordinate [mm]",
    ylabel="`y` coordinate [mm]",
    legend=:bottom_right
)

add_line(chart, X, Y1, mark=:circle, line_width=0.5, label="`sin(x)`")
add_line(chart, X, Y2, mark=:utriangle, color=:royalblue, label="`cos(x)`")
add_scatter(chart, X, Y3, mark=:square, color=:green, line_style=:dash, label="`sin(x) cos(x)`")

add_annotation(chart, Annotation(" `-2 + bold(A)^(1/2) - (A^2/B_3^2 + x^2/y^2 )^2 + f_2^2/g_2^2`", 0.3, 0.9, target=[1.5, 1], alignment=:top))
add_annotation(chart, Annotation("A `+ sigma`", 0.1, 0.2, target=[4.0, 0], alignment=:top))
add_annotation(chart, Annotation("`A_2^2 F_2^2 T_2^2`", 0.15, 0.3))
add_annotation(chart, Annotation("`a_2^2 f_2^2 a^2 f t_2^2`", 0.15, 0.4))
add_annotation(chart, Annotation("Legend overlap", 0.88, 0.08, alignment=:right))

save(chart, "chart.pdf")
@test isfile("chart.pdf")