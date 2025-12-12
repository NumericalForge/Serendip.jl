using Serendip

X = collect(0:0.4:8)
Y1 = sin.(X)
Y2 = cos.(X)
Y3 = sin.(X) .* cos.(X)

chart = Chart(
    xlabel=L"$x$ coordinate",
    ylabel=L"$y$ coordinate",
    legend=:bottom_right
)

add_line(chart, X, Y1, mark=:circle, line_width=0.5, label=L"\sin(x)")
add_line(chart, X, Y2, mark=:utriangle, color=:royalblue, label=L"\cos(x)")
add_line(chart, X, Y3, mark=:square, color=:green, line_style=:dash, label=L"\sin(x)\, \cos(x)")

add_annotation(chart, Annotation(L"Serendip $\frac{x}{y}$", 0.5, 0.9, target=[1.5, 1], text_alignment=:top))
add_annotation(chart, Annotation(L"Serendip $\frac{x}{y}$", 0.7, 0.2, target=[4, 0], text_alignment=:top))


save(chart, "chart.pdf")