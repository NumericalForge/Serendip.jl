using Serendip

X = collect(1:5)
Y = [1.5, 2.2, 0.1, 1.1, 2.8]

chart = Chart(
    xlabel=t"Category",
    ylabel=t"Value",
    legend=:top_left
)

add_bar(chart, X, Y;
    color=:steelblue,
    label=t"bar series"
)

save(chart, "bar-chart.pdf")
