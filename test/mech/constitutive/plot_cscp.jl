using Serendip

table = DataTable("cscp.dat")


chart = Chart(
    xlabel = L"\xi",
    ylabel = L"\rho",
    xmult = 1e-3,
)

rθ = table.r[2]
α = 1.55

ξi = table.xi_i[1]
f̄c = table.fcb[1]

ξa = table.xi_a[1]
ξb = table.xi_b[1]
k = table.kappa[1]
rξ(ξ) = Serendip.spow((ξb-ξ)/f̄c, 1/α)

function rc(ξ::Float64)
    ξ>=ξi && return 1.0
    ξ<ξa  && return 0.0
    return √(1 - ((ξi-ξ)/(ξi-ξa))^2)
end

ξ0 = range(ξa, ξb, 400)
ρ0 = [ rθ*rc(ξ)*rξ(ξ)*k for ξ in ξ0 ]

ξa = table.xi_a[end]
ξb = table.xi_b[end]
k = table.kappa[end]


ξ1 = range(ξa, ξb, 400)
ρ1 = [ rθ*rc(ξ)*rξ(ξ)*k for ξ in ξ1]


series = [
    DataSeries(ξ0, ρ0, color=:red),
    DataSeries(ξ1, ρ1, color=:blue),
    DataSeries(table.xi, table.rho, marker=:circle),
]

add_series(chart, series)
save(chart, "chart.pdf")

# chart2

chart = Chart(
    xlabel = L"\varepsilon_{yy}\times 1000",
    ylabel = L"\sigma_{yy}",
    xmult = 1e3,
)


series = [
    DataSeries(table.eyy, table.syy, marker=:circle),
]

add_series(chart, series)
save(chart, "chart2.pdf")

# chart3

chart = Chart(
    xlabel = L"\varepsilon_{xx}\times 1000",
    ylabel = L"\sigma_{xx}",
    xmult = 1e3,
)


series = [
    DataSeries(table.eyy, table.sxx, marker=:circle),
]

add_series(chart, series)
save(chart, "chart3.pdf")