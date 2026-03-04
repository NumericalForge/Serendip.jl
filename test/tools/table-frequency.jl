using Serendip
using Test

f0 = 7.5
t = collect(range(0.0, 1.2, length=1601))
u = 0.2*sin.(2*pi*f0.*t) .+ 0.01*cos.(2*pi*2.0.*t)

table1 = DataTable(t=t, uz=u)
f1 = get_frequency(table1, :uz)
@test isapprox(f1, f0; rtol=0.02)

table2 = DataTable(["t", "uz"], [t, u])
f2 = get_frequency(table2, "uz")
@test isapprox(f2, f0; rtol=0.02)

f3 = get_frequency(table2, "uz", extrema=:maxima)
@test isapprox(f3, f0; rtol=0.03)
