using Serendip

function mirrored_curve(σx::Vector{Float64}, σy::Vector{Float64})
    x = Float64[]
    y = Float64[]

    for i in length(σx):-1:1
        i == 1 && isapprox(σx[i], σy[i]; atol=1e-8, rtol=1e-8) && continue
        push!(x, σy[i])
        push!(y, σx[i])
    end

    append!(x, σx)
    append!(y, σy)
    return x, y
end


envelope = DataTable(joinpath(@__DIR__, "escp-biaxial-envelope.dat"))

f0 = envelope["f0"][1]
σx = envelope["σxx"] ./ f0
σy = envelope["σyy"] ./ f0

curve_x, curve_y = mirrored_curve(collect(σx), collect(σy))

diag = range(-1.15, 0.15, 100)

chart = Chart(
    xlabel=t"$macron(σ)_1 / f_0$",
    ylabel=t"$macron(σ)_2 / f_0$",
    legend=:bottom_right,
    xlimits=[-1.4, 0.2],
    ylimits=[-1.4, 0.2],
    aspect_ratio=:equal,
)
add_line(chart, diag, diag, color=:gray, line_style=:dash, label="symmetry line")
add_line(chart, curve_x, curve_y, color=:black, label="ESCP envelope")
add_scatter(chart, σx, σy, color=:red, mark=:square, label="evaluated angles")
add_scatter(chart, σy[2:end], σx[2:end], color=:red, mark=:square, label="")
save(chart, joinpath(@__DIR__, "escp-biaxial-stress-plane.pdf"))
