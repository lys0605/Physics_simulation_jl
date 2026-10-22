using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
include(joinpath(@__DIR__, "..", "src", "mechanics", "cyclic_pursuit.jl"))
include(joinpath(@__DIR__, "..", "src", "plot_utils.jl"))


φ = range(0, 4π; length=400)

fig, ax = plot_function1d(radial_path_angle, φ)

save_path = joinpath(@__DIR__, "..", "output", "plots", "05_cyclic_pursuit_path_r_against_phi_n=$N")
savefig()

function plot_cyclic_pursuit_path_polar(;n=4, a=1)
    fig  = Figure()
    ax   = Axis(fig[1,1]; xlabel=L"x", ylabel=L"y", aspect=DataAspect())
    initial_angles = range(0, 2π/n*(n-1); length=n)
    initial_points_xs = Float64[]
    initial_points_ys = Float64[]
    for φ₀ in initial_angles
        φ = range(0+φ₀, 4π+φ₀; length=400)
        xs = radial_path_angle.(φ ;n=n, a=a, φ₀=φ₀) .* cos.(φ)
        ys = radial_path_angle.(φ ;n=n, a=a, φ₀=φ₀) .* sin.(φ)
        push!(initial_points_xs, xs[1])
        push!(initial_points_ys, ys[1])
        lines!(ax, xs, ys)
    end
    push!(initial_points_xs, initial_points_xs[1])
    push!(initial_points_ys, initial_points_ys[1])
    lines!(ax, initial_points_xs, initial_points_ys; color=:black, linewidth=1.5)
    return (; fig, ax)
end

plot = plot_cyclic_pursuit_path_polar(n=8)
plot.fig


