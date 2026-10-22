using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
include(joinpath(@__DIR__, "..", "src", "mechanics", "projectile.jl"))
include(joinpath(@__DIR__, "..", "src", "plot_utils.jl"))

sol = simulate_projectile(v0=30.0, angle_deg=55.0, drag=0.02, tspan=(0.0, 6.0))

ts        = range(sol.t[1], sol.t[end]; length=240)
positions = [Point2f(sol(t)[1], sol(t)[2]) for t in ts]

ymax    = maximum(p[2] for p in positions) + 2
outpath = joinpath(@__DIR__, "..", "output", "videos", "01_projectile.mp4")
animate_trail(positions, outpath;
              title  = "Projectile motion with quadratic drag",
              xlabel = "x (m)", ylabel = "y (m)",
              ylims  = (0.0, ymax))

@info "wrote $outpath"
