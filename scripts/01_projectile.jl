using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
include(joinpath(@__DIR__, "..", "src", "mechanics", "projectile.jl"))

sol = simulate_projectile(v0=30.0, angle_deg=55.0, drag=0.02,
                          tspan=(0.0, 6.0))

ts = range(sol.t[1], sol.t[end]; length=240)
positions = [Point2f(sol(t)[1], sol(t)[2]) for t in ts]

xmin = minimum(p[1] for p in positions) - 2
xmax = maximum(p[1] for p in positions) + 2
ymin = 0.0
ymax = maximum(p[2] for p in positions) + 2

fig = Figure(size=(960, 540))
ax = Axis(fig[1, 1];
          xlabel="x (m)", ylabel="y (m)",
          title="Projectile motion with quadratic drag",
          limits=(xmin, xmax, ymin, ymax))

hlines!(ax, [0.0]; color=:gray, linewidth=2)

trail = Observable(Point2f[])
pt    = Observable(Point2f(positions[1]))

lines!(ax, trail; color=:dodgerblue, linewidth=3)
scatter!(ax, pt; color=:crimson, markersize=20)

outpath = joinpath(@__DIR__, "..", "output", "videos", "01_projectile.mp4")
record(fig, outpath, eachindex(ts); framerate=30) do i
    push!(trail[], positions[i]); notify(trail)
    pt[] = positions[i]
end

@info "wrote $outpath"
