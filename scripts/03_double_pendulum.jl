using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
include(joinpath(@__DIR__, "..", "src", "mechanics", "double_pendulum.jl"))

L1, L2 = 1.0, 1.0
sol = simulate_double_pendulum(θ1₀=2.2, θ2₀=2.6, L1=L1, L2=L2,
                               tspan=(0.0, 15.0))

ts = range(sol.t[1], sol.t[end]; length=450)

R = L1 + L2 + 0.2
fig = Figure(size=(720, 720))
ax = Axis(fig[1, 1];
          title="Double pendulum (chaotic)",
          limits=(-R, R, -R, R), aspect=DataAspect())
hidedecorations!(ax); hidespines!(ax)

rods = Observable([Point2f(0, 0), Point2f(0, 0), Point2f(0, 0)])
trail_pts = Point2f[]
trail = Observable(Point2f[])
bobs = Observable([Point2f(0, 0), Point2f(0, 0)])

lines!(ax, trail; color=:crimson, linewidth=1.5)
lines!(ax, rods; color=:black, linewidth=3)
scatter!(ax, bobs; color=[:steelblue, :crimson], markersize=[22, 22])

outpath = joinpath(@__DIR__, "..", "output", "videos", "03_double_pendulum.mp4")
record(fig, outpath, eachindex(ts); framerate=30) do i
    θ1, _, θ2, _ = sol(ts[i])
    p1 = Point2f(L1 * sin(θ1), -L1 * cos(θ1))
    p2 = p1 + Point2f(L2 * sin(θ2), -L2 * cos(θ2))
    rods[] = [Point2f(0, 0), p1, p2]
    bobs[] = [p1, p2]
    push!(trail_pts, p2)
    length(trail_pts) > 250 && popfirst!(trail_pts)
    trail[] = copy(trail_pts)
end

@info "wrote $outpath"
