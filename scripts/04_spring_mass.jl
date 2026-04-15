using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
include(joinpath(@__DIR__, "..", "src", "mechanics", "spring_mass.jl"))

sol = simulate_spring_mass(x0=1.5, v0=0.0, m=1.0, k=4.0, c=0.15,
                           tspan=(0.0, 20.0))

ts = range(sol.t[1], sol.t[end]; length=400)
xs = [sol(t)[1] for t in ts]
vs = [sol(t)[2] for t in ts]

fig = Figure(size=(1200, 540))
ax1 = Axis(fig[1, 1];
           xlabel="position", ylabel="",
           title="Damped harmonic oscillator",
           limits=(-2.5, 2.5, -0.6, 0.6))
hideydecorations!(ax1)

ax2 = Axis(fig[1, 2];
           xlabel="x", ylabel="v",
           title="Phase space",
           limits=(-2.5, 2.5, -3.5, 3.5), aspect=DataAspect())

vlines!(ax1, [-2.3]; color=:gray, linewidth=3)
mass   = Observable(Point2f(xs[1], 0.0))
spring = Observable([Point2f(-2.3, 0.0), Point2f(xs[1], 0.0)])
lines!(ax1, spring; color=:slategray, linewidth=3)
scatter!(ax1, mass; color=:crimson, markersize=32)

phase = Observable(Point2f[])
head  = Observable(Point2f(xs[1], vs[1]))
lines!(ax2, phase; color=:dodgerblue, linewidth=2)
scatter!(ax2, head; color=:crimson, markersize=16)

outpath = joinpath(@__DIR__, "..", "output", "videos", "04_spring_mass.mp4")
record(fig, outpath, eachindex(ts); framerate=30) do i
    mass[]   = Point2f(xs[i], 0.0)
    spring[] = [Point2f(-2.3, 0.0), Point2f(xs[i], 0.0)]
    push!(phase[], Point2f(xs[i], vs[i])); notify(phase)
    head[] = Point2f(xs[i], vs[i])
end

@info "wrote $outpath"
