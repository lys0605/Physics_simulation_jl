using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
include(joinpath(@__DIR__, "..", "src", "mechanics", "pendulum.jl"))

L = 1.0
θ0 = 1.2   # ~69° — small-angle approx starts to break down
sol = simulate_pendulum(θ0=θ0, L=L, tspan=(0.0, 10.0))

ts = range(sol.t[1], sol.t[end]; length=300)

fig = Figure(size=(1200, 540))
ax1 = Axis(fig[1, 1];
           xlabel="x (m)", ylabel="y (m)",
           title="Pendulum: nonlinear (blue) vs small-angle (orange)",
           limits=(-1.4, 1.4, -1.4, 0.4), aspect=DataAspect())
ax2 = Axis(fig[1, 2];
           xlabel="t (s)", ylabel="θ (rad)",
           title="Angle vs time",
           limits=(ts[1], ts[end], -1.5, 1.5))

rod_nl = Observable([Point2f(0, 0), Point2f(0, 0)])
rod_sa = Observable([Point2f(0, 0), Point2f(0, 0)])
bob_nl = Observable(Point2f(0, 0))
bob_sa = Observable(Point2f(0, 0))

lines!(ax1, rod_nl; color=:dodgerblue, linewidth=3)
lines!(ax1, rod_sa; color=:orange, linewidth=3, linestyle=:dash)
scatter!(ax1, bob_nl; color=:dodgerblue, markersize=22)
scatter!(ax1, bob_sa; color=:orange, markersize=18)

θ_nl = Observable(Point2f[])
θ_sa = Observable(Point2f[])
lines!(ax2, θ_nl; color=:dodgerblue, linewidth=2)
lines!(ax2, θ_sa; color=:orange, linewidth=2, linestyle=:dash)

outpath = joinpath(@__DIR__, "..", "output", "videos", "02_pendulum.mp4")
record(fig, outpath, eachindex(ts); framerate=30) do i
    t = ts[i]
    θ1 = sol(t)[1]
    θ2 = small_angle(t; θ0=θ0, L=L)

    p1 = Point2f(L * sin(θ1), -L * cos(θ1))
    p2 = Point2f(L * sin(θ2), -L * cos(θ2))
    rod_nl[] = [Point2f(0, 0), p1]
    rod_sa[] = [Point2f(0, 0), p2]
    bob_nl[] = p1
    bob_sa[] = p2

    push!(θ_nl[], Point2f(t, θ1)); notify(θ_nl)
    push!(θ_sa[], Point2f(t, θ2)); notify(θ_sa)
end

@info "wrote $outpath"
