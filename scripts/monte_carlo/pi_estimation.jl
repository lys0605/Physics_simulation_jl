using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Estimating π by hit-or-miss Monte Carlo ──────────────────────────────────
#
# Sample (x, y) uniformly in the unit square. The fraction inside the
# inscribed quarter-circle estimates π/4:
#
#       fraction inside  ≈  π/4   ⟹   π ≈ 4 · fraction
#
# Two demonstrations:
#   (a) one fixed run with N=2000 points (visualisation)
#   (b) error vs N (loglog) showing the universal 1/√N scaling

Random.seed!(2)

# ── (a) Single run, visualised ────────────────────────────────────────────────
N_vis = 2000
xs    = rand(N_vis);  ys = rand(N_vis)
inside = xs.^2 .+ ys.^2 .≤ 1.0
π_est  = 4 * count(inside) / N_vis

# ── (b) Error scaling: many runs at increasing N ─────────────────────────────
Ns       = [10^k for k in 1:6]                  # 10, 100, … 10⁶
trials   = 200
errs_med = Float64[];  errs_qup = Float64[];  errs_qlo = Float64[]
for N in Ns
    errs = Float64[]
    for _ in 1:trials
        x = rand(N);  y = rand(N)
        push!(errs, abs(4 * count(x.^2 .+ y.^2 .≤ 1.0) / N - π))
    end
    push!(errs_med, median(errs))
    push!(errs_qlo, quantile(errs, 0.16))
    push!(errs_qup, quantile(errs, 0.84))
end

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1100, 500))
Label(fig[0, :], "Estimating π by Monte Carlo  —  hit-or-miss sampling";
      fontsize=14, font=:bold)

ax1 = Axis(fig[1,1]; aspect=DataAspect(),
           title="(a) $N_vis samples,  π ≈ $(round(π_est, digits=4))",
           xlabel="x", ylabel="y", limits=(0,1,0,1))
scatter!(ax1, xs[inside],  ys[inside];  color=:dodgerblue, markersize=3)
scatter!(ax1, xs[.!inside], ys[.!inside]; color=:tomato,    markersize=3)
θ = range(0, π/2; length=200)
lines!(ax1, cos.(θ), sin.(θ); color=:black, linewidth=2)

ax2 = Axis(fig[1,2]; xscale=log10, yscale=log10,
           title="(b)  |π̂ − π| vs N    (slope = −1/2)",
           xlabel="number of samples N", ylabel="error")
band!(ax2, Ns, errs_qlo, errs_qup; color=(:dodgerblue, 0.25))
scatter!(ax2, Ns, errs_med; color=:dodgerblue, markersize=10)
lines!(ax2, Ns, errs_med;  color=:dodgerblue, linewidth=2)
# Reference 1/√N
ref = errs_med[1] .* sqrt.(Ns[1] ./ Ns)
lines!(ax2, Ns, ref; color=:black, linestyle=:dash, label="1/√N")
axislegend(ax2; position=:rt)

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/pi_estimation.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
