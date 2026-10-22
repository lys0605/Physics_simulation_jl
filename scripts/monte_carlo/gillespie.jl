using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Gillespie's Stochastic Simulation Algorithm (SSA, 1976) ──────────────────
#
# Exact stochastic simulation of a chemical / biological reaction network.
# Demonstrated on the Lotka-Volterra predator-prey reactions:
#
#   R ──────►  R + R          (rate  α · R)            prey birth
#   R + F ──►  F + F          (rate  β · R · F)        predation
#   F ──────►  ∅              (rate  γ · F)            predator death
#
# At each step:
#   1. compute propensities  a_i  for every reaction
#   2. sample  Δt  ~  Exp(Σ a_i)
#   3. sample reaction k with probability  a_k / Σ a_i
#   4. update populations, advance time
#
# Outputs:
#   (a) one trajectory R(t), F(t)
#   (b) phase portrait F(R)
#   (c) ensemble of 30 trajectories — variation due to demographic noise

const α = 1.0          # prey birth rate
const β = 0.005        # predation rate
const γ = 0.6          # predator death rate

function gillespie_lv(R0, F0, T_max; rng=Random.default_rng())
    R = R0;  F = F0;  t = 0.0
    ts = [t];  Rs = [R];  Fs = [F]
    while t < T_max && R > 0 && F > 0
        a1 = α * R                # birth
        a2 = β * R * F            # predation
        a3 = γ * F                # death
        a0 = a1 + a2 + a3
        a0 == 0 && break
        τ = -log(rand(rng)) / a0
        u = rand(rng) * a0
        if u < a1
            R += 1
        elseif u < a1 + a2
            R -= 1;  F += 1
        else
            F -= 1
        end
        t += τ
        push!(ts, t);  push!(Rs, R);  push!(Fs, F)
    end
    ts, Rs, Fs
end

Random.seed!(2026)
ts, Rs, Fs = gillespie_lv(120, 40, 30.0)

# Ensemble
runs = [gillespie_lv(120, 40, 30.0) for _ in 1:30]

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1300, 720))
Label(fig[0, :], "Stochastic Lotka-Volterra via Gillespie SSA  (α=$α, β=$β, γ=$γ)";
      fontsize=14, font=:bold)

ax1 = Axis(fig[1, 1:2]; title="(a) populations over time", xlabel="t",
           ylabel="count", limits=(0, 30, 0, nothing))
lines!(ax1, ts, Rs; color=:dodgerblue, linewidth=1.6, label="rabbits R(t)")
lines!(ax1, ts, Fs; color=:crimson,    linewidth=1.6, label="foxes F(t)")
axislegend(ax1; position=:rt)

ax2 = Axis(fig[1, 3]; title="(b) phase portrait R–F", xlabel="R", ylabel="F",
           aspect=DataAspect())
lines!(ax2, Rs, Fs; color=:purple, linewidth=1.2)
scatter!(ax2, [Rs[1]], [Fs[1]]; color=:black, markersize=10, marker=:star5)

ax3 = Axis(fig[2, 1:3]; title="(c) ensemble of 30 trajectories",
           xlabel="t", ylabel="R(t)")
for (ts_r, Rs_r, _) in runs
    lines!(ax3, ts_r, Rs_r; color=(:dodgerblue, 0.35), linewidth=0.8)
end
xlims!(ax3, 0, 30)

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/gillespie.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
