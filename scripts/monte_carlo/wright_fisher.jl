using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# Simple Binomial(N, p) sampler — sum of N Bernoullis (fine for N ≲ 10⁴)
rand_binomial(N, p) = count(rand() < p for _ in 1:N)

# ── Wright-Fisher model — genetic drift in a finite population ───────────────
#
# N diploid individuals at a locus with two alleles A (frequency p) and a
# (frequency 1−p).  Each generation, draw N gametes independently with the
# parent allele frequencies; the new allele count is Binomial(N, p).
#
#   p_{t+1} | p_t  ~  (1/N) Binomial(N, p_t)
#
# With selection coefficient s (advantage of A), use the biased frequency
#   p_t' = p_t (1 + s) / (1 + s p_t),
# and sample Binomial(N, p_t').
#
# Outputs:
#   (a) 30 neutral trajectories  — drift, occasional fixation/extinction
#   (b) 30 trajectories with s = 0.02 — selection sweep
#   (c) fixation probability vs initial frequency p_0
#       neutral:  P_fix = p_0
#       weak sel: P_fix = (1 − e^{−2 s p_0}) / (1 − e^{−2 s})    (Kimura)

Random.seed!(2026)
const N = 200             # population size

function wf_trajectory(N, p0, T; s=0.0)
    p = zeros(Float64, T + 1)
    p[1] = p0
    for t in 1:T
        p_eff = p[t] * (1 + s) / (1 + s * p[t])
        p_eff = clamp(p_eff, 0.0, 1.0)
        p[t+1] = rand_binomial(N, p_eff) / N
    end
    p
end

T = 200
trajs_neutral   = [wf_trajectory(N, 0.5, T; s=0.0)  for _ in 1:30]
trajs_selection = [wf_trajectory(N, 0.05, T; s=0.04) for _ in 1:30]

# ── Fixation probability sweep ───────────────────────────────────────────────
function fixation_prob(N, p0; s=0.0, T_max=10_000, n_trial=400)
    fixed = 0
    for _ in 1:n_trial
        p = p0
        for _ in 1:T_max
            p_eff = p * (1 + s) / (1 + s*p)
            p = rand_binomial(N, clamp(p_eff, 0.0, 1.0)) / N
            (p == 0.0 || p == 1.0) && break
        end
        fixed += p ≈ 1.0
    end
    fixed / n_trial
end

p0_grid = range(0.05, 0.95; length=8)
P_neutral_sim  = [fixation_prob(N, p; s=0.0) for p in p0_grid]
P_selection_sim = [fixation_prob(N, p; s=0.02) for p in p0_grid]

# Kimura's formula for arbitrary selection (s > 0)
kimura(p, N, s) = s == 0 ? p :
                   (1 - exp(-2 * 2 * N * s * p / (2 * N))) /
                   (1 - exp(-2 * 2 * N * s / (2 * N)))   # = (1−e^{−2sNp/N}) / (1−e^{−2s})

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1300, 700))
Label(fig[0, :], "Wright-Fisher Genetic Drift  (N = $N)";
      fontsize=14, font=:bold)

ax1 = Axis(fig[1,1]; title="(a) neutral drift, p₀ = 0.5",
           xlabel="generation", ylabel="frequency of A")
for tr in trajs_neutral
    lines!(ax1, 0:T, tr; color=(:dodgerblue, 0.45), linewidth=0.8)
end
ylims!(ax1, 0, 1)

ax2 = Axis(fig[1,2]; title="(b) selection s = 0.04, p₀ = 0.05",
           xlabel="generation", ylabel="frequency of A")
for tr in trajs_selection
    lines!(ax2, 0:T, tr; color=(:crimson, 0.55), linewidth=0.8)
end
ylims!(ax2, 0, 1)

ax3 = Axis(fig[1,3]; title="(c) fixation probability",
           xlabel="initial frequency p₀", ylabel="P_fix",
           limits=(0, 1, 0, 1))
# Neutral theory: P_fix = p_0
lines!(ax3, 0:0.01:1, 0:0.01:1; color=:black, linestyle=:dash, label="P=p₀ (neutral)")
scatter!(ax3, p0_grid, P_neutral_sim;  color=:dodgerblue, markersize=10,
         label="neutral, simulated")
# With selection s = 0.02
ps_grid = range(0.01, 0.99; length=50)
lines!(ax3, ps_grid, kimura.(ps_grid, N, 0.02);
       color=:crimson, label="Kimura, s = 0.02")
scatter!(ax3, p0_grid, P_selection_sim; color=:crimson, markersize=10,
         marker=:diamond)
axislegend(ax3; position=:lt, labelsize=9)

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/wright_fisher.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
