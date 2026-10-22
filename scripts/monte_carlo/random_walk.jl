using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Random walks — building block of MC, diffusion, and beyond ───────────────
#
# Two demos:
#   (a) sample 2-D trajectories from a few independent walkers
#   (b) mean-square displacement ⟨x²(t)⟩ vs t for 1-D walks — should be linear
#       in t with slope 1 (per step), confirming diffusive scaling.

Random.seed!(11)

# ── (a) 2-D trajectories ──────────────────────────────────────────────────────
N_steps   = 4000
N_walkers = 6
trajs2d   = [zeros(Float64, N_steps + 1, 2) for _ in 1:N_walkers]
dirs      = ([1,0], [-1,0], [0,1], [0,-1])
for w in 1:N_walkers
    for t in 1:N_steps
        d = rand(dirs)
        trajs2d[w][t+1, :] = trajs2d[w][t, :] .+ d
    end
end

# ── (b) MSD in 1-D ────────────────────────────────────────────────────────────
N_msd     = 2000
N_avg     = 5000
positions = zeros(Int, N_avg, N_msd + 1)
for w in 1:N_avg, t in 1:N_msd
    positions[w, t+1] = positions[w, t] + rand((-1, 1))
end
ts  = 0:N_msd
msd = [mean(positions[:, t+1].^2) for t in ts]

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1200, 520))
Label(fig[0, :], "Random Walks  —  trajectories and ⟨x²⟩ ∼ t";
      fontsize=14, font=:bold)

ax1 = Axis(fig[1,1]; aspect=DataAspect(),
           title="(a)  Six 2-D walkers, $N_steps steps each",
           xlabel="x", ylabel="y")
cols = cgrad(:viridis, N_walkers; categorical=true)
for (i, tr) in enumerate(trajs2d)
    lines!(ax1, tr[:,1], tr[:,2]; color=cols[i], linewidth=1.0)
end
scatter!(ax1, [0.0], [0.0]; color=:black, markersize=12, marker=:star5)

ax2 = Axis(fig[1,2];
           title="(b)  ⟨x²(t)⟩ over $N_avg one-dimensional walks",
           xlabel="t (steps)", ylabel="⟨x²⟩")
lines!(ax2, ts, msd; color=:dodgerblue, linewidth=2, label="simulation")
lines!(ax2, ts, Float64.(ts); color=:black, linestyle=:dash, linewidth=1.5,
       label="t  (theory)")
axislegend(ax2; position=:lt)

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/random_walk.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
