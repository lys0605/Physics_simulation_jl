using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Drossel-Schwabl forest-fire cellular automaton (1992) ────────────────────
#
# States per cell:  0 = empty,  1 = tree,  2 = burning
# Rules (synchronous update each step):
#   • burning  → empty
#   • tree     → burning  if any neighbour is burning
#   • tree     → burning  with prob f  (lightning)
#   • empty    → tree     with prob g
#
# Self-organised critical: cluster sizes follow a power-law distribution
# at the (f/g → 0) limit.  This script animates the evolution and plots
# the cluster-size distribution averaged over the steady state.

Random.seed!(11)

const L  = 200
const G  = 0.02     # tree growth probability
const F  = 5e-5     # lightning probability
const T_BURN_IN = 400
const T_STEADY  = 80
const FRAMERATE = 30

state = zeros(Int8, L, L)

# ── Update rule ──────────────────────────────────────────────────────────────
function ff_step!(state)
    L = size(state, 1)
    new_state = copy(state)
    @inbounds for i in 1:L, j in 1:L
        if state[i, j] == 2
            new_state[i, j] = 0
        elseif state[i, j] == 1
            burning = false
            for (di, dj) in ((1,0),(-1,0),(0,1),(0,-1))
                ni, nj = mod1(i+di, L), mod1(j+dj, L)
                if state[ni, nj] == 2
                    burning = true; break
                end
            end
            new_state[i, j] = burning ? 2 : (rand() < F ? 2 : 1)
        else  # empty
            new_state[i, j] = rand() < G ? 1 : 0
        end
    end
    state .= new_state
end

# ── Cluster sizes of tree component (flood fill) ─────────────────────────────
function tree_cluster_sizes(state)
    L = size(state, 1)
    labels = zeros(Int, L, L)
    sizes  = Int[]
    next = 1
    for i in 1:L, j in 1:L
        if state[i, j] == 1 && labels[i, j] == 0
            stack = [(i, j)];  labels[i, j] = next;  sz = 0
            while !isempty(stack)
                ci, cj = pop!(stack);  sz += 1
                for (di, dj) in ((1,0),(-1,0),(0,1),(0,-1))
                    ni, nj = mod1(ci+di, L), mod1(cj+dj, L)
                    if state[ni, nj] == 1 && labels[ni, nj] == 0
                        labels[ni, nj] = next; push!(stack, (ni, nj))
                    end
                end
            end
            push!(sizes, sz); next += 1
        end
    end
    sizes
end

# ── Burn-in (no animation) ───────────────────────────────────────────────────
@info "burning in for $T_BURN_IN steps…"
for _ in 1:T_BURN_IN
    ff_step!(state)
end

# Collect cluster sizes from steady state
@info "collecting cluster statistics…"
all_sizes = Int[]
for _ in 1:200
    ff_step!(state)
    append!(all_sizes, tree_cluster_sizes(state))
end

# ── Animate ──────────────────────────────────────────────────────────────────
set_theme!(theme_dark())
fig = Figure(size=(1000, 540))
Label(fig[0, :], "Drossel-Schwabl Forest Fire   (L=$L,  g=$G,  f=$F)";
      fontsize=14, font=:bold)

ax = Axis(fig[1,1]; aspect=DataAspect(),
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false,
          title="state map  (green=tree, red=burning)")
state_obs = Observable(copy(state))
heatmap!(ax, 1:L, 1:L, state_obs;
         colormap=cgrad([:black, :green, :red]), colorrange=(0, 2))

ax_h = Axis(fig[1,2]; xscale=log10, yscale=log10,
            title="cluster-size distribution (steady state)",
            xlabel="size s", ylabel="N(s)")

# Build histogram
bins  = 10 .^ range(log10(1), log10(maximum(all_sizes)+1); length=30)
counts = zeros(Int, length(bins)-1)
for s in all_sizes
    k = searchsortedlast(bins, s)
    if 1 ≤ k ≤ length(counts)
        counts[k] += 1
    end
end
centres = sqrt.(bins[1:end-1] .* bins[2:end])
mask    = counts .> 0
scatter!(ax_h, centres[mask], counts[mask]; color=:cyan, markersize=8)

# Power-law reference τ ≈ 2 (theoretical exponent ~2 in 2-D limit)
ref_x = centres[mask][1:min(end, 10)]
ref_y = counts[mask][1] .* (ref_x ./ ref_x[1]).^(-2)
lines!(ax_h, ref_x, ref_y; color=:white, linestyle=:dash, label="N(s) ∝ s⁻²")
axislegend(ax_h; position=:rt, labelsize=9)

mkpath(joinpath(@__DIR__, "../../output/videos"))
out = joinpath(@__DIR__, "../../output/videos/forest_fire.mp4")

record(fig, out, 1:T_STEADY; framerate=FRAMERATE) do _
    ff_step!(state)
    state_obs[] = copy(state)
end

@info "wrote $out"
