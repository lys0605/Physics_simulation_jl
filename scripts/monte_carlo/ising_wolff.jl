using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/monte_carlo/mc_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Wolff cluster algorithm for the 2-D Ising model ──────────────────────────
#
# Builds a Fortuin-Kasteleyn cluster of aligned spins via random bond
# activation with probability  p = 1 − exp(−2β J)  and flips it as a unit.
# Detailed balance is preserved at the cluster level.
#
# Near T_c, autocorrelation collapses from τ ∼ L² (Metropolis) to τ ∼ L^0.2.
# This script compares τ_int for ⟨m²⟩ between Metropolis and Wolff at several
# temperatures, with the dramatic gap visible right at T_c.

@inline neighbours(i, j, L) = (
    (mod1(i+1, L), j), (mod1(i-1, L), j),
    (i, mod1(j+1, L)), (i, mod1(j-1, L)),
)

# ── Wolff sweep: flip one cluster per call ───────────────────────────────────
function wolff_step!(σ::Matrix{Int8}, β::Float64)
    L = size(σ, 1)
    p = 1 - exp(-2β)
    i0, j0 = rand(1:L), rand(1:L)
    s0 = σ[i0, j0]
    stack = [(i0, j0)]
    visited = falses(L, L);  visited[i0, j0] = true
    σ[i0, j0] = -s0
    cluster_size = 1
    while !isempty(stack)
        i, j = pop!(stack)
        for (ni, nj) in neighbours(i, j, L)
            if !visited[ni, nj] && σ[ni, nj] == s0 && rand() < p
                visited[ni, nj] = true
                σ[ni, nj] = -s0
                cluster_size += 1
                push!(stack, (ni, nj))
            end
        end
    end
    cluster_size
end

# ── Metropolis sweep (single-spin flips per site) ────────────────────────────
@inline function neighbour_sum(σ, i, j, L)
    σ[mod1(i+1,L), j] + σ[mod1(i-1,L), j] + σ[i, mod1(j+1,L)] + σ[i, mod1(j-1,L)]
end

function metropolis_sweep!(σ::Matrix{Int8}, β::Float64)
    L = size(σ, 1)
    @inbounds for _ in 1:L*L
        i = rand(1:L); j = rand(1:L)
        ΔE = 2 * σ[i,j] * neighbour_sum(σ, i, j, L)
        if ΔE ≤ 0 || rand() < exp(-β*ΔE)
            σ[i, j] = -σ[i, j]
        end
    end
end

# ── Run comparison ────────────────────────────────────────────────────────────
function run_chain(updater!, L, T; warmup=2000, nmeas=4000, seed=10)
    Random.seed!(seed)
    σ = rand(Int8[-1, 1], L, L)
    β = 1/T
    for _ in 1:warmup
        updater!(σ, β)
    end
    Ms = Float64[]
    for _ in 1:nmeas
        updater!(σ, β)
        push!(Ms, sum(σ)^2 / (L*L)^2)
    end
    Ms
end

L  = 24
Ts = [2.0, 2.15, 2.269, 2.4, 2.8]
τ_metro = Float64[];  τ_wolff = Float64[]
for T in Ts
    Ms_m = run_chain((σ, β) -> metropolis_sweep!(σ, β), L, T)
    Ms_w = run_chain((σ, β) -> wolff_step!(σ, β),       L, T)
    push!(τ_metro, tau_int(Ms_m))
    push!(τ_wolff, tau_int(Ms_w))
end

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1100, 500))
Label(fig[0, :], "Metropolis vs Wolff — integrated autocorrelation τ_int(⟨m²⟩)";
      fontsize=14, font=:bold)

ax1 = Axis(fig[1,1]; xlabel="T", ylabel="τ_int (sweeps / cluster moves)",
           title="Autocorrelation at L = $L", yscale=log10)
scatter!(ax1, Ts, τ_metro; color=:crimson,    markersize=12, label="Metropolis")
lines!  (ax1, Ts, τ_metro; color=:crimson,    linewidth=2)
scatter!(ax1, Ts, τ_wolff; color=:dodgerblue, markersize=12, label="Wolff")
lines!  (ax1, Ts, τ_wolff; color=:dodgerblue, linewidth=2)
vlines! (ax1, [2.269];     color=(:gray, 0.5), linestyle=:dash)
text!(ax1, 2.275, maximum(τ_metro)*0.9; text="T_c", color=:gray)
axislegend(ax1; position=:rt)

ax2 = Axis(fig[1,2]; xlabel="T", ylabel="speed-up τ_metro / τ_wolff",
           title="Wolff acceleration factor", yscale=log10)
ratios = τ_metro ./ τ_wolff
scatter!(ax2, Ts, ratios; color=:purple, markersize=12)
lines!  (ax2, Ts, ratios; color=:purple, linewidth=2)
vlines! (ax2, [2.269]; color=(:gray, 0.5), linestyle=:dash)

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/ising_wolff.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
