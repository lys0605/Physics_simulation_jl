using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/monte_carlo/mc_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── 2-D Ising model — Metropolis single-spin-flip ────────────────────────────
#
#   H = -J Σ_<ij> σ_i σ_j     with J = 1, periodic boundaries
#   Onsager critical temperature  T_c = 2J / ln(1 + √2) ≈ 2.269  (J = k_B = 1)
#
# Outputs:
#   (a) energy per site E/N  vs T
#   (b) absolute magnetisation ⟨|m|⟩ vs T
#   (c) susceptibility χ = β·N·(⟨m²⟩−⟨|m|⟩²) vs T  →  peak locates T_c
#   (d) specific heat  C_V = β² · N · Var(E_per_site) vs T  →  peak at T_c
#   (e) three snapshots at T = 1.8, T_c, 3.0

# ── Energy / move helpers ─────────────────────────────────────────────────────
@inline function neighbour_sum(σ, i, j, L)
    σ[mod1(i+1, L), j] + σ[mod1(i-1, L), j] +
    σ[i, mod1(j+1, L)] + σ[i, mod1(j-1, L)]
end

function total_energy(σ::Matrix{Int8})
    L = size(σ, 1);  E = 0
    @inbounds for i in 1:L, j in 1:L
        # Only count right and down neighbours (avoid double counting)
        E -= σ[i, j] * (σ[mod1(i+1, L), j] + σ[i, mod1(j+1, L)])
    end
    E
end

function metropolis_sweep!(σ::Matrix{Int8}, β::Float64, exp_dE::Dict{Int,Float64})
    L = size(σ, 1);  ΔE_total = 0
    @inbounds for _ in 1:L*L
        i = rand(1:L);  j = rand(1:L)
        s = σ[i, j];   nn = neighbour_sum(σ, i, j, L)
        ΔE = 2 * s * nn
        if ΔE ≤ 0 || rand() < exp_dE[ΔE]
            σ[i, j] = -s
            ΔE_total += ΔE
        end
    end
    ΔE_total
end

function run_at_T(L, T; warmup=2000, nmeas=4000, seed=2026)
    Random.seed!(seed)
    σ = rand(Int8[-1, 1], L, L)
    β = 1 / T
    exp_dE = Dict(ΔE => exp(-β*ΔE) for ΔE in (-8,-4,0,4,8))
    E = Float64(total_energy(σ))

    for _ in 1:warmup
        E += metropolis_sweep!(σ, β, exp_dE)
    end
    Es, Ms = Float64[], Float64[]
    for _ in 1:nmeas
        E += metropolis_sweep!(σ, β, exp_dE)
        push!(Es, E)
        push!(Ms, sum(σ))
    end
    σ, Es, Ms
end

# ── Sweep T ───────────────────────────────────────────────────────────────────
L   = 24
Ts  = vcat(range(1.5, 2.1; length=8), range(2.15, 2.4; length=10),
           range(2.45, 3.2; length=8))
E_means = Float64[];  M_means = Float64[]
chis    = Float64[];  cvs     = Float64[]

for T in Ts
    _, Es, Ms = run_at_T(L, T; warmup=1500, nmeas=4000)
    N = L*L
    push!(E_means, mean(Es) / N)
    push!(M_means, mean(abs.(Ms)) / N)
    β = 1/T
    push!(chis, β/N * (mean(Ms.^2) - mean(abs.(Ms))^2))
    push!(cvs,  β^2/N * (mean(Es.^2) - mean(Es)^2))
end

# ── Snapshots ────────────────────────────────────────────────────────────────
T_c_onsager = 2 / log(1 + sqrt(2))
snap_Ts = [1.8, T_c_onsager, 3.0]
snaps   = [first(run_at_T(L, T; warmup=2000, nmeas=1)) for T in snap_Ts]

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1400, 900))
Label(fig[0, :], "2-D Ising Metropolis  (L = $L)  —  Onsager  T_c ≈ $(round(T_c_onsager, digits=3))";
      fontsize=14, font=:bold)

ax1 = Axis(fig[1,1]; title="(a) E / N", xlabel="T", ylabel="energy per site")
lines!(ax1, Ts, E_means; color=:dodgerblue, linewidth=2)
scatter!(ax1, Ts, E_means; color=:dodgerblue, markersize=6)
vlines!(ax1, [T_c_onsager]; color=(:gray, 0.6), linestyle=:dash)

ax2 = Axis(fig[1,2]; title="(b) ⟨|m|⟩", xlabel="T", ylabel="magnetisation")
lines!(ax2, Ts, M_means; color=:crimson, linewidth=2)
scatter!(ax2, Ts, M_means; color=:crimson, markersize=6)
vlines!(ax2, [T_c_onsager]; color=(:gray, 0.6), linestyle=:dash)

ax3 = Axis(fig[1,3]; title="(c) χ  (susceptibility)", xlabel="T", ylabel="χ")
lines!(ax3, Ts, chis; color=:darkorange, linewidth=2)
scatter!(ax3, Ts, chis; color=:darkorange, markersize=6)
vlines!(ax3, [T_c_onsager]; color=(:gray, 0.6), linestyle=:dash)

ax4 = Axis(fig[2,1]; title="(d) Cᵥ  (specific heat)", xlabel="T", ylabel="C_V")
lines!(ax4, Ts, cvs; color=:seagreen, linewidth=2)
scatter!(ax4, Ts, cvs; color=:seagreen, markersize=6)
vlines!(ax4, [T_c_onsager]; color=(:gray, 0.6), linestyle=:dash)

for (k, (T, σ)) in enumerate(zip(snap_Ts, snaps))
    ax = Axis(fig[2, k+1]; aspect=DataAspect(),
              title="(e.$k) T = $(round(T, digits=3))",
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false)
    heatmap!(ax, 1:L, 1:L, σ; colormap=:RdBu_r, colorrange=(-1, 1))
end

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/ising_metropolis.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
