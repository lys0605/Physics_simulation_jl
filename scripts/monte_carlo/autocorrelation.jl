using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/monte_carlo/mc_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Autocorrelation, integrated time, and binning errors ─────────────────────
#
# Reuses the Metropolis Ising chain at three temperatures (subcritical,
# critical, supercritical).  For each chain we show:
#
#   (a) the time series of the magnetisation
#   (b) the autocorrelation function ρ(t)
#   (c) the binning curve  σ(b)  vs block size b — plateau ≈ correct stderr
#
# At T = T_c the series is highly autocorrelated; naive √N error bars under-
# estimate by a factor √(2·τ_int).

@inline function neighbour_sum(σ, i, j, L)
    σ[mod1(i+1,L), j] + σ[mod1(i-1,L), j] + σ[i, mod1(j+1,L)] + σ[i, mod1(j-1,L)]
end

function metropolis_sweep!(σ::Matrix{Int8}, β::Float64)
    L = size(σ, 1)
    @inbounds for _ in 1:L*L
        i = rand(1:L); j = rand(1:L)
        ΔE = 2 * σ[i,j] * neighbour_sum(σ, i, j, L)
        if ΔE ≤ 0 || rand() < exp(-β*ΔE)
            σ[i,j] = -σ[i,j]
        end
    end
end

function magnetisation_series(L, T; warmup=3000, nmeas=20_000, seed=33)
    Random.seed!(seed)
    σ = rand(Int8[-1, 1], L, L)
    β = 1/T
    for _ in 1:warmup
        metropolis_sweep!(σ, β)
    end
    Ms = Float64[]
    for _ in 1:nmeas
        metropolis_sweep!(σ, β)
        push!(Ms, sum(σ) / (L*L))
    end
    Ms
end

L = 24
Ts = [2.0, 2.269, 2.8]
labels = ["T = 2.0 (ordered)", "T = T_c = 2.269", "T = 2.8 (disordered)"]
series = [magnetisation_series(L, T) for T in Ts]

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1400, 900))
Label(fig[0, :], "Autocorrelation and block-binning errors  —  Metropolis Ising,  L = $L";
      fontsize=14, font=:bold)
colors = [:dodgerblue, :crimson, :seagreen]

for (k, (T, lab, m, c)) in enumerate(zip(Ts, labels, series, colors))
    ax = Axis(fig[1, k]; title="(a.$k) m(t)  —  $lab", xlabel="sweep", ylabel="m")
    lines!(ax, 1:length(m), m; color=c, linewidth=0.6)
    ylims!(ax, -1.0, 1.0)

    ρ  = autocorrelation(m; cmax=400)
    τ  = round(tau_int(m); digits=1)
    ax2 = Axis(fig[2, k]; title="(b.$k)  ρ(t),  τ_int ≈ $τ",
               xlabel="lag t", ylabel="ρ(t)", limits=(0, 400, -0.2, 1.05))
    lines!(ax2, 0:(length(ρ)-1), ρ; color=c, linewidth=1.5)
    hlines!(ax2, [0.0]; color=:black, linewidth=0.5)

    bs, σs = block_curve(m)
    naive  = std(m) / sqrt(length(m))
    ax3 = Axis(fig[3, k]; title="(c.$k)  σ(block size)",
               xlabel="block size b", ylabel="σ_block",
               xscale=log10, yscale=log10)
    scatter!(ax3, bs, σs; color=c, markersize=8)
    lines!  (ax3, bs, σs; color=c, linewidth=1.5)
    hlines! (ax3, [naive]; color=:black, linestyle=:dash,
             label="naive σ/√N")
    axislegend(ax3; position=:lt, labelsize=9)
end

Label(fig[4, :], "Plateau of σ(b) is the correct standard error.\n" *
                 "Naive σ/√N under-estimates by √(2 τ_int) when τ_int > 1.",
      fontsize=10, color=:gray)
rowsize!(fig.layout, 4, Fixed(40))

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/autocorrelation.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
