using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/monte_carlo/mc_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Variational Monte Carlo for the 1-D harmonic oscillator ──────────────────
#
#   H = −½ ∂_x² + ½ x²          (ℏ = ω = m = 1)
#   trial Ψ_α(x) = exp(−α x²)
#
# Local energy E_L(x) = (HΨ)/Ψ:
#   Ψ' = −2α x Ψ
#   Ψ'' = (4α²x² − 2α) Ψ
#   E_L = α + (½ − 2α²) x²
#
# At α = ½ the variance vanishes and E_L = ½ everywhere — the exact ground
# state.  Plot ⟨E_L⟩(α), σ(E_L)(α), histogram of |Ψ|² for optimal α, and the
# stochastic-gradient descent of α toward the optimum.

E_local(α, x)   = α + (0.5 - 2α^2) * x^2
log_ψ²(α, x)    = -2α * x^2

# Sample |Ψ|² by Metropolis on x
function vmc_chain(α; nwalkers=500, nsteps=4000, σ_prop=1.0, warmup=500)
    xs = randn(nwalkers)
    El_acc = Float64[]
    grad_acc = Float64[]                # gradient of ⟨E⟩ wrt α via SR-free form
    for step in 1:(warmup + nsteps)
        xs_new = xs .+ σ_prop .* (rand(nwalkers) .* 2 .- 1)
        Δlog   = log_ψ²(α, xs_new) .- log_ψ²(α, xs)
        accept = rand(nwalkers) .< exp.(Δlog)
        xs = ifelse.(accept, xs_new, xs)
        if step > warmup
            EL = E_local.(α, xs)
            append!(El_acc, EL)
            # d log ψ²/dα = −2 x²;  gradient of ⟨E⟩ ≈ 2(⟨E·O⟩ − ⟨E⟩⟨O⟩),  O = ∂_α log ψ
            O = -xs.^2
            append!(grad_acc, EL .* O .- mean(EL) .* O)
        end
    end
    El_acc, grad_acc, xs
end

# ── Sweep α and gather ⟨E⟩, σ_E ──────────────────────────────────────────────
αs = range(0.20, 0.95; length=20)
E_means = Float64[];  E_vars = Float64[]
for α in αs
    Random.seed!(42)
    El, _, _ = vmc_chain(α; nwalkers=400, nsteps=2000, warmup=300)
    push!(E_means, mean(El))
    push!(E_vars,  var(El))
end

# ── Stochastic optimisation of α ─────────────────────────────────────────────
Random.seed!(2026)
α_traj = [0.30]
α      = α_traj[1]
η      = 0.04
for it in 1:50
    El, grads, _ = vmc_chain(α; nwalkers=300, nsteps=600, warmup=200)
    g = 2 * mean(grads)
    α = α - η * g
    push!(α_traj, α)
end

# Optimal-α walker distribution
Random.seed!(7)
_, _, xs_opt = vmc_chain(0.5; nwalkers=4000, nsteps=400, warmup=300)

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1300, 850))
Label(fig[0, :], "Variational MC for the harmonic oscillator  (exact:  E₀ = ½, α* = ½)";
      fontsize=14, font=:bold)

ax1 = Axis(fig[1,1]; title="(a) ⟨E_L⟩(α)", xlabel="α", ylabel="energy")
lines!(ax1, αs, E_means; color=:dodgerblue, linewidth=2)
scatter!(ax1, αs, E_means; color=:dodgerblue, markersize=7)
hlines!(ax1, [0.5]; color=:black, linestyle=:dash, label="exact ½")
vlines!(ax1, [0.5]; color=(:gray, 0.5), linestyle=:dash)
axislegend(ax1; position=:lt)

ax2 = Axis(fig[1,2]; title="(b) Var(E_L)(α)  →  zero-variance principle",
           xlabel="α", ylabel="Var(E_L)")
lines!(ax2, αs, E_vars; color=:crimson, linewidth=2)
scatter!(ax2, αs, E_vars; color=:crimson, markersize=7)
vlines!(ax2, [0.5]; color=(:gray, 0.5), linestyle=:dash)

ax3 = Axis(fig[1,3]; title="(c) Histogram of |Ψ_½(x)|²",
           xlabel="x", ylabel="density")
hist!(ax3, xs_opt; bins=80, color=(:dodgerblue, 0.5), strokecolor=:dodgerblue,
      strokewidth=0.5, normalization=:pdf)
xg = range(-4, 4; length=200)
lines!(ax3, xg, exp.(-xg.^2) ./ sqrt(π/2); color=:crimson, linewidth=2,
       label="exact |Ψ_½|² ∝ e^{-x²}")
axislegend(ax3; position=:rt)

ax4 = Axis(fig[2, 1:3]; title="(d) Stochastic-gradient descent of α",
           xlabel="iteration", ylabel="α")
lines!(ax4, 0:(length(α_traj)-1), α_traj; color=:purple, linewidth=2)
scatter!(ax4, 0:(length(α_traj)-1), α_traj; color=:purple, markersize=6)
hlines!(ax4, [0.5]; color=:black, linestyle=:dash, label="α* = ½")
axislegend(ax4; position=:rt)

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/vmc_harmonic.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
