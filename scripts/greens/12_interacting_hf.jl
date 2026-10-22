using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, LinearAlgebra
include(joinpath(@__DIR__, "../../src/greens/qwz_model.jl"))
include(joinpath(@__DIR__, "../../src/greens/bulk_greens.jl"))
include(joinpath(@__DIR__, "../../src/greens/interacting_greens.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Hartree-Fock renormalization of the QWZ Chern insulator ───────────────────
#
# Interaction: H_int = U Σ_R  n_{R,1} n_{R,2}  (inter-orbital Hubbard)
#
# Because the QWZ model has inversion symmetry (k_x → −k_x), the Fock
# off-diagonal ρ₁₂ vanishes identically.  The only surviving HF term is
# the diagonal Hartree shift:
#
#   Σ^HF = diag(U ρ₂₂, U ρ₁₁)
#
# This is equivalent to renormalising the topological mass:
#
#   m_eff(U) = m + U(ρ₂₂ − ρ₁₁)/2
#
# Self-consistently, ρ₁₁ > ρ₂₂ (orbital 1 preferentially occupied for m = −1),
# and the imbalance grows with U.  Once m_eff < −2 the system crosses the
# topological phase boundary and C flips from −1 to 0.  The critical value
# is U_c ≈ 2.46 for m = −1.
#
# Panels:
#   (a–c) Spectral function A_HF(k,ω) at U = 0, 2.2, 3.0
#   (d)   Chern number C vs U (0 → 4)
#   (e)   Orbital occupancies ρ₁₁, ρ₂₂ vs U
#   (f)   Effective mass m_eff vs U with phase boundary at m = −2

const m_top   = QWZModel(-1.0)
const η       = 0.05
const Nk_hf   = 80      # BZ mesh for SCF density matrix
const Nk_chern = 60     # BZ mesh for Chern number
const Nk_dos   = 80     # BZ mesh for DOS / spectral function
const U_show  = [0.0, 2.2, 3.0]

# ── SCF for display U values ──────────────────────────────────────────────────
@info "HF self-consistent loops (display panels)…"
hf_show = [(U, hf_selfconsistent(m_top, U; Nk=Nk_hf)) for U in U_show]

# ── fine U sweep for Chern number and density matrix ─────────────────────────
U_sweep = range(0.0, 4.0; length=41)
@info "HF sweep over $(length(U_sweep)) U values…"
hf_sweep = [hf_selfconsistent(m_top, U; Nk=Nk_hf) for U in U_sweep]
Cs       = [chern_number_hf(m_top, Σ; Nk=Nk_chern) for (Σ, _) in hf_sweep]
ρ11_vec  = [real(ρ[1,1]) for (_, ρ) in hf_sweep]
ρ22_vec  = [real(ρ[2,2]) for (_, ρ) in hf_sweep]
meff_vec = [-1.0 + U*(ρ[2,2] - ρ[1,1])/2 |> real for (U, (_, ρ)) in zip(U_sweep, hf_sweep)]

# ── spectral function maps along Γ→X→M→Γ ─────────────────────────────────────
@info "Spectral function maps…"
kxs, kys, tick_pos, tick_lbl = kpath_points(50)
ωs  = range(-4.5, 4.5; length=220)
k_idx = 1:length(kxs)
A_maps = [interacting_spectral_map(m_top, Σ, kxs, kys, ωs; η)
          for (_, (Σ, _)) in hf_show]

# HF band energies along path (for overlay)
function hf_bands(model, Σ_HF, kxs, kys)
    [sort(real.(eigvals(bloch_hamiltonian(model, kxs[i], kys[i]) + Σ_HF)))
     for i in eachindex(kxs)]
end
bands_show = [hf_bands(m_top, Σ, kxs, kys) for (_, (Σ, _)) in hf_show]

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1400, 900))

# ── (a–c) spectral heatmaps ───────────────────────────────────────────────────
panel_labels = ["(a)", "(b)", "(c)"]
for (col, (U, (Σ, ρ), A, bands)) in
        enumerate(zip(U_show,
                      [hf_show[i][2] for i in eachindex(hf_show)],
                      A_maps, bands_show))
    C = chern_number_hf(m_top, Σ; Nk=Nk_chern)
    ax = Axis(fig[1, col];
              title     = "$(panel_labels[col])  U = $U  (C = $C)",
              xlabel    = col == 2 ? "k" : "",
              ylabel    = col == 1 ? "ω" : "",
              limits    = (1, length(kxs), -4.5, 4.5),
              xticks    = (tick_pos, tick_lbl),
              titlesize = 12,
              yticklabelsvisible = col == 1)

    heatmap!(ax, k_idx, ωs, A'; colormap = :hot, colorrange = (0, 2.5))
    E_lo = [b[1] for b in bands];  E_hi = [b[2] for b in bands]
    lines!(ax, k_idx, E_lo; color = :dodgerblue, linewidth = 1.5, linestyle = :dash)
    lines!(ax, k_idx, E_hi; color = :dodgerblue, linewidth = 1.5, linestyle = :dash)
    for tp in tick_pos
        vlines!(ax, tp; color = (:white, 0.35), linewidth = 1)
    end
end

# ── (d) Chern number vs U ─────────────────────────────────────────────────────
ax_c = Axis(fig[2, 1];
            title   = "(d) Chern number vs U",
            xlabel  = "U",
            ylabel  = "C",
            yticks  = ([-1, 0], ["-1", "0"]),
            limits  = (0, 4, -1.4, 0.4),
            titlesize = 12)
vlines!(ax_c, [2.46]; color = (:gray, 0.5), linewidth = 1.5, linestyle = :dash)
scatter!(ax_c, collect(U_sweep), Float64.(Cs); color = :crimson, markersize = 8)
lines!(ax_c, collect(U_sweep), Float64.(Cs); color = :crimson, linewidth = 2)
text!(ax_c, 2.46, -0.2; text = "U_c ≈ 2.46", color = :gray,
      align = (:left, :bottom), fontsize = 10)

# ── (e) orbital occupancies ρ₁₁, ρ₂₂ vs U ───────────────────────────────────
ax_rho = Axis(fig[2, 2];
              title   = "(e) Orbital occupancies vs U",
              xlabel  = "U",
              ylabel  = "occupancy",
              limits  = (0, 4, 0, 1.05),
              titlesize = 12)
vlines!(ax_rho, [2.46]; color = (:gray, 0.5), linewidth = 1.5, linestyle = :dash)
lines!(ax_rho, collect(U_sweep), ρ11_vec; color = :dodgerblue, linewidth = 2,
       label = "ρ₁₁  (orbital 1)")
lines!(ax_rho, collect(U_sweep), ρ22_vec; color = :crimson,    linewidth = 2,
       label = "ρ₂₂  (orbital 2)")
hlines!(ax_rho, [0.5]; color = :gray, linewidth = 1, linestyle = :dot)
axislegend(ax_rho; position = :rc, labelsize = 10)

# ── (f) effective mass m_eff vs U ─────────────────────────────────────────────
ax_m = Axis(fig[2, 3];
            title   = "(f) Effective topological mass m_eff(U)",
            xlabel  = "U",
            ylabel  = "m_eff = m + U(ρ₂₂ − ρ₁₁)/2",
            titlesize = 12)
vlines!(ax_m, [2.46];  color = (:gray, 0.5), linewidth = 1.5, linestyle = :dash)
hlines!(ax_m, [-2.0];  color = :crimson,     linewidth = 1.5, linestyle = :dash,
        label = "phase boundary  m = −2")
lines!(ax_m, collect(U_sweep), meff_vec; color = :dodgerblue, linewidth = 2,
       label = "m_eff(U)")
text!(ax_m, 2.5, -1.7; text = "C = −1", color = :dodgerblue, fontsize = 10)
text!(ax_m, 2.6, -2.7; text = "C = 0",  color = :gray,       fontsize = 10)
axislegend(ax_m; position = :lb, labelsize = 10)

Label(fig[0, :],
      "QWZ Chern insulator — Hartree-Fock renormalization from U n₁n₂ interaction  (m = −1)";
      fontsize = 14, font = :bold)

outpath = joinpath(@__DIR__, "../../output/plots/12_interacting_hf.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
