using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie
include(joinpath(@__DIR__, "../../src/greens/qwz_model.jl"))
include(joinpath(@__DIR__, "../../src/greens/bulk_greens.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Bulk spectral function of the QWZ Chern insulator ─────────────────────────
#
# Periodic BC in both x and y → k is a good quantum number.
# We compute:
#   (a) Band structure E±(k) along Γ→X→M→Γ
#   (b) Spectral function A(k,ω) = −Im Tr G(k,ω)/π  (heatmap)
#   (c) Density of states ρ(ω) (integrated over BZ)
#   (d) Chern number C (printed; shown as annotation)
#
# Topological phase: m = −1  (Chern number C = +1)

const m_top  = -1.0
const m_triv = +3.0
const η      = 0.04

models = [QWZModel(m_top), QWZModel(m_triv)]
titles = ["Topological  (m = $m_top,  C = $(chern_number(QWZModel(m_top))))",
          "Trivial       (m = $m_triv,  C = $(chern_number(QWZModel(m_triv))))"]

# ── k-path and ω grid ─────────────────────────────────────────────────────────
kxs, kys, tick_pos, tick_lbl = kpath_points(70)
Nk   = length(kxs)
ωs   = range(-4.0, 4.0; length=300)
k_idx = 1:Nk

# ── precompute ────────────────────────────────────────────────────────────────
function compute_all(model)
    # Band energies along path
    E_lower = [band_energies(model, kxs[i], kys[i])[1] for i in k_idx]
    E_upper = [band_energies(model, kxs[i], kys[i])[2] for i in k_idx]
    # Spectral function heatmap  (Nk × Nω)
    A = [spectral_function_fast(model, kxs[i], kys[i], ω; η)
         for i in k_idx, ω in ωs]
    # DOS
    ρ = dos(model, ωs; Nk=120, η)
    E_lower, E_upper, A, ρ
end

@info "Computing spectral functions…"
results = map(compute_all, models)

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1400, 640))

for (col, (title, (E₋, E₊, A, ρ))) in enumerate(zip(titles, results))
    # ── band structure + spectral heatmap (shared k-axis)
    ga  = fig[1, 2col-1] = GridLayout()

    ax_A = Axis(ga[1, 1];
                title       = title,
                ylabel      = "ω",
                limits      = (1, Nk, -4.0, 4.0),
                xticks      = (tick_pos, tick_lbl),
                xticklabelsize = 13,
                yticklabelsize = 11,
                titlesize   = 12)

    heatmap!(ax_A, k_idx, ωs, A';
             colormap = :hot, colorrange = (0, 2.5))
    lines!(ax_A, k_idx, E₋; color = :dodgerblue,  linewidth = 1.5, linestyle = :dash)
    lines!(ax_A, k_idx, E₊; color = :dodgerblue,  linewidth = 1.5, linestyle = :dash)

    # Vertical lines at high-symmetry points
    for tp in tick_pos
        vlines!(ax_A, tp; color = (:white, 0.4), linewidth = 1)
    end

    # ── DOS panel to the right (shared ω-axis)
    ax_D = Axis(ga[1, 2];
                ylabel      = "ω",
                xlabel      = "DOS",
                yticklabelsvisible = false,
                limits      = (0, nothing, -4.0, 4.0),
                width       = 90)
    linkaxes!(ax_A, ax_D)
    lines!(ax_D, ρ, collect(ωs); color = :orangered, linewidth = 2)
    hlines!(ax_D, [0.0]; color = :gray, linewidth = 0.8, linestyle = :dot)

    colsize!(ga, 2, Fixed(90))
end

Label(fig[0, :], "QWZ bulk spectral function  A(k,ω) = −Im Tr G(k,ω) / π";
      fontsize = 15, font = :bold)

outpath = joinpath(@__DIR__, "../../output/plots/09_bulk_spectral.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
