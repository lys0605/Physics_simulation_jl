using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, LinearAlgebra
include(joinpath(@__DIR__, "../../src/greens/qwz_model.jl"))
include(joinpath(@__DIR__, "../../src/greens/surface_greens.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Semi-infinite boundary condition: surface Green's function ─────────────────
#
# Semi-infinite geometry: periodic BC in y, open (hard-wall) at x = 0, and
# the system extends to x → +∞.  There is no finite-size truncation.
#
# The surface spectral function
#
#   A_s(k_y, ω) = −(1/π) Im Tr G_s(k_y, ω)
#
# isolates the boundary contribution.  For a topological system it shows
# sharply defined in-gap edge modes — their "hardness" (no finite-size
# hybridisation gap) distinguishes them from the ribbon edge states.
#
# Panels:
#   (a)  A_s(k_y, ω) — topological phase m = −1
#   (b)  A_s(k_y, ω) — trivial phase m = 3
#   (c)  ‖α‖ per López Sancho iteration for two η values
#   (d)  A_s vs ω at k_y = 0 : comparison naive vs López Sancho

const Nky = 250
const Nω  = 250
const η   = 0.025
const kys = range(-π, π; length = Nky)
const ωs  = range(-3.5, 3.5; length = Nω)

m_top  = QWZModel(-1.0)
m_triv = QWZModel( 3.0)

# ── compute surface spectral maps ─────────────────────────────────────────────
@info "Surface GF — topological phase…"
A_top  = surface_spectral_map(m_top,  kys, ωs; η)
@info "Surface GF — trivial phase…"
A_triv = surface_spectral_map(m_triv, kys, ωs; η)

# ── convergence curve ─────────────────────────────────────────────────────────
conv_small = lopezSancho_convergence(m_top, 0.0, 0.0; η=0.02)
conv_large = lopezSancho_convergence(m_top, 0.0, 0.0; η=0.2)

# ── naive vs López Sancho at k_y = 0 ─────────────────────────────────────────
@info "Naive iteration…"
A_naive = [surface_spectral(m_top, 0.0, ω; η, tol=1e-6,
                             max_iter=10_000,
                             # dispatch to naive
                             ) for ω in ωs]

# Reuse surface_spectral (which calls López Sancho) for comparison
A_ls = [begin
    G = surface_green(m_top, 0.0, ω; η)
    -imag(tr(G)) / π
end for ω in ωs]

# naive reference
A_naive2 = [begin
    G = surface_green_naive(m_top, 0.0, ω; η=0.1, max_iter=500)
    -imag(tr(G)) / π
end for ω in ωs]

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1300, 780))

# (a) topological surface spectral function
ax_a = Axis(fig[1, 1];
            title  = "(a) Topological  m = −1  (C = +1)",
            xlabel = "k_y",
            ylabel = "ω",
            xticks = ([-π, 0, π], ["-π", "0", "π"]),
            limits = (-π, π, -3.5, 3.5),
            titlesize = 12)
heatmap!(ax_a, kys, ωs, A_top; colormap = :hot, colorrange = (0, 2))
text!(ax_a, 0.0, 0.05; text = "edge modes", color = :cyan,
      align = (:center, :bottom), fontsize = 11)

# (b) trivial surface spectral function
ax_b = Axis(fig[1, 2];
            title  = "(b) Trivial  m = 3  (C = 0)",
            xlabel = "k_y",
            ylabel = "ω",
            xticks = ([-π, 0, π], ["-π", "0", "π"]),
            limits = (-π, π, -3.5, 3.5),
            titlesize = 12)
heatmap!(ax_b, kys, ωs, A_triv; colormap = :hot, colorrange = (0, 2))

# (c) convergence of ‖α‖
ax_c = Axis(fig[2, 1];
            title   = "(c) López Sancho convergence — ‖α_n‖ per iteration",
            xlabel  = "iteration n",
            ylabel  = "‖α_n‖",
            yscale  = log10,
            titlesize = 12)
lines!(ax_c, 1:length(conv_small), conv_small; color = :dodgerblue, linewidth=2,
       label = "η = 0.02")
lines!(ax_c, 1:length(conv_large), conv_large; color = :crimson,    linewidth=2,
       label = "η = 0.2")
hlines!(ax_c, [1e-8]; color = :gray, linestyle = :dash, linewidth = 1,
        label = "tolerance")
axislegend(ax_c; position = :rt)

# (d) naive vs López Sancho slice at k_y = 0
ax_d = Axis(fig[2, 2];
            title   = "(d) A_s(k_y=0, ω) : naive vs López Sancho  (η = 0.025 / 0.1)",
            xlabel  = "ω",
            ylabel  = "A_s",
            titlesize = 12)
lines!(ax_d, ωs, A_ls;     color = :dodgerblue, linewidth = 2,   label = "López Sancho η=0.025")
lines!(ax_d, ωs, A_naive2; color = :crimson,    linewidth = 1.5, label = "Naive iter. η=0.1",
       linestyle = :dash)
axislegend(ax_d; position = :rt)

Label(fig[0, :], "Surface Green's function — semi-infinite BC in x, periodic in y";
      fontsize = 14, font = :bold)

outpath = joinpath(@__DIR__, "../../output/plots/11_surface_greens.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
