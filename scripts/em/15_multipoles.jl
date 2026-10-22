using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie
include(joinpath(@__DIR__, "../../src/em/em_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Electrostatic multipoles ──────────────────────────────────────────────────
#
# Four panels showing the electric potential (heat map) and field lines
# for successively higher-order multipoles:
#   (a) Monopole    (+1 charge)
#   (b) Dipole      (+1, −1) separated along y
#   (c) Quadrupole  four charges with zero monopole and dipole moment
#   (d) Octupole    eight charges with zero monopole, dipole, and quadrupole

# ── Charge configurations ─────────────────────────────────────────────────────

d = 0.7   # charge separation (half-distance)

configs = [
    # (label, charges, positions)
    ("(a)  Monopole",
     [1.0],
     [(0.0, 0.0)]),

    ("(b)  Dipole",
     [1.0, -1.0],
     [(0.0,  d), (0.0, -d)]),

    ("(c)  Quadrupole",
     [1.0, 1.0, -1.0, -1.0],
     [( d, 0.0), (-d, 0.0), (0.0,  d), (0.0, -d)]),

    ("(d)  Octupole",
     [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0],
     [(d*cos(k*π/4), d*sin(k*π/4)) for k in 0:7]),
]

# ── Grid ──────────────────────────────────────────────────────────────────────

N  = 300
L  = 2.0
xs = range(-L, L; length=N)
ys = range(-L, L; length=N)

# ── Figure ────────────────────────────────────────────────────────────────────

set_theme!(theme_dark())
fig = Figure(size = (1200, 1100))
Label(fig[0, :], "Electrostatic Multipoles — Potential and Field Lines";
      fontsize = 15, font = :bold)

colsize!(fig.layout, 1, Relative(0.5))
colsize!(fig.layout, 2, Relative(0.5))

for (idx, (label, qs, pos)) in enumerate(configs)
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1

    ax = Axis(fig[row, col];
              title   = label,
              aspect  = DataAspect(),
              limits  = (-L, L, -L, L),
              xticksvisible = false, yticksvisible = false,
              xticklabelsvisible = false, yticklabelsvisible = false,
              titlesize = 13)

    # Potential heatmap
    V = [potential_2d(qs, pos, x, y; reg=0.04) for x in xs, y in ys]

    # Clip extremes so color scale isn't dominated by near-field divergence
    Vclip = clamp.(V, -5.0, 5.0)

    heatmap!(ax, xs, ys, Vclip;
             colormap = :RdBu_r, colorrange = (-5, 5), interpolate = true)

    # Field lines via streamplot
    function field!(pos_in, q_in, p, t)
        Ex, Ey = efield_2d(q_in, pos_in, p[1], p[2]; reg=0.06)
        Point2f(Ex, Ey)
    end
    streamplot!(ax,
                (p, t) -> field!(pos, qs, p, t),
                xs, ys;
                colormap  = [:white],
                linewidth = 1.0,
                density   = 1.0,
                arrow_size = 7)

    # Mark charge positions
    for (q, (px, py)) in zip(qs, pos)
        color  = q > 0 ? :red : :blue
        marker = q > 0 ? :circle : :xcross
        scatter!(ax, [px], [py]; color, marker, markersize = 14, strokewidth = 1,
                 strokecolor = :white)
    end
end

# ── Save ──────────────────────────────────────────────────────────────────────

mkpath(joinpath(@__DIR__, "../../output/plots"))
outpath = joinpath(@__DIR__, "../../output/plots/15_multipoles.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
