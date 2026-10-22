using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, LinearAlgebra
include(joinpath(@__DIR__, "../../src/lattice/lattice_utils.jl"))
include(joinpath(@__DIR__, "../../src/lattice/voronoi.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Voronoi diagrams and their condensed-matter use ───────────────────────────
#
# A Voronoi diagram of a discrete point set S partitions space into cells
#   V(p) = { x : |x − p| ≤ |x − q|  ∀ q ∈ S \ {p} }
#
# Each cell is convex; the cell boundary between p and q is a segment of the
# perpendicular bisector of pq.  When S is a Bravais lattice, every cell is
# congruent and centered on a lattice site — the Wigner-Seitz (WS) cell.
#
# This script demonstrates the algorithm and the WS-Voronoi connection:
#
#   (a) Random points: irregular Voronoi cells of varying shapes.
#   (b) Square lattice:   WS = square (rotated 0°).
#   (c) Triangular lattice: WS = regular hexagon.
#   (d) Centered rectangular: WS is hexagon whose proportions depend on the
#       primitive vectors — useful for showing that WS shape is a function
#       of the lattice geometry, not its symbol.

# ── Common drawing helper ─────────────────────────────────────────────────────

function draw_cell!(ax, poly; color, fill=(color, 0.25), linewidth=1.6)
    isempty(poly) && return
    xs = [v[1] for v in poly]; push!(xs, poly[1][1])
    ys = [v[2] for v in poly]; push!(ys, poly[1][2])
    poly!(ax, xs, ys; color=fill, strokecolor=color, strokewidth=linewidth)
end

# ── Panel layout ──────────────────────────────────────────────────────────────

set_theme!(theme_dark())
fig = Figure(size = (1500, 980))
Label(fig[0, :], "Voronoi Diagrams  —  General Point Sets and the Wigner-Seitz Cell";
      fontsize = 15, font = :bold)

# ── (a) Voronoi of random points ──────────────────────────────────────────────
ax_a = Axis(fig[1, 1];
            title="(a)  Voronoi tessellation of 30 random points",
            aspect=DataAspect(),
            limits=(-5, 5, -5, 5),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false,
            titlesize=12)

Random.seed!(7)
N_rand = 30
random_pts = [(rand()*8 - 4, rand()*8 - 4) for _ in 1:N_rand]
random_pts_real = [(Float64(p[1]), Float64(p[2])) for p in random_pts]

for p in random_pts_real
    cell = voronoi_cell_2d(p, random_pts_real; bbox=(-5.0, 5.0, -5.0, 5.0))
    draw_cell!(ax_a, cell; color=:cyan, fill=(:cyan, 0.08), linewidth=1.2)
end
scatter!(ax_a, [p[1] for p in random_pts_real], [p[2] for p in random_pts_real];
         color=:gold, markersize=8)

# ── (b–d) Lattice WS cells via Voronoi ────────────────────────────────────────

function wigner_seitz_panel!(ax, name; n_neighbors=4, span_factor=2.5)
    a1, a2, _, _, label = bravais_2d(name)
    # Generate lattice neighbors (origin excluded)
    pts = enumerate_2d(a1, a2, n_neighbors)
    pts_nz = [p for p in pts if !(p[1] == 0 && p[2] == 0)]

    span = max(norm(a1), norm(a2)) * span_factor
    xlims!(ax, -span, span);  ylims!(ax, -span, span)

    # Draw all Voronoi cells (each lattice site gets a WS cell)
    for p in pts
        neighbors = [q for q in pts if !(q[1] == p[1] && q[2] == p[2])]
        cell = voronoi_cell_2d(p, neighbors; bbox=(-2span, 2span, -2span, 2span))
        draw_cell!(ax, cell; color=(:cyan, 0.55), fill=(:cyan, 0.05), linewidth=1.2)
    end

    # Highlight the cell at the origin
    central = ws_cell_2d(pts_nz; bbox=(-2span, 2span, -2span, 2span))
    draw_cell!(ax, central; color=:darkorange, fill=(:darkorange, 0.20), linewidth=2.5)

    # Lattice points
    scatter!(ax, [p[1] for p in pts], [p[2] for p in pts];
             color=:white, markersize=7)
    scatter!(ax, [0.0], [0.0]; color=:gold, markersize=11)

    # Primitive vector arrows
    arrows!(ax, [0,0], [0,0], [a1[1], a2[1]], [a1[2], a2[2]];
            color=:gold, linewidth=1.5, arrowsize=10)
    label
end

ax_b = Axis(fig[1, 2]; aspect=DataAspect(),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false, titlesize=12)
lab_b = wigner_seitz_panel!(ax_b, :square)
ax_b.title[] = "(b)  Square lattice  →  WS cell = square"

ax_c = Axis(fig[1, 3]; aspect=DataAspect(),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false, titlesize=12)
lab_c = wigner_seitz_panel!(ax_c, :hexagonal)
ax_c.title[] = "(c)  Triangular lattice  →  WS = regular hexagon"

ax_d = Axis(fig[2, 1]; aspect=DataAspect(),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false, titlesize=12)
lab_d = wigner_seitz_panel!(ax_d, :centered_rectangular)
ax_d.title[] = "(d)  Centered rectangular  →  WS = elongated hexagon"

# ── (e) Bisector-construction zoom ────────────────────────────────────────────
ax_e = Axis(fig[2, 2];
            title="(e)  How the WS cell is built  —  perpendicular bisectors",
            aspect=DataAspect(),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false,
            titlesize=12)

a1, a2, _, _, _ = bravais_2d(:hexagonal)
pts_e = enumerate_2d(a1, a2, 2)
pts_e_nz = [p for p in pts_e if !(p[1]==0 && p[2]==0)]

# 6 nearest neighbours (hexagonal)
neighbors_sorted = sort(pts_e_nz; by = q -> q[1]^2 + q[2]^2)
nearest = neighbors_sorted[1:6]

# Draw bisectors
for q in nearest
    # midpoint
    mx, my = q[1]/2, q[2]/2
    # perpendicular direction
    px, py = -q[2], q[1];  L = sqrt(px^2 + py^2);  px /= L; py /= L
    lines!(ax_e, [mx - 1.3*px, mx + 1.3*px], [my - 1.3*py, my + 1.3*py];
           color=(:cyan, 0.65), linewidth=1.2, linestyle=:dash)
    # connecting segment to neighbour
    lines!(ax_e, [0.0, q[1]], [0.0, q[2]];
           color=(:white, 0.35), linewidth=0.9)
end

# WS cell
ws = ws_cell_2d(pts_e_nz; bbox=(-5.0, 5.0, -5.0, 5.0))
draw_cell!(ax_e, ws; color=:darkorange, fill=(:darkorange, 0.25), linewidth=2.5)

scatter!(ax_e, [p[1] for p in pts_e], [p[2] for p in pts_e];
         color=:white, markersize=8)
scatter!(ax_e, [0.0], [0.0]; color=:gold, markersize=12)

xlims!(ax_e, -1.6, 1.6);  ylims!(ax_e, -1.6, 1.6)

# ── (f) Notes ─────────────────────────────────────────────────────────────────
ax_info = Axis(fig[2, 3];
               limits=(0,1,0,1),
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)
text!(ax_info, 0.02, 0.98;
      text = "Construction (Wigner-Seitz)\n\n" *
             "  1. From site R₀, list neighbours {Rᵢ}.\n" *
             "  2. For each Rᵢ, the perpendicular\n" *
             "     bisector of  R₀-Rᵢ  is a half-plane\n" *
             "     boundary:  Rᵢ · x ≤ |Rᵢ|²/2  with\n" *
             "     R₀ placed at origin.\n" *
             "  3. WS = ⋂ᵢ (these half-planes).\n" *
             "  4. Only the nearest few shells contribute;\n" *
             "     in 2-D triangular lattice, only the\n" *
             "     6 nearest neighbours are needed.\n\n" *
             "Each WS cell:\n" *
             "  • holds exactly one lattice site\n" *
             "  • tiles space without gaps or overlap\n" *
             "  • is centro-symmetric\n" *
             "  • is the unique most-symmetric primitive cell",
      color=:white, fontsize=10, align=(:left, :top))

mkpath(joinpath(@__DIR__, "../../output/plots"))
outpath = joinpath(@__DIR__, "../../output/plots/voronoi.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
