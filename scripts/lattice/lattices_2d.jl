using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie
include(joinpath(@__DIR__, "../../src/lattice/lattice_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── 2-D Bravais lattices: primitive vs conventional unit cells ────────────────
#
# Five 2-D Bravais lattices, drawn with:
#   • a grid of lattice sites
#   • primitive vectors  a₁, a₂  (arrows)
#   • the parallelogram of the primitive cell (orange)
#   • the conventional cell (cyan) — different from primitive only for the
#     centered rectangular lattice
#
# For the centered-rectangular case the conventional cell holds two lattice
# sites (one at a corner, one at the cell center), which is why the
# primitive cell is preferred for counting but the conventional cell is
# preferred for displaying the lattice's symmetry.

# ── Geometry helpers ─────────────────────────────────────────────────────────

function parallelogram(a1, a2; origin=(0.0, 0.0))
    ox, oy = origin
    [(ox, oy),
     (ox + a1[1],          oy + a1[2]),
     (ox + a1[1] + a2[1],  oy + a1[2] + a2[2]),
     (ox + a2[1],          oy + a2[2])]
end

function draw_cell!(ax, verts; color, label, linestyle=:solid, linewidth=2.5)
    xs = [v[1] for v in verts]; push!(xs, verts[1][1])
    ys = [v[2] for v in verts]; push!(ys, verts[1][2])
    lines!(ax, xs, ys; color, linewidth, linestyle, label)
end

function draw_vectors!(ax, origin, vecs, colors; lengths=nothing)
    ox, oy = origin
    for (v, c) in zip(vecs, colors)
        arrows!(ax, [ox], [oy], [v[1]], [v[2]];
                color = c, linewidth = 1.8, arrowsize = 12)
    end
end

# ── Configurations ───────────────────────────────────────────────────────────

configs = [
    (:oblique,              "Oblique"),
    (:rectangular,          "Rectangular"),
    (:centered_rectangular, "Centered Rectangular"),
    (:square,               "Square"),
    (:hexagonal,            "Hexagonal (Triangular)"),
]

# ── Figure ────────────────────────────────────────────────────────────────────

set_theme!(theme_dark())
fig = Figure(size = (1500, 920))
Label(fig[0, :], "Two-Dimensional Bravais Lattices  —  Primitive vs Conventional Cells";
      fontsize = 15, font = :bold)

for (i, (name, _)) in enumerate(configs)
    row = (i - 1) ÷ 3 + 1
    col = (i - 1) % 3 + 1
    a1, a2, A1, A2, label = bravais_2d(name)

    ax = Axis(fig[row, col];
              title  = "($("abcde"[i]))  $label",
              aspect = DataAspect(),
              titlesize = 12,
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false)

    # Lattice points
    pts = enumerate_2d(a1, a2, 4)
    xs = [p[1] for p in pts];  ys = [p[2] for p in pts]
    scatter!(ax, xs, ys; color = :white, markersize = 7)

    # Conventional cell first (so primitive overlays)
    conv = parallelogram(A1, A2)
    prim = parallelogram(a1, a2)
    same_cell = a1 == A1 && a2 == A2
    if !same_cell
        draw_cell!(ax, conv; color=(:cyan, 0.85), label="conventional",
                   linewidth=2.5)
        # extra basis site inside centered-rect conventional cell
        center = (A1 .+ A2) ./ 2
        scatter!(ax, [center[1]], [center[2]];
                 color=:cyan, markersize=10, marker=:diamond,
                 strokewidth=1, strokecolor=:white)
    end
    draw_cell!(ax, prim;
               color  = (:darkorange, 0.95),
               label  = same_cell ? "primitive = conventional" : "primitive",
               linestyle = same_cell ? :solid : :dashdot,
               linewidth = 2.5)

    # Primitive vectors arrows from origin
    arrows!(ax, [0, 0], [0, 0], [a1[1], a2[1]], [a1[2], a2[2]];
            color = :gold, linewidth = 2, arrowsize = 10)
    text!(ax, a1[1]*0.55, a1[2]*0.55 - 0.18; text="a₁",
          color=:gold, fontsize=11)
    text!(ax, a2[1]*0.55 + 0.05, a2[2]*0.55; text="a₂",
          color=:gold, fontsize=11)

    axislegend(ax; position = :rb, labelsize = 9, framevisible=true,
                bgcolor=(:black, 0.4))

    # Crop to a clean region
    span = max(norm(a1), norm(a2)) * 2.6
    xlims!(ax, -span, span); ylims!(ax, -span, span)
end

# Caption panel
ax_info = Axis(fig[2, 3];
               limits = (0, 1, 0, 1),
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)
text!(ax_info, 0.05, 0.88;
      text = "Conventions:\n" *
             "  white dots  =  lattice sites\n" *
             "  gold arrows =  primitive vectors a₁, a₂\n" *
             "  orange      =  primitive cell\n" *
             "  cyan        =  conventional cell\n\n" *
             "Five distinct point-group symmetries:\n" *
             "  C₂ (oblique), C₂ᵥ (rect, c-rect),\n" *
             "  C₄ᵥ (square), C₆ᵥ (hex)\n\n" *
             "Only the centered rectangular lattice\n" *
             "has a conventional cell ≠ primitive cell.\n" *
             "The extra site sits at the cell center.",
      color = :white, fontsize = 11, align = (:left, :top))

mkpath(joinpath(@__DIR__, "../../output/plots"))
outpath = joinpath(@__DIR__, "../../output/plots/lattices_2d.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
