using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, LinearAlgebra
include(joinpath(@__DIR__, "../../src/lattice/lattice_utils.jl"))
include(joinpath(@__DIR__, "../../src/lattice/voronoi.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Wigner-Seitz cells and Brillouin zones ────────────────────────────────────
#
# Three layers of comparison:
#
#   row 1   2-D square        |  WS (real space)  |  BZ (reciprocal)
#   row 2   2-D triangular    |  WS (real space)  |  BZ (reciprocal)
#   row 3   3-D SC            |  WS = cube        |  BZ = cube
#   row 4   3-D BCC           |  WS conv          WS prim     BZ conv     BZ prim
#                                                  (rhombic    cube       truncated
#                                                   dodecah.)              octahedron)
#   row 5   3-D FCC           |  similar set
#
# For cubic systems the "conventional" cell is the simple cube of edge a.
# Its reciprocal is also a simple cube of edge 2π/a — that gives the BZ of
# a hypothetical simple-cubic lattice with the same edge.  The physically
# relevant BZ is built from the primitive reciprocal basis.

# ── Helpers ──────────────────────────────────────────────────────────────────

function draw_polygon!(ax, poly; stroke, fill, linewidth=2.0)
    isempty(poly) && return
    xs = [v[1] for v in poly]; push!(xs, poly[1][1])
    ys = [v[2] for v in poly]; push!(ys, poly[1][2])
    poly!(ax, xs, ys; color=fill, strokecolor=stroke, strokewidth=linewidth)
end

function gather_neighbors_3d(a1, a2, a3, n)
    pts = enumerate_3d(a1, a2, a3, n)
    pts_nz = [p for p in pts if !(p[1]==0 && p[2]==0 && p[3]==0)]
    # Keep only neighbours within a radius safely larger than the WS extent
    dmin = minimum(norm.([a1, a2, a3]))
    [p for p in pts_nz if norm(collect(p)) ≤ 2.5*dmin]
end

function draw_polyhedron!(ax, vertices, faces; edgecolor=:white, fillcolor=(:cyan, 0.10),
                          linewidth=1.5)
    # Faces as transparent polygons; edges via linesegments
    for face in faces
        length(face) < 3 && continue
        xs = [vertices[i][1] for i in face]; push!(xs, vertices[face[1]][1])
        ys = [vertices[i][2] for i in face]; push!(ys, vertices[face[1]][2])
        zs = [vertices[i][3] for i in face]; push!(zs, vertices[face[1]][3])
        # Outline
        lines!(ax, xs, ys, zs; color=edgecolor, linewidth)
    end
    # Optional: scatter vertices for clarity
    if !isempty(vertices)
        vxs = [v[1] for v in vertices]
        vys = [v[2] for v in vertices]
        vzs = [v[3] for v in vertices]
        scatter!(ax, vxs, vys, vzs; color=edgecolor, markersize=4)
    end
end

# ── Figure ────────────────────────────────────────────────────────────────────

set_theme!(theme_dark())
fig = Figure(size = (1700, 1500))
Label(fig[0, :], "Wigner-Seitz Cells (real space)  and  Brillouin Zones (reciprocal space)";
      fontsize = 15, font = :bold)

# ─────────────────────── 2-D: square lattice  ────────────────────────────────

row = 1
a1, a2, _, _, _ = bravais_2d(:square)
b1, b2 = reciprocal_2d(a1, a2)

# Real-space WS
ax = Axis(fig[row, 1]; aspect=DataAspect(),
          title="(a)  Square lattice  —  WS in real space",
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false, titlesize=12)
pts = enumerate_2d(a1, a2, 3)
pts_nz = [p for p in pts if !(p[1]==0 && p[2]==0)]
scatter!(ax, [p[1] for p in pts], [p[2] for p in pts]; color=:white, markersize=7)
ws = ws_cell_2d(pts_nz)
draw_polygon!(ax, ws; stroke=:darkorange, fill=(:darkorange,0.20), linewidth=2.5)
arrows!(ax, [0,0],[0,0],[a1[1],a2[1]],[a1[2],a2[2]]; color=:gold, linewidth=1.5)
xlims!(ax, -2.5, 2.5); ylims!(ax, -2.5, 2.5)

# Reciprocal-space BZ
ax = Axis(fig[row, 2]; aspect=DataAspect(),
          title="(b)  Square — BZ = WS of reciprocal lattice",
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false, titlesize=12)
rpts = NTuple{2,Float64}[]
for i in -3:3, j in -3:3
    p = i .* b1 .+ j .* b2
    push!(rpts, (p[1], p[2]))
end
rpts_nz = [p for p in rpts if !(p[1]==0 && p[2]==0)]
scatter!(ax, [p[1] for p in rpts], [p[2] for p in rpts]; color=:white, markersize=7)
bz = ws_cell_2d(rpts_nz)
draw_polygon!(ax, bz; stroke=:cyan, fill=(:cyan,0.20), linewidth=2.5)
arrows!(ax, [0,0],[0,0],[b1[1],b2[1]],[b1[2],b2[2]]; color=:gold, linewidth=1.5)
span = 2π * 1.6
xlims!(ax, -span, span); ylims!(ax, -span, span)

# Info
ax_info = Axis(fig[row, 3]; limits=(0,1,0,1),
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)
text!(ax_info, 0.02, 0.92;
      text="Square lattice (a = 1):\n" *
           "  a₁ = x̂,  a₂ = ŷ\n" *
           "  b₁ = 2π x̂,  b₂ = 2π ŷ\n\n" *
           "Real-space WS:  square of edge 1\n" *
           "Reciprocal BZ:  square of edge 2π\n\n" *
           "Self-dual:  the reciprocal of a\n" *
           "square lattice is again square.",
      color=:white, fontsize=10, align=(:left,:top))

# ─────────────────────── 2-D: triangular lattice ─────────────────────────────

row = 2
a1, a2, _, _, _ = bravais_2d(:hexagonal)
b1, b2 = reciprocal_2d(a1, a2)

ax = Axis(fig[row, 1]; aspect=DataAspect(),
          title="(c)  Triangular lattice  —  WS = regular hexagon",
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false, titlesize=12)
pts = enumerate_2d(a1, a2, 3)
pts_nz = [p for p in pts if !(p[1]==0 && p[2]==0)]
scatter!(ax, [p[1] for p in pts], [p[2] for p in pts]; color=:white, markersize=7)
ws = ws_cell_2d(pts_nz)
draw_polygon!(ax, ws; stroke=:darkorange, fill=(:darkorange,0.20), linewidth=2.5)
arrows!(ax, [0,0],[0,0],[a1[1],a2[1]],[a1[2],a2[2]]; color=:gold, linewidth=1.5)
xlims!(ax, -2.5, 2.5); ylims!(ax, -2.5, 2.5)

ax = Axis(fig[row, 2]; aspect=DataAspect(),
          title="(d)  Triangular — BZ (rotated 30°, also hexagon)",
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false, titlesize=12)
rpts = NTuple{2,Float64}[]
for i in -3:3, j in -3:3
    p = i .* b1 .+ j .* b2
    push!(rpts, (p[1], p[2]))
end
rpts_nz = [p for p in rpts if !(p[1]==0 && p[2]==0)]
scatter!(ax, [p[1] for p in rpts], [p[2] for p in rpts]; color=:white, markersize=7)
bz = ws_cell_2d(rpts_nz)
draw_polygon!(ax, bz; stroke=:cyan, fill=(:cyan,0.20), linewidth=2.5)
arrows!(ax, [0,0],[0,0],[b1[1],b2[1]],[b1[2],b2[2]]; color=:gold, linewidth=1.5)
span = norm(b1) * 1.6
xlims!(ax, -span, span); ylims!(ax, -span, span)

ax_info = Axis(fig[row, 3]; limits=(0,1,0,1),
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)
text!(ax_info, 0.02, 0.92;
      text="Triangular lattice (a = 1):\n" *
           "  a₁ = (1, 0)\n" *
           "  a₂ = (1/2, √3/2)\n\n" *
           "Reciprocal: rotated by 30°,\n" *
           "  |b| = 4π/(a√3)\n\n" *
           "Both WS and BZ are regular hexagons\n" *
           "(triangular lattice is self-dual up to\n" *
           "rotation and rescaling).",
      color=:white, fontsize=10, align=(:left,:top))

# ─────────────────────── 3-D row helper ──────────────────────────────────────

function plot_3d_cell!(grid, vertices, faces, cell_vecs, lat_neighbors;
                       title, edgecolor=:cyan, span=nothing)
    ax = Axis3(grid;
               title, aspect = :data,
               elevation=0.18π, azimuth=0.32π,
               protrusions=(8,8,8,8),
               xlabelvisible=false, ylabelvisible=false, zlabelvisible=false,
               xticksvisible=false, yticksvisible=false, zticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false,
               zticklabelsvisible=false,
               titlesize=11)
    # Neighbouring lattice points
    nxs = [p[1] for p in lat_neighbors]
    nys = [p[2] for p in lat_neighbors]
    nzs = [p[3] for p in lat_neighbors]
    scatter!(ax, nxs, nys, nzs; color=(:white, 0.5), markersize=8)
    scatter!(ax, [0.0],[0.0],[0.0]; color=:gold, markersize=12)

    # Cell vectors arrows
    if cell_vecs !== nothing
        v1, v2, v3 = cell_vecs
        arrows!(ax, [0,0,0],[0,0,0],[0,0,0],
                [v1[1], v2[1], v3[1]],
                [v1[2], v2[2], v3[2]],
                [v1[3], v2[3], v3[3]];
                color=:gold, linewidth=1.5, arrowsize=0.04)
    end
    # Polyhedron
    draw_polyhedron!(ax, vertices, faces;
                     edgecolor, fillcolor=(edgecolor, 0.10), linewidth=2.0)

    if span !== nothing
        xlims!(ax, -span, span); ylims!(ax, -span, span); zlims!(ax, -span, span)
    end
    ax
end

# ─────────────────────── 3-D SC ──────────────────────────────────────────────

row = 3
a = 1.0
sc = bravais_3d(:sc; a)
a1, a2, a3 = sc.prim
b1, b2, b3 = reciprocal_3d(a1, a2, a3)

# Real WS (cube)
neigh_real = gather_neighbors_3d(a1, a2, a3, 2)
verts, faces = ws_cell_3d(neigh_real)
# Sample lattice points for context
sample_real = [p for p in neigh_real if norm(collect(p)) ≤ 1.2]

plot_3d_cell!(fig[row, 1], verts, faces, (a1, a2, a3), sample_real;
              title="(e)  SC — real-space WS  (cube)",
              edgecolor=:darkorange, span=0.9)

# Reciprocal BZ — same shape, different size
neigh_recip = gather_neighbors_3d(b1, b2, b3, 2)
verts_b, faces_b = ws_cell_3d(neigh_recip)
sample_b = [p for p in neigh_recip if norm(collect(p)) ≤ 1.2*norm(b1)]

plot_3d_cell!(fig[row, 2], verts_b, faces_b, (b1, b2, b3), sample_b;
              title="(f)  SC — BZ  (cube of edge 2π/a)",
              edgecolor=:cyan, span=norm(b1)*0.9)

ax_info = Axis(fig[row, 3]; limits=(0,1,0,1),
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)
text!(ax_info, 0.02, 0.92;
      text="Simple cubic (SC):\n" *
           "  conventional = primitive\n" *
           "  a₁,₂,₃ = a x̂, a ŷ, a ẑ\n" *
           "  b₁,₂,₃ = (2π/a) x̂, ŷ, ẑ\n\n" *
           "WS and BZ are both cubes.\n" *
           "Reciprocal of SC is SC scaled by 2π/a²;\n" *
           "the BZ edge is 2π/a.",
      color=:white, fontsize=10, align=(:left,:top))

# ─────────────────────── 3-D BCC ─────────────────────────────────────────────

row = 4
bcc = bravais_3d(:bcc; a)
a1, a2, a3 = bcc.prim
A1, A2, A3 = bcc.conv          # cubic edge vectors
b1, b2, b3 = reciprocal_3d(a1, a2, a3)             # primitive reciprocal (FCC)
B1, B2, B3 = reciprocal_3d(A1, A2, A3)             # conventional reciprocal (SC)

# Real WS — primitive (truncated octahedron)
neigh = gather_neighbors_3d(a1, a2, a3, 2)
verts_p, faces_p = ws_cell_3d(neigh)
sample = [p for p in neigh if norm(collect(p)) ≤ 1.4*norm(a1)]

plot_3d_cell!(fig[row, 1], verts_p, faces_p, (a1, a2, a3), sample;
              title="(g)  BCC — real-space WS (truncated octahedron)",
              edgecolor=:darkorange, span=norm(a1)*1.4)

# Conventional cubic cell — shown as wireframe with body-center site
ax_conv = Axis3(fig[row, 2]; aspect=:data,
                title="(h)  BCC — conventional cubic cell  (a = 1)",
                elevation=0.18π, azimuth=0.32π,
                protrusions=(8,8,8,8),
                xlabelvisible=false, ylabelvisible=false, zlabelvisible=false,
                xticksvisible=false, yticksvisible=false, zticksvisible=false,
                xticklabelsvisible=false, yticklabelsvisible=false,
                zticklabelsvisible=false,
                titlesize=11)
# cubic frame
let v = [(i,j,k) .* a for i in 0:1, j in 0:1, k in 0:1] |> vec,
    es = [(1,2),(1,3),(1,5),(2,4),(2,6),(3,4),(3,7),
          (4,8),(5,6),(5,7),(6,8),(7,8)]
    for (i, j) in es
        lines!(ax_conv, [v[i][1], v[j][1]], [v[i][2], v[j][2]], [v[i][3], v[j][3]];
               color=:cyan, linewidth=2)
    end
    # corners
    cx = [p[1] for p in v]; cy = [p[2] for p in v]; cz = [p[3] for p in v]
    scatter!(ax_conv, cx, cy, cz; color=:white, markersize=14)
    # body center
    scatter!(ax_conv, [0.5*a], [0.5*a], [0.5*a]; color=:gold, markersize=20)
end
xlims!(ax_conv, -0.2, 1.2); ylims!(ax_conv, -0.2, 1.2); zlims!(ax_conv, -0.2, 1.2)

# Primitive BZ (rhombic dodecahedron) and conventional BZ (cube 2π/a)
neigh_b = gather_neighbors_3d(b1, b2, b3, 2)
verts_bp, faces_bp = ws_cell_3d(neigh_b)
sample_b = [p for p in neigh_b if norm(collect(p)) ≤ 1.4*norm(b1)]

plot_3d_cell!(fig[row, 3], verts_bp, faces_bp, (b1, b2, b3), sample_b;
              title="(i)  BCC — primitive BZ (rhombic dodecahedron)",
              edgecolor=:cyan, span=norm(b1)*1.2)

# Conventional cubic BZ
neigh_B = gather_neighbors_3d(B1, B2, B3, 2)
verts_Bp, faces_Bp = ws_cell_3d(neigh_B)
sample_B = [p for p in neigh_B if norm(collect(p)) ≤ 1.2*norm(B1)]

plot_3d_cell!(fig[row, 4], verts_Bp, faces_Bp, (B1, B2, B3), sample_B;
              title="(j)  BCC — conventional cubic BZ (cube, 2π/a)",
              edgecolor=:cyan, span=norm(B1)*0.9)

# ─────────────────────── 3-D FCC ─────────────────────────────────────────────

row = 5
fcc = bravais_3d(:fcc; a)
a1, a2, a3 = fcc.prim
A1, A2, A3 = fcc.conv
b1, b2, b3 = reciprocal_3d(a1, a2, a3)             # reciprocal = BCC
B1, B2, B3 = reciprocal_3d(A1, A2, A3)             # conventional reciprocal = SC

# Real WS — primitive (rhombic dodecahedron)
neigh = gather_neighbors_3d(a1, a2, a3, 2)
verts_p, faces_p = ws_cell_3d(neigh)
sample = [p for p in neigh if norm(collect(p)) ≤ 1.4*norm(a1)]

plot_3d_cell!(fig[row, 1], verts_p, faces_p, (a1, a2, a3), sample;
              title="(k)  FCC — real-space WS (rhombic dodecahedron)",
              edgecolor=:darkorange, span=norm(a1)*1.4)

# Conventional cubic cell — wireframe with face-center sites
ax_conv = Axis3(fig[row, 2]; aspect=:data,
                title="(l)  FCC — conventional cubic cell  (a = 1)",
                elevation=0.18π, azimuth=0.32π,
                protrusions=(8,8,8,8),
                xlabelvisible=false, ylabelvisible=false, zlabelvisible=false,
                xticksvisible=false, yticksvisible=false, zticksvisible=false,
                xticklabelsvisible=false, yticklabelsvisible=false,
                zticklabelsvisible=false,
                titlesize=11)
let v = [(i,j,k) .* a for i in 0:1, j in 0:1, k in 0:1] |> vec,
    es = [(1,2),(1,3),(1,5),(2,4),(2,6),(3,4),(3,7),
          (4,8),(5,6),(5,7),(6,8),(7,8)]
    for (i, j) in es
        lines!(ax_conv, [v[i][1], v[j][1]], [v[i][2], v[j][2]], [v[i][3], v[j][3]];
               color=:cyan, linewidth=2)
    end
    cx = [p[1] for p in v]; cy = [p[2] for p in v]; cz = [p[3] for p in v]
    scatter!(ax_conv, cx, cy, cz; color=:white, markersize=14)
    faces_fcc = [t .* a for t in [(0.5,0.5,0.0),(0.5,0.5,1.0),(0.5,0.0,0.5),
                                    (0.5,1.0,0.5),(0.0,0.5,0.5),(1.0,0.5,0.5)]]
    fx = [p[1] for p in faces_fcc]; fy = [p[2] for p in faces_fcc]; fz = [p[3] for p in faces_fcc]
    scatter!(ax_conv, fx, fy, fz; color=:tomato, markersize=20)
end
xlims!(ax_conv, -0.2, 1.2); ylims!(ax_conv, -0.2, 1.2); zlims!(ax_conv, -0.2, 1.2)

# Primitive BZ — truncated octahedron
neigh_b = gather_neighbors_3d(b1, b2, b3, 2)
verts_bp, faces_bp = ws_cell_3d(neigh_b)
sample_b = [p for p in neigh_b if norm(collect(p)) ≤ 1.4*norm(b1)]

plot_3d_cell!(fig[row, 3], verts_bp, faces_bp, (b1, b2, b3), sample_b;
              title="(m)  FCC — primitive BZ (truncated octahedron)",
              edgecolor=:cyan, span=norm(b1)*1.2)

# Conventional cubic BZ
neigh_B = gather_neighbors_3d(B1, B2, B3, 2)
verts_Bp, faces_Bp = ws_cell_3d(neigh_B)
sample_B = [p for p in neigh_B if norm(collect(p)) ≤ 1.2*norm(B1)]

plot_3d_cell!(fig[row, 4], verts_Bp, faces_Bp, (B1, B2, B3), sample_B;
              title="(n)  FCC — conventional cubic BZ (cube, 2π/a)",
              edgecolor=:cyan, span=norm(B1)*0.9)

# Final note row
ax_info = Axis(fig[6, :]; limits=(0,1,0,1),
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)
text!(ax_info, 0.01, 0.92;
      text = "Reciprocal duality of cubic lattices:    SC ↔ SC,    BCC ↔ FCC,    FCC ↔ BCC\n" *
             "The physically meaningful Brillouin zone is the WS cell of the *primitive* reciprocal lattice.\n" *
             "Using the conventional cubic basis (with its larger cell) gives a smaller reciprocal cube — \n" *
             "useful for band-folding diagrams but not the true first Brillouin zone.",
      color=:white, fontsize=11, align=(:left,:top))
rowsize!(fig.layout, 6, Fixed(80))

mkpath(joinpath(@__DIR__, "../../output/plots"))
outpath = joinpath(@__DIR__, "../../output/plots/brillouin.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
