using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, LinearAlgebra
include(joinpath(@__DIR__, "../../src/lattice/lattice_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── 3-D Bravais lattices: simple cubic, BCC, FCC ──────────────────────────────
#
# Each row shows one lattice with two views:
#   (left)  conventional cubic cell — natural axes, multiple sites inside
#   (right) primitive cell          — single site, parallelepiped shape
#
# Sites are drawn shared-correctly: corners with their occupancy fraction
# (corner = 1/8, face = 1/2, body = 1).  The total count inside the
# conventional cube equals 1 (SC), 2 (BCC), 4 (FCC).

# ── Drawing helpers ──────────────────────────────────────────────────────────

const CUBE_EDGES = let
    v = [(0,0,0),(1,0,0),(1,1,0),(0,1,0),
         (0,0,1),(1,0,1),(1,1,1),(0,1,1)]
    edges = [(1,2),(2,3),(3,4),(4,1),
             (5,6),(6,7),(7,8),(8,5),
             (1,5),(2,6),(3,7),(4,8)]
    [(v[i], v[j]) for (i,j) in edges]
end

function draw_box!(ax, A1, A2, A3; color, linewidth=1.8, linestyle=:solid)
    # Eight corners of the parallelepiped spanned by A1,A2,A3
    corners = [i.*A1 .+ j.*A2 .+ k.*A3 for i in 0:1, j in 0:1, k in 0:1]
    corners = vec(corners)
    edges = [(1,2),(1,3),(1,5),(2,4),(2,6),(3,4),(3,7),
             (4,8),(5,6),(5,7),(6,8),(7,8)]
    pts = NTuple{2,NTuple{3,Float64}}[]
    for (i, j) in edges
        push!(pts, (Tuple(corners[i]), Tuple(corners[j])))
    end
    for (p, q) in pts
        lines!(ax, [p[1], q[1]], [p[2], q[2]], [p[3], q[3]];
               color, linewidth, linestyle)
    end
end

function lattice_sites_in_box(a1, a2, a3, n=3)
    sites = NTuple{3,Float64}[]
    for i in -n:n, j in -n:n, k in -n:n
        p = i .* a1 .+ j .* a2 .+ k .* a3
        push!(sites, (p[1], p[2], p[3]))
    end
    sites
end

# Show only sites inside (or on the boundary of) the conventional box [0,a]³
function sites_in_cube(sites, a; tol=1e-6)
    [p for p in sites if -tol ≤ p[1] ≤ a+tol &&
                          -tol ≤ p[2] ≤ a+tol &&
                          -tol ≤ p[3] ≤ a+tol]
end

# ── Figure ────────────────────────────────────────────────────────────────────

set_theme!(theme_dark())
fig = Figure(size = (1300, 1200))
Label(fig[0, :], "Three-Dimensional Cubic Bravais Lattices  —  Conventional vs Primitive Cells";
      fontsize = 15, font = :bold)

a = 1.0
configs = [(:sc, "Simple Cubic (SC)"),
           (:bcc, "Body-Centered Cubic (BCC)"),
           (:fcc, "Face-Centered Cubic (FCC)")]

for (row_idx, (name, label)) in enumerate(configs)
    info = bravais_3d(name; a)
    (a1, a2, a3)   = info.prim
    (A1, A2, A3)   = info.conv

    # Sites in conventional cell (corners + face/body centers)
    conv_sites = [Float64.(b) .* a for b in info.basis_conv]   # corner site
    # Add the 7 other corners of the cubic cell
    corners = [Float64.((i,j,k)) .* a for i in 0:1, j in 0:1, k in 0:1] |> vec
    # Faces (FCC)
    face_sites = NTuple{3,Float64}[]
    if name === :fcc
        face_sites = [t .* a for t in [(0.5, 0.5, 0.0), (0.5, 0.5, 1.0),
                                         (0.5, 0.0, 0.5), (0.5, 1.0, 0.5),
                                         (0.0, 0.5, 0.5), (1.0, 0.5, 0.5)]]
    end
    body = name === :bcc ? [(0.5, 0.5, 0.5) .* a] : NTuple{3,Float64}[]

    # ── Conventional view ────────────────────────────────────────────────────
    ax_c = Axis3(fig[row_idx, 1];
                 title = "($("abc"[row_idx]).1)  $label — conventional cell",
                 aspect = :data,
                 elevation = 0.18π, azimuth = 0.32π,
                 protrusions = (8,8,8,8),
                 xlabelvisible=false, ylabelvisible=false, zlabelvisible=false,
                 xticksvisible=false, yticksvisible=false, zticksvisible=false,
                 xticklabelsvisible=false, yticklabelsvisible=false,
                 zticklabelsvisible=false,
                 titlesize = 12)

    # Draw cubic box
    draw_box!(ax_c, [a,0,0], [0,a,0], [0,0,a]; color=:cyan, linewidth=2)

    # Corner sites (small)
    cxs = [p[1] for p in corners]; cys = [p[2] for p in corners]; czs = [p[3] for p in corners]
    scatter!(ax_c, cxs, cys, czs; color=:white, markersize=18, strokewidth=1.5,
             strokecolor=:black)

    # Face sites (FCC)
    if !isempty(face_sites)
        fxs = [p[1] for p in face_sites]; fys = [p[2] for p in face_sites]; fzs = [p[3] for p in face_sites]
        scatter!(ax_c, fxs, fys, fzs; color=:tomato, markersize=22, strokewidth=1.5,
                 strokecolor=:black)
    end
    # Body site (BCC)
    if !isempty(body)
        bxs = [p[1] for p in body]; bys = [p[2] for p in body]; bzs = [p[3] for p in body]
        scatter!(ax_c, bxs, bys, bzs; color=:gold, markersize=24, strokewidth=1.5,
                 strokecolor=:black)
    end

    # Lattice axes
    arrows!(ax_c, [0,0,0], [0,0,0], [0,0,0],
            [a*1.15, 0, 0], [0, a*1.15, 0], [0, 0, a*1.15];
            color = :gold, linewidth = 1.5, arrowsize = 0.05)

    span = a * 1.4
    xlims!(ax_c, -0.2, span); ylims!(ax_c, -0.2, span); zlims!(ax_c, -0.2, span)

    # ── Primitive view ───────────────────────────────────────────────────────
    ax_p = Axis3(fig[row_idx, 2];
                 title = "($("abc"[row_idx]).2)  $label — primitive cell",
                 aspect = :data,
                 elevation = 0.18π, azimuth = 0.32π,
                 protrusions = (8,8,8,8),
                 xlabelvisible=false, ylabelvisible=false, zlabelvisible=false,
                 xticksvisible=false, yticksvisible=false, zticksvisible=false,
                 xticklabelsvisible=false, yticklabelsvisible=false,
                 zticklabelsvisible=false,
                 titlesize = 12)

    # Ghost cubic frame for reference
    draw_box!(ax_p, [a,0,0], [0,a,0], [0,0,a];
              color=(:cyan, 0.25), linewidth=1.0, linestyle=:dot)

    # Primitive parallelepiped
    draw_box!(ax_p, a1, a2, a3; color=:darkorange, linewidth=2.5)

    # Lattice sites nearby
    sites = lattice_sites_in_box(a1, a2, a3, 3)
    keep  = [p for p in sites if -0.4 ≤ p[1] ≤ a*1.4 &&
                                  -0.4 ≤ p[2] ≤ a*1.4 &&
                                  -0.4 ≤ p[3] ≤ a*1.4]
    sxs = [p[1] for p in keep]; sys = [p[2] for p in keep]; szs = [p[3] for p in keep]
    scatter!(ax_p, sxs, sys, szs; color=:white, markersize=16, strokewidth=1.5,
             strokecolor=:black)

    # Primitive vectors as arrows from the origin
    arrows!(ax_p, [0,0,0], [0,0,0], [0,0,0],
            [a1[1], a2[1], a3[1]],
            [a1[2], a2[2], a3[2]],
            [a1[3], a2[3], a3[3]];
            color = :gold, linewidth = 2, arrowsize = 0.04)

    xlims!(ax_p, -0.2, span); ylims!(ax_p, -0.2, span); zlims!(ax_p, -0.2, span)
end

# Caption panel
ax_info = Axis(fig[4, 1:2];
               limits = (0, 1, 0, 1),
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)
text!(ax_info, 0.02, 0.98;
      text = "Primitive vectors  (a is the conventional cubic edge)\n" *
             "  SC :  a₁ = a x̂,   a₂ = a ŷ,   a₃ = a ẑ          (V_p = a³)\n" *
             "  BCC:  a₁ = (a/2)(−x̂+ŷ+ẑ),  a₂ = (a/2)(x̂−ŷ+ẑ),  a₃ = (a/2)(x̂+ŷ−ẑ)   (V_p = a³/2)\n" *
             "  FCC:  a₁ = (a/2)(ŷ+ẑ),     a₂ = (a/2)(x̂+ẑ),    a₃ = (a/2)(x̂+ŷ)        (V_p = a³/4)\n\n" *
             "Number of lattice sites in the conventional cubic cell:\n" *
             "  SC = 1   (corners only:  8 × ⅛)\n" *
             "  BCC = 2  (corners + body center)\n" *
             "  FCC = 4  (corners + 6 × ½ face sites)",
      color=:white, fontsize=10, align=(:left, :top), font=:regular)
rowsize!(fig.layout, 4, Fixed(140))

mkpath(joinpath(@__DIR__, "../../output/plots"))
outpath = joinpath(@__DIR__, "../../output/plots/lattices_3d.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
