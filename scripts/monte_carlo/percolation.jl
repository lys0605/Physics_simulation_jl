using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── 2-D site percolation on the square lattice ───────────────────────────────
#
# Occupy each site independently with probability p ∈ [0, 1].  Find
# connected clusters via flood fill.  The largest cluster size, divided by
# L², is the order parameter for a continuous percolation transition.
#
#   p_c(square, site) ≈ 0.5927
#
# Outputs:
#   (a) three lattice snapshots at p = 0.5, p_c, 0.7
#   (b) order-parameter curve  ⟨largest cluster / L²⟩ vs p, for several L

# ── Flood-fill cluster labelling ──────────────────────────────────────────────
function cluster_labels(grid::BitMatrix)
    L = size(grid, 1)
    labels = zeros(Int, L, L)
    sizes  = Int[]
    next = 1
    for i in 1:L, j in 1:L
        if grid[i, j] && labels[i, j] == 0
            stack = [(i, j)]
            labels[i, j] = next
            sz = 0
            while !isempty(stack)
                ci, cj = pop!(stack)
                sz += 1
                for (di, dj) in ((1,0),(-1,0),(0,1),(0,-1))
                    ni, nj = ci+di, cj+dj
                    if 1 ≤ ni ≤ L && 1 ≤ nj ≤ L && grid[ni,nj] && labels[ni,nj] == 0
                        labels[ni, nj] = next
                        push!(stack, (ni, nj))
                    end
                end
            end
            push!(sizes, sz)
            next += 1
        end
    end
    labels, sizes
end

# ── Sweep p for several L ────────────────────────────────────────────────────
Ls       = [32, 64, 128]
ps       = range(0.4, 0.75; length=20)
n_trial  = 12
curves   = Dict{Int, Vector{Float64}}()
for L in Ls
    means = Float64[]
    for p in ps
        max_sizes = Int[]
        for _ in 1:n_trial
            g = rand(L, L) .< p
            _, sizes = cluster_labels(g)
            push!(max_sizes, isempty(sizes) ? 0 : maximum(sizes))
        end
        push!(means, mean(max_sizes) / (L*L))
    end
    curves[L] = means
end

# ── Snapshots ─────────────────────────────────────────────────────────────────
Random.seed!(7)
L_snap = 128;  p_c = 0.5927
snap_ps = [0.5, p_c, 0.7]
snaps   = NTuple{2, Any}[]
for p in snap_ps
    g = rand(L_snap, L_snap) .< p
    labels, sizes = cluster_labels(BitMatrix(g))
    # color-by-cluster (largest highlighted)
    img = zeros(Float64, L_snap, L_snap)
    if !isempty(sizes)
        biggest = argmax(sizes)
        for i in 1:L_snap, j in 1:L_snap
            img[i, j] = labels[i, j] == 0 ? 0.0 :
                        labels[i, j] == biggest ? 1.0 : 0.4
        end
    end
    push!(snaps, (p, img))
end

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(theme_light())
fig = Figure(size=(1300, 760))
Label(fig[0, :], "2-D Site Percolation on the Square Lattice  (p_c ≈ 0.5927)";
      fontsize=14, font=:bold)

ax = Axis(fig[1, 1:3]; title="Largest-cluster fraction vs p",
          xlabel="occupation probability p",
          ylabel="⟨max cluster⟩ / L²")
cols = (:dodgerblue, :crimson, :seagreen)
for (i, L) in enumerate(Ls)
    lines!(ax, ps, curves[L]; color=cols[i], linewidth=2, label="L = $L")
    scatter!(ax, ps, curves[L]; color=cols[i], markersize=7)
end
vlines!(ax, [p_c]; color=(:gray, 0.6), linestyle=:dash)
text!(ax, p_c+0.005, 0.4; text="p_c", color=:gray)
axislegend(ax; position=:lt)

for (k, (p, img)) in enumerate(snaps)
    ax2 = Axis(fig[2, k]; aspect=DataAspect(),
               title="p = $(round(p, digits=3))",
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)
    heatmap!(ax2, 1:L_snap, 1:L_snap, img;
             colormap=cgrad([:white, :lightblue, :crimson]),
             colorrange=(0, 1))
end

mkpath(joinpath(@__DIR__, "../../output/plots"))
out = joinpath(@__DIR__, "../../output/plots/percolation.png")
savefig(out, fig; px_per_unit=2)
@info "wrote $out"
