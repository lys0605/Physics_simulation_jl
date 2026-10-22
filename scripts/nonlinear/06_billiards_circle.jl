using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, DynamicalBilliards
include(joinpath(@__DIR__, "../../src/nonlinear/billiards.jl"))

# ── Integrable billiard: 100-gon approximation of the circular billiard ───────
#
# In the circular billiard the tangential velocity component p = sin α is a
# conserved quantity (angular momentum / R).  Each orbit stays on a horizontal
# line in the Poincaré section (s, p) — the signature of an integrable system.
#
# Three orbits with increasing |p| are evolved simultaneously.  In the left
# panel each chord is drawn as it is added; the accumulation of chords reveals
# the "caustic" — an inner circle tangent to every chord.  The right panel
# shows the Poincaré map building up one point per bounce; the near-horizontal
# scatter of each orbit confirms integrability.

bd = make_circle(1.0; n_sides = 100)
bdy = boundary_outline(bd)

N  = 300        # collisions per orbit
φs = (1.1, 0.7, 0.35)   # launch angles → three distinct rotation numbers
clrs = (:dodgerblue, :crimson, :seagreen)

# Pre-compute trajectories and boundary maps (DynamicalBilliards mutates the
# particle in place, so we create fresh Particle objects for each call).
pos0 = (0.3, 0.0)
trajs = [begin
           p = Particle(pos0..., φ)
           ct, poss, _ = evolve!(p, bd, N)
           poss
         end for φ in φs]

bmaps = [begin
           p = Particle(pos0..., φ)
           bm, _ = boundarymap(p, bd, N)
           bm
         end for φ in φs]

_, ints = let p = Particle(pos0..., φs[1])
    boundarymap(p, bd, 1)
end

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1100, 520))

ax_bd = Axis(fig[1, 1];
             title  = "Circular billiard — invariant caustics",
             xlabel = "x", ylabel = "y",
             aspect = DataAspect(),
             limits = (-1.25, 1.25, -1.25, 1.25))

ax_ps = Axis(fig[1, 2];
             title  = "Poincaré section  (s,  sin α)",
             xlabel = "arc-length  s",
             ylabel = "sin α",
             limits = (0.0, ints[end], -1.05, 1.05))

# static boundary
lines!(ax_bd, bdy; color = :black, linewidth = 2)

# Observables: one linesegments series per orbit (chord accumulation)
segs  = [Observable(Point2f[]) for _ in φs]
pss   = [Observable(Point2f[]) for _ in φs]
balls = [Observable(Point2f[poss[1]]) for poss in trajs]

labels = ["orbit 1", "orbit 2", "orbit 3"]
for (i, (clr, lbl)) in enumerate(zip(clrs, labels))
    linesegments!(ax_bd, segs[i];  color = (clr, 0.45), linewidth = 1.0)
    scatter!(ax_bd,     balls[i];  color = clr,         markersize = 10)
    scatter!(ax_ps,     pss[i];    color = (clr, 0.7),  markersize = 3, label = lbl)
end
axislegend(ax_ps; position = :rb, framevisible = true)

outpath = joinpath(@__DIR__, "../../output/videos/06_billiards_circle.mp4")
mkpath(dirname(outpath))

record(fig, outpath, 1:N; framerate = 24) do k
    for i in 1:3
        a = Point2f(trajs[i][k])
        b = Point2f(trajs[i][k+1])
        push!(segs[i][], a, b);  notify(segs[i])
        balls[i][] = [b]
        push!(pss[i][], Point2f(bmaps[i][k]));  notify(pss[i])
    end
end

@info "wrote $outpath"
