using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, DynamicalBilliards, Statistics
include(joinpath(@__DIR__, "../../src/nonlinear/billiards.jl"))

# ── Bunimovich stadium — chaos, ergodicity, and sensitive dependence ──────────
#
# The stadium billiard (two semicircular caps + two straight walls) is one of
# the first rigorously proven chaotic billiards (Bunimovich 1979).  Its
# positive Lyapunov exponent means that nearby trajectories diverge
# exponentially.  There are no invariant curves except at the boundary of phase
# space (p = ±1, grazing orbits) — the system is ergodic and mixing.
#
# This script produces two panels:
#   Left  — a single long trajectory; the stadium fills uniformly (ergodic)
#   Right — two particles launched with Δφ = 10⁻⁵ rad; their separation grows
#           exponentially and they become visually uncorrelated within ~30 bounces.

bd  = billiard_stadium()
bdy = boundary_outline(bd)

# Compute the positive Lyapunov exponent (averaged over many orbits).
λ_avg = mean(lyapunovspectrum(randominside(bd), bd, 5_000)[1] for _ in 1:10)
@info "Bunimovich stadium  λ₁ ≈ $(round(λ_avg; digits=3))"

# ── ergodic trajectory (left panel) ──────────────────────────────────────────
N_long = 600
p_erg = randominside(bd)
_, poss_erg, _ = evolve!(p_erg, bd, N_long)
pts_erg = dense_trajectory(poss_erg; pts_per_seg = 8)   # 4800 path points

# ── divergence pair (right panel) ────────────────────────────────────────────
N_div  = 80
φ0     = 0.8
Δφ     = 1e-5
p_a    = Particle(0.5, 0.5, φ0)
p_b    = Particle(0.5, 0.5, φ0 + Δφ)

_, poss_a, _ = evolve!(p_a, bd, N_div)
_, poss_b, _ = evolve!(p_b, bd, N_div)

pts_a = dense_trajectory(poss_a; pts_per_seg = 10)
pts_b = dense_trajectory(poss_b; pts_per_seg = 10)
n_anim = length(pts_a)   # == length(pts_b) since same N_div, pts_per_seg

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1200, 550))

ax_erg = Axis(fig[1, 1];
              title  = "Ergodic trajectory  ($N_long bounces)",
              xlabel = "x", ylabel = "y",
              aspect = DataAspect())

ax_div = Axis(fig[1, 2];
              title  = "Sensitive dependence on initial conditions  (Δφ₀ = 10⁻⁵)",
              xlabel = "x", ylabel = "y",
              aspect = DataAspect())

for ax in (ax_erg, ax_div)
    lines!(ax, bdy; color = :black, linewidth = 2)
end

# Left: draw the full ergodic trajectory as a static coloured path.
# Colour by segment index so the trajectory "ages" visually.
n_seg = N_long * 8
xs = [pt[1] for pt in pts_erg]
ys = [pt[2] for pt in pts_erg]
lines!(ax_erg, xs, ys;
       color     = 1:length(pts_erg),
       colormap  = :plasma,
       linewidth = 0.6,
       alpha     = 0.8)

# Right: animate the two trajectories diverging.
trail_a = Observable(Point2f[])
trail_b = Observable(Point2f[])
ball_a  = Observable(pts_a[1:1])
ball_b  = Observable(pts_b[1:1])

lines!(ax_div, trail_a; color = (:dodgerblue, 0.7), linewidth = 1.5, label = "orbit A")
lines!(ax_div, trail_b; color = (:crimson,    0.7), linewidth = 1.5, label = "orbit B  (Δφ₀=10⁻⁵)")
scatter!(ax_div, ball_a; color = :dodgerblue, markersize = 10)
scatter!(ax_div, ball_b; color = :crimson,    markersize = 10)
axislegend(ax_div; position = :rt, framevisible = true)

outpath = joinpath(@__DIR__, "../../output/videos/07_billiards_stadium.mp4")
mkpath(dirname(outpath))

record(fig, outpath, 1:n_anim; framerate = 30) do k
    push!(trail_a[], pts_a[k]);  notify(trail_a)
    push!(trail_b[], pts_b[k]);  notify(trail_b)
    ball_a[] = [pts_a[k]]
    ball_b[] = [pts_b[k]]
end

@info "wrote $outpath"
