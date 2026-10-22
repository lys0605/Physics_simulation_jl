# Thin helpers on top of DynamicalBilliards.jl for the billiards tutorial scripts.
# The heavy lifting (collision geometry, Lyapunov spectrum, boundary maps) is all
# done by DynamicalBilliards.jl — this file only provides:
#   • make_circle(R; n_sides)          — polygon approximation of a circle billiard
#   • boundary_outline(bd)             — closed Vector{Point2f} for plotting
#   • dense_trajectory(poss; pps)      — interpolate collision sequence for animation
#
# Every script should `include` this file (not `using`) so it picks up DynamicalBilliards
# and CairoMakie implicitly; those packages must already be loaded in the script.

using DynamicalBilliards

# ─────────────────────────── billiard constructors ────────────────────────────

"""
    make_circle(R=1.0; n_sides=100) → Billiard

Regular `n_sides`-gon inscribed in a circle of radius `R`.  For n ≥ 100 this is
an excellent approximation to the integrable circular billiard: angular momentum
is conserved to O(1/n²) and the Poincaré section shows near-horizontal invariant
curves instead of the scattered cloud produced by chaotic billiards.
"""
make_circle(R::Float64 = 1.0; n_sides::Int = 100) = billiard_polygon(n_sides, R)

# ─────────────────────────── boundary drawing ─────────────────────────────────

function _obs_outline(w::Union{InfiniteWall, FiniteWall, PeriodicWall})
    [Point2f(w.sp), Point2f(w.ep)]
end

function _obs_outline(s::Semicircle)
    # facedir points INTO the billiard, so the arc faces −facedir
    θ0 = atan(-s.facedir[2], -s.facedir[1])
    θs = range(θ0 - π/2, θ0 + π/2; length=60)
    [Point2f(s.c[1] + s.r*cos(θ), s.c[2] + s.r*sin(θ)) for θ in θs]
end

function _obs_outline(d::Union{Disk, Antidot, RandomDisk})
    θs = range(0, 2π; length=120)
    [Point2f(d.c[1] + d.r*cos(θ), d.c[2] + d.r*sin(θ)) for θ in θs]
end

function _obs_outline(e::Ellipse)
    θs = range(0, 2π; length=200)
    [Point2f(e.c[1] + e.a*cos(θ), e.c[2] + e.b*sin(θ)) for θ in θs]
end

"""
    boundary_outline(bd) → Vector{Point2f}

Closed polygon tracing all obstacles in `bd` in order.  For billiards where the
outer boundary is stored first (stadium, rectangle, polygon) this traces the
perimeter.  For the Sinai billiard (Disk first, then rectangle walls) use
`outer_boundary` / `inner_obstacles` instead.
"""
function boundary_outline(bd::Billiard)
    pts = Point2f[]
    for obs in bd
        append!(pts, _obs_outline(obs))
    end
    isempty(pts) || push!(pts, pts[1])
    pts
end

"""
    outer_boundary(bd), inner_obstacles(bd)

For billiards with inner scatterers (Sinai-type): split the boundary into
the outer perimeter and a vector of (center, radius) tuples for circular
inner obstacles.
"""
function outer_boundary(bd::Billiard)
    pts = Point2f[]
    for obs in bd
        obs isa Union{Disk, Antidot, RandomDisk} && continue
        append!(pts, _obs_outline(obs))
    end
    isempty(pts) || push!(pts, pts[1])
    pts
end

function inner_obstacles(bd::Billiard)
    [(obs.c, obs.r) for obs in bd if obs isa Union{Disk, Antidot, RandomDisk}]
end

# ─────────────────────────── trajectory helpers ───────────────────────────────

"""
    dense_trajectory(poss; pts_per_seg=20) → Vector{Point2f}

Expand a collision sequence `poss` (as returned by `evolve!`, including the
initial position at index 1) into a dense path by linearly interpolating each
straight-line segment.  The particle travels in a straight line between
collisions, so this is exact — not an approximation.
"""
function dense_trajectory(poss; pts_per_seg::Int = 20)
    pts = Point2f[]
    sizehint!(pts, length(poss) * pts_per_seg)
    for i in 1:length(poss)-1
        p0 = Point2f(poss[i]);  p1 = Point2f(poss[i+1])
        for k in 0:pts_per_seg-1
            α = k / pts_per_seg
            push!(pts, Point2f((1-α)*p0[1]+α*p1[1], (1-α)*p0[2]+α*p1[2]))
        end
    end
    push!(pts, Point2f(poss[end]))
    pts
end
