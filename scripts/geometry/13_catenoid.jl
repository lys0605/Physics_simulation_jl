using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie

# ── Catenoid: the minimal surface of revolution ────────────────────────────────
#
# Rotating the catenary  r(z) = a cosh(z/a)  about the z-axis gives the
# catenoid — the unique (up to scaling) minimal surface of revolution with
# mean curvature H = 0 everywhere.
#
# For two coaxial rings of equal radius R at z = ±h, the neck parameter a
# must satisfy the boundary condition
#
#   R = a cosh(h/a)                                             (*)
#
# which is a transcendental equation.  For a given R the function
# a cosh(h/a) is non-monotone in a and has a minimum value
#
#   R_min = h cosh(u_c) / u_c      where  u_c tanh(u_c) = 1,  u_c ≈ 1.1997
#
# When h/R > (h/R)_c ≡ u_c / cosh(u_c) ≈ 0.6627, equation (*) has no real
# solution: the soap film cannot form a catenoid and snaps into two flat discs.
#
# The animation stretches the rings from h ≈ 0 toward and past this threshold,
# showing the narrowing waist and the final snap.

# ── Critical-ratio solver ─────────────────────────────────────────────────────

# u_c satisfies u_c tanh(u_c) = 1  (Newton's method)
const u_c = let u = 1.2
    for _ in 1:60
        fu  = u * tanh(u) - 1
        dfu = tanh(u) + u * (1 - tanh(u)^2)
        u  -= fu / dfu
    end
    u
end
const CRIT = u_c / cosh(u_c)      # critical h/R ≈ 0.6627

"""
    solve_catenoid(R, h) → a

Neck radius a (outer / stable solution) for rings of radius R at ±h.
Returns NaN when h/R ≥ CRIT (no catenoid exists).
Bisects in (h/u_c, R) where the root is bracketed.
"""
function solve_catenoid(R::Real, h::Real; tol::Real = 1e-13)
    h < 1e-14 && return float(R)
    h / R >= CRIT && return NaN
    aL, aR = h / u_c, R * (1 - 1e-12)
    for _ in 1:120
        am = (aL + aR) / 2
        (am * cosh(h / am) - R < 0) ? (aL = am) : (aR = am)
        (aR - aL) < tol && break
    end
    (aL + aR) / 2
end

# ── Mesh builders ─────────────────────────────────────────────────────────────

const Nθ    = 80
const Nz    = 80
const θgrid = range(0, 2π; length = Nθ + 1)   # +1 closes the seam

# Catenoid surface: (θ, z) → (a cosh(z/a) cosθ,  a cosh(z/a) sinθ,  z)
function cat_mesh(a::Real, h::Real)
    zg = range(-h, h; length = Nz)
    xs = [a * cosh(z / a) * cos(θ) for θ in θgrid, z in zg]
    ys = [a * cosh(z / a) * sin(θ) for θ in θgrid, z in zg]
    Zs = [Float64(z)               for θ in θgrid, z in zg]
    xs, ys, Zs
end

# Flat disc at height z0
function disc_mesh(R::Real, z0::Real; Nr::Int = 30)
    rg = range(0, R; length = Nr)
    xs = [r * cos(θ) for θ in θgrid, r in rg]
    ys = [r * sin(θ) for θ in θgrid, r in rg]
    Zs = fill(Float64(z0), Nθ + 1, Nr)
    xs, ys, Zs
end

# ── Animation schedule ────────────────────────────────────────────────────────

const R = 1.0
const h_c = CRIT * R          # ≈ 0.6627

# Approach: 160 frames stretching toward snap; 40 frames settling after snap
h_pre  = Vector(range(0.06, h_c * 0.994;  length = 160))
h_post = Vector(range(h_c * 1.006, h_c * 1.30; length = 40))
h_all  = vcat(h_pre, h_post)
N_pre  = length(h_pre)

@info "Catenoid: CRIT = $(round(CRIT, digits=5))  u_c = $(round(u_c, digits=5))"
@info "Frames: $(length(h_all))  (pre-snap: $N_pre, post-snap: $(length(h_post)))"

# ── Figure setup ──────────────────────────────────────────────────────────────

set_theme!(theme_dark())
fig = Figure(size = (1100, 780))

# 3-D main view
ax3 = Axis3(fig[1:3, 1];
            aspect           = :data,
            limits           = (-1.7, 1.7, -1.7, 1.7, -1.0, 1.0),
            elevation        = 0.26π,
            azimuth          = 0.38π,
            protrusions      = (10, 10, 10, 10),
            xlabelvisible    = false,
            ylabelvisible    = false,
            zlabelvisible    = false,
            xticksvisible    = false,
            yticksvisible    = false,
            zticksvisible    = false,
            xticklabelsvisible = false,
            yticklabelsvisible = false,
            zticklabelsvisible = false)

# Side panel: a/R vs h/R
ax_neck = Axis(fig[1, 2];
               title    = "neck radius  a / R",
               xlabel   = "separation  h / R",
               ylabel   = "a / R",
               limits   = (0.0, CRIT * 1.12, 0.0, 1.05),
               titlesize = 11, xlabelsize = 10, ylabelsize = 10)

let hc = range(0.001, h_c * 0.9999; length = 300)
    lines!(ax_neck, hc ./ R, [solve_catenoid(R, h) / R for h in hc];
           color = :cyan, linewidth = 1.5)
end
vlines!(ax_neck, [CRIT]; color = :tomato, linewidth = 1.5, linestyle = :dash)
text!(ax_neck, CRIT * 1.02, 0.15;
      text = "snap\nh/R=$(round(CRIT,digits=4))",
      color = :tomato, fontsize = 8, align = (:left, :bottom))

# Text info box
ax_info = Axis(fig[2:3, 2];
               limits = (0, 1, 0, 1),
               leftspinevisible   = false, rightspinevisible  = false,
               topspinevisible    = false, bottomspinevisible = false,
               xticksvisible      = false, yticksvisible      = false,
               xticklabelsvisible = false, yticklabelsvisible = false)

Label(fig[0, :], "Catenoid  —  Soap Film Stretching";
      fontsize = 15, font = :bold)
colsize!(fig.layout, 2, Fixed(230))

# ── Observables ───────────────────────────────────────────────────────────────

h_obs = Observable(h_all[1])
a_obs = @lift solve_catenoid(R, $h_obs)

# Catenoid mesh (zeros when snapped — invisible anyway)
cat_data = @lift begin
    h, a = $h_obs, $a_obs
    isnan(a) ? (zeros(Nθ+1, Nz), zeros(Nθ+1, Nz), zeros(Nθ+1, Nz)) :
               cat_mesh(a, h)
end
cat_xs = @lift $cat_data[1]
cat_ys = @lift $cat_data[2]
cat_zs = @lift $cat_data[3]

# Post-snap flat discs
top_data = @lift disc_mesh(R,  $h_obs)
bot_data = @lift disc_mesh(R, -$h_obs)
top_xs = @lift $top_data[1]; top_ys = @lift $top_data[2]; top_zs = @lift $top_data[3]
bot_xs = @lift $bot_data[1]; bot_ys = @lift $bot_data[2]; bot_zs = @lift $bot_data[3]

# End-ring z coordinates (x, y are constant)
ring_cos = R .* cos.(θgrid)
ring_sin = R .* sin.(θgrid)
ring_z_top = @lift fill($h_obs,  Nθ + 1)
ring_z_bot = @lift fill(-$h_obs, Nθ + 1)

# Info text
info_str = @lift begin
    h, a = $h_obs, $a_obs
    hr = round(h / R, digits = 3)
    if isnan(a)
        "h/R = $hr\n\n⚡  SNAPPED\n\n(two flat discs)"
    else
        "h / R  =  $hr\n" *
        "h_c/ R =  $(round(CRIT, digits=4))\n\n" *
        "neck a/R =  $(round(a/R, digits=4))\n" *
        "waist    =  $(round(a, digits=4)) R"
    end
end

# Progress dot on the curve
dot_h = @lift [$h_obs / R]
dot_a = @lift [isnan($a_obs) ? 0.0 : $a_obs / R]

# ── Plot objects ──────────────────────────────────────────────────────────────

# Catenoid surface (plasma = purple→yellow mapped to z)
p_cat = surface!(ax3, cat_xs, cat_ys, cat_zs;
                 colormap   = :ice,
                 colorrange = (-h_c, h_c),
                 shading    = true)

# Flat discs (post-snap)
p_top = surface!(ax3, top_xs, top_ys, top_zs;
                 color   = (:steelblue, 0.88),
                 shading = true)
p_bot = surface!(ax3, bot_xs, bot_ys, bot_zs;
                 color   = (:steelblue, 0.88),
                 shading = true)

# Initial visibility: catenoid on, discs off
p_top.visible[] = false
p_bot.visible[] = false

# Toggle on snap event
on(a_obs) do a
    snp = isnan(a)
    p_cat.visible[] = !snp
    p_top.visible[] =  snp
    p_bot.visible[] =  snp
end

# End rings (always visible, track h)
lines!(ax3, ring_cos, ring_sin, ring_z_top; color = :white, linewidth = 2.5)
lines!(ax3, ring_cos, ring_sin, ring_z_bot; color = :white, linewidth = 2.5)

# Progress dot
scatter!(ax_neck, dot_h, dot_a; color = :yellow, markersize = 10)

# Info text
text!(ax_info, 0.08, 0.88; text = info_str, color = :white,
      fontsize = 12, align = (:left, :top))

# ── Record ────────────────────────────────────────────────────────────────────

mkpath(joinpath(@__DIR__, "../../output/videos"))
outpath = joinpath(@__DIR__, "../../output/videos/13_catenoid.mp4")

record(fig, outpath, h_all; framerate = 30) do h
    h_obs[] = h
end

@info "wrote $outpath"
