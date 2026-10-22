using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie
include(joinpath(@__DIR__, "../../src/em/em_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Radiation from an oscillating Hertzian dipole  p(t) = p₀ cos(ωt) ẑ ───────
#
# The animation shows a snapshot of the retarded E_z field in the xz-plane as
# the dipole oscillates, revealing the characteristic sin²θ wavefront pattern
# and the transition from near-field (1/r²) to radiation (1/r) behavior.
#
# Side panel: polar power pattern  dP/dΩ ∝ sin²θ
#
# Convention: c = 1, k = ω = 2π (wavelength λ = 1 in these units)

const k_rad  = 2π          # wave-number = ω (c=1)
const N_grid = 300         # spatial resolution
const L      = 3.5         # domain half-size in wavelengths
const N_frames = 90        # frames per period
const framerate = 30

xs_r = range(-L, L; length=N_grid)
zs_r = range(-L, L; length=N_grid)

# Pre-compute Ez snapshot for a given phase φ = ωt
function ez_snapshot(φ)
    [hertzian_fields(x, z, φ / k_rad; k=k_rad, p0=1.0, reg=0.15)[2]
     for x in xs_r, z in zs_r]
end

# ── Polar radiation pattern (static) ─────────────────────────────────────────

θ_pol = range(0, 2π; length=500)
r_pol = sin.(θ_pol).^2

# ── Figure setup ──────────────────────────────────────────────────────────────

set_theme!(theme_dark())
fig = Figure(size = (1100, 720))

Label(fig[0, :], "Hertzian Dipole Radiation  p(t) = cos(ωt) ẑ  (c = 1, λ = 1)";
      fontsize = 14, font = :bold)

# Main 2D field panel
ax_field = Axis(fig[1:2, 1];
                aspect  = DataAspect(),
                limits  = (-L, L, -L, L),
                title   = "E_z(x, z, t)  in the xz-plane",
                xlabel  = "x / λ",
                ylabel  = "z / λ",
                titlesize = 12)

# Polar pattern panel
ax_polar = Axis(fig[1, 2];
                title   = "Power pattern  dP/dΩ ∝ sin²θ",
                aspect  = DataAspect(),
                limits  = (-1.25, 1.25, -1.25, 1.25),
                xticksvisible=false, yticksvisible=false,
                xticklabelsvisible=false, yticklabelsvisible=false,
                titlesize = 12)

# Info panel
ax_info = Axis(fig[2, 2];
               limits = (0, 1, 0, 1),
               leftspinevisible=false, rightspinevisible=false,
               topspinevisible=false, bottomspinevisible=false,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)

colsize!(fig.layout, 2, Fixed(250))

# ── Static polar pattern ──────────────────────────────────────────────────────

# Fill the pattern
θ_fill = range(0, 2π; length=500)
rfill  = sin.(θ_fill).^2
poly!(ax_polar, rfill .* cos.(θ_fill), rfill .* sin.(θ_fill);
      color=(:darkorange, 0.35))
lines!(ax_polar, rfill .* cos.(θ_fill), rfill .* sin.(θ_fill);
       color=:darkorange, linewidth=2.5)

# Dipole axis indicator
arrows!(ax_polar, [0.0], [-1.1], [0.0], [2.2];
        arrowsize=10, color=:white, linewidth=1.5)
text!(ax_polar, 0.06, 1.05; text="z (dipole)", color=:white, fontsize=9)

# Angle annotations
for θ_ann in [0, 30, 60, 90, 120, 150]
    θ_r = deg2rad(Float64(θ_ann))
    r_ann = sin(θ_r)^2
    tx = r_ann * sin(θ_r) + 0.05;  tz = r_ann * cos(θ_r)
    text!(ax_polar, tx, tz; text="$(θ_ann)°", color=(:white, 0.6), fontsize=7)
end

# ── Info text ─────────────────────────────────────────────────────────────────

text!(ax_info, 0.05, 0.90;
      text = "Near field:  E ∝ 1/r³\nInduction:  E ∝ 1/r²\nRadiation: E ∝ 1/r\n\n" *
             "kλ = 2π  (λ = 1)\nPower: P ∝ (kd)² (ωp₀)²\n\n" *
             "Color: E_z  (red = +, blue = −)",
      color=:white, fontsize=10, align=(:left, :top))

# ── Observables ───────────────────────────────────────────────────────────────

phase_obs   = Observable(0.0)
Ezmat_obs   = @lift ez_snapshot($phase_obs)
phase_text  = @lift "ωt = $(round($phase_obs / π, digits=2))π"

# Clamp for color scale
Ez_clipped  = @lift clamp.($Ezmat_obs, -0.8, 0.8)

# ── Plot objects ──────────────────────────────────────────────────────────────

hm = heatmap!(ax_field, xs_r, zs_r, Ez_clipped;
              colormap=:RdBu_r, colorrange=(-0.8, 0.8), interpolate=true)
Colorbar(fig[1:2, 3], hm; label="E_z", labelsize=10)

# Dipole marker
scatter!(ax_field, [0.0], [0.0]; color=:yellow, markersize=14, marker=:star5)

# Ring at r = λ to show near/far boundary
θ_ring = range(0, 2π; length=300)
lines!(ax_field, cos.(θ_ring), sin.(θ_ring);
       color=(:white, 0.25), linewidth=1, linestyle=:dash)
text!(ax_field, 0.72, 0.72; text="r=λ", color=(:white, 0.5), fontsize=8)

# Phase text
text!(ax_field, -L + 0.1, L - 0.2; text=phase_text,
      color=:white, fontsize=11, align=(:left, :top))

# ── Record ────────────────────────────────────────────────────────────────────

phases = range(0, 2π; length=N_frames+1)[1:N_frames]

mkpath(joinpath(@__DIR__, "../../output/videos"))
outpath = joinpath(@__DIR__, "../../output/videos/17_radiation.mp4")

record(fig, outpath, phases; framerate) do φ
    phase_obs[] = φ
end

@info "wrote $outpath"
