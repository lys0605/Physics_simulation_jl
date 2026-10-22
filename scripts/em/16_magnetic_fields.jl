using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie
include(joinpath(@__DIR__, "../../src/em/em_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Magnetic field configurations ─────────────────────────────────────────────
#
# Panels (all in the xz-plane):
#   (a) Magnetic dipole (analytical m = ẑ)
#   (b) Single circular current loop (Biot-Savart, R=1, I=1)
#   (c) Helmholtz coil: two loops at z = ±R/2 (uniform field region visible)
#   (d) Anti-Helmholtz (Maxwell coil): loops at z = ±R/2 with opposite currents
#       → linear gradient, used in magnetic traps

# ── Grid ──────────────────────────────────────────────────────────────────────

N  = 60       # streamplot grid (too dense → slow)
L  = 2.5
xs_g = range(-L, L; length=N)
zs_g = range(-L, L; length=N)

# ── Magnitude helper ─────────────────────────────────────────────────────────

clip_field(Bx, Bz, Bmax=5.0) = Point2f(clamp(Bx, -Bmax, Bmax),
                                         clamp(Bz, -Bmax, Bmax))

# ── Figure ────────────────────────────────────────────────────────────────────

set_theme!(theme_dark())
fig = Figure(size = (1200, 1100))
Label(fig[0, :], "Magnetic Field Configurations (xz-plane)";
      fontsize = 15, font = :bold)

# ── (a) Magnetic dipole ───────────────────────────────────────────────────────
ax_a = Axis(fig[1, 1];
            title  = "(a)  Magnetic dipole  m = ẑ",
            aspect = DataAspect(), limits = (-L, L, -L, L),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false,
            titlesize=13)

B_dip_mag = [let (Bx, Bz) = bfield_dipole_xz(1.0, x, z)
                 sqrt(Bx^2 + Bz^2) end
             for x in xs_g, z in zs_g]
heatmap!(ax_a, xs_g, zs_g, log1p.(clamp.(B_dip_mag, 0, 20)');
         colormap=:magma, interpolate=true)
streamplot!(ax_a,
            (p, _) -> begin
                Bx, Bz = bfield_dipole_xz(1.0, p[1], p[2])
                Point2f(Bx, Bz)
            end,
            xs_g, zs_g;
            colormap=[:white], linewidth=1.0, density=1.2, arrow_size=8)
scatter!(ax_a, [0.0], [0.0]; color=:yellow, markersize=10, marker=:star5)
text!(ax_a, 0.1, 0.15; text="m", color=:yellow, fontsize=11)

# ── (b) Single circular loop ──────────────────────────────────────────────────
R_loop = 1.0
ax_b = Axis(fig[1, 2];
            title  = "(b)  Circular loop  R=1, I=1",
            aspect = DataAspect(), limits = (-L, L, -L, L),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false,
            titlesize=13)

# For xz-plane, evaluate at (x, 0, z) → (Bx, Bz)
B_loop_mag = [let (Bx, Bz) = bfield_loop_xz_v2(R_loop, 1.0, x, z)
                  sqrt(Bx^2 + Bz^2) end
              for x in xs_g, z in zs_g]
heatmap!(ax_b, xs_g, zs_g, log1p.(clamp.(B_loop_mag, 0, 20)');
         colormap=:magma, interpolate=true)
streamplot!(ax_b,
            (p, _) -> begin
                Bx, Bz = bfield_loop_xz_v2(R_loop, 1.0, p[1], p[2])
                Point2f(Bx, Bz)
            end,
            xs_g, zs_g;
            colormap=[:white], linewidth=1.0, density=1.2, arrow_size=8)
# Show loop cross-section
scatter!(ax_b, [-R_loop, R_loop], [0.0, 0.0]; color=:cyan, markersize=12,
         marker=:circle, strokewidth=1.5, strokecolor=:white)
text!(ax_b, -R_loop-0.2, -0.3; text="⊗", color=:cyan, fontsize=14)
text!(ax_b,  R_loop+0.05, -0.3; text="⊙", color=:cyan, fontsize=14)

# ── (c) Helmholtz coil ────────────────────────────────────────────────────────
R_hh = 1.0
z_hh = R_hh / 2
ax_c = Axis(fig[2, 1];
            title  = "(c)  Helmholtz coil  (same current, z = ±R/2)",
            aspect = DataAspect(), limits = (-L, L, -L, L),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false,
            titlesize=13)

function b_helmholtz(x, z)
    Bx1, Bz1 = bfield_loop_xz_v2(R_hh, 1.0, x, z - z_hh)
    Bx2, Bz2 = bfield_loop_xz_v2(R_hh, 1.0, x, z + z_hh)
    Bx1 + Bx2, Bz1 + Bz2
end
B_hh_mag = [let (Bx, Bz) = b_helmholtz(x, z); sqrt(Bx^2+Bz^2) end
            for x in xs_g, z in zs_g]
heatmap!(ax_c, xs_g, zs_g, log1p.(clamp.(B_hh_mag, 0, 20)');
         colormap=:magma, interpolate=true)
streamplot!(ax_c, (p, _) -> Point2f(b_helmholtz(p[1], p[2])...),
            xs_g, zs_g; colormap=[:white], linewidth=1.0, density=1.2, arrow_size=8)
for z0 in [-z_hh, z_hh]
    scatter!(ax_c, [-R_hh, R_hh], [z0, z0]; color=:cyan, markersize=12,
             strokewidth=1.5, strokecolor=:white)
end
# Annotate uniform-field region
poly!(ax_c, Rect(-0.4, -0.4, 0.8, 0.8); color=(:yellow, 0.10), strokecolor=:yellow,
      strokewidth=1.5, linestyle=:dash)
text!(ax_c, 0.05, 0.5; text="uniform", color=:yellow, fontsize=9)

# ── (d) Anti-Helmholtz (Maxwell) coil ────────────────────────────────────────
ax_d = Axis(fig[2, 2];
            title  = "(d)  Anti-Helmholtz (opposite currents, z = ±R/2)",
            aspect = DataAspect(), limits = (-L, L, -L, L),
            xticksvisible=false, yticksvisible=false,
            xticklabelsvisible=false, yticklabelsvisible=false,
            titlesize=13)

function b_antihh(x, z)
    Bx1, Bz1 = bfield_loop_xz_v2(R_hh,  1.0, x, z - z_hh)
    Bx2, Bz2 = bfield_loop_xz_v2(R_hh, -1.0, x, z + z_hh)
    Bx1 + Bx2, Bz1 + Bz2
end
B_ahh_mag = [let (Bx, Bz) = b_antihh(x, z); sqrt(Bx^2+Bz^2) end
             for x in xs_g, z in zs_g]
heatmap!(ax_d, xs_g, zs_g, log1p.(clamp.(B_ahh_mag, 0, 20)');
         colormap=:magma, interpolate=true)
streamplot!(ax_d, (p, _) -> Point2f(b_antihh(p[1], p[2])...),
            xs_g, zs_g; colormap=[:white], linewidth=1.0, density=1.2, arrow_size=8)
for (z0, sgn) in [(-z_hh, "+I"), (z_hh, "−I")]
    scatter!(ax_d, [-R_hh, R_hh], [z0, z0]; color=:cyan, markersize=12,
             strokewidth=1.5, strokecolor=:white)
    text!(ax_d, R_hh+0.15, z0+0.05; text=sgn, color=:cyan, fontsize=9)
end
scatter!(ax_d, [0.0], [0.0]; color=:red, markersize=10, marker=:xcross)
text!(ax_d, 0.1, 0.15; text="B=0", color=:red, fontsize=9)

# ── Save ──────────────────────────────────────────────────────────────────────

mkpath(joinpath(@__DIR__, "../../output/plots"))
outpath = joinpath(@__DIR__, "../../output/plots/16_magnetic_fields.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
