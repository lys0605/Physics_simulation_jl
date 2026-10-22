using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie
include(joinpath(@__DIR__, "../../src/em/em_utils.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Fresnel equations: reflection and transmission at a planar interface ──────
#
# Panels:
#   (a) R_s, R_p vs θ_i for glass→air (n1=1.5, n2=1.0): shows TIR at θ_c
#   (b) R_s, R_p vs θ_i for air→glass (n1=1.0, n2=1.5): shows Brewster angle
#   (c) Snell's law θ_t vs θ_i for three n2 values
#   (d) Phase shift δ_s, δ_p vs θ_i for air→glass
#   (e) Wavefront diagram at θ_i = 35° (air→glass)

# ── Common angle grid ─────────────────────────────────────────────────────────

θs = range(0, π/2; length=500)

# ── Media sets ───────────────────────────────────────────────────────────────

n_air   = 1.0
n_glass = 1.5
n_water = 1.33
n_dense = 2.0

θ_B  = brewster_angle(n_air, n_glass)
θ_c  = critical_angle(n_glass, n_air)

@info "Brewster angle (air→glass): $(round(rad2deg(θ_B), digits=2))°"
@info "Critical angle (glass→air): $(round(rad2deg(θ_c), digits=2))°"

# ── Compute Fresnel coefficients ──────────────────────────────────────────────

function sweep(n1, n2, θs)
    Rs = [fresnel(n1, n2, θ).Rs for θ in θs]
    Rp = [fresnel(n1, n2, θ).Rp for θ in θs]
    Ts = [fresnel(n1, n2, θ).Ts for θ in θs]
    Tp = [fresnel(n1, n2, θ).Tp for θ in θs]
    (; Rs, Rp, Ts, Tp)
end

sw_ag = sweep(n_air,   n_glass, θs)
sw_ga = sweep(n_glass, n_air,   θs)
θ_deg = rad2deg.(θs)

# Snell's law: θ_t vs θ_i for different n2
snell_t(n1, n2, θi) = begin
    st = (n1/n2)*sin(θi)
    st >= 1.0 ? NaN : rad2deg(asin(st))
end
θt_glass = [snell_t(n_air, n_glass, θ) for θ in θs]
θt_water = [snell_t(n_air, n_water, θ) for θ in θs]
θt_dense = [snell_t(n_air, n_dense, θ) for θ in θs]

# Phase shifts: arg(r) = 0 or π
rs_phase(n1, n2, θ) = fresnel(n1, n2, θ).rs < 0 ? 180.0 : 0.0
rp_phase(n1, n2, θ) = fresnel(n1, n2, θ).rp < 0 ? 180.0 : 0.0
δs_ag = [rs_phase(n_air, n_glass, θ) for θ in θs]
δp_ag = [rp_phase(n_air, n_glass, θ) for θ in θs]

# ── Figure ────────────────────────────────────────────────────────────────────

set_theme!(theme_light())
fig = Figure(size = (1300, 920))

Label(fig[0, :], "Fresnel Equations — Reflection and Transmission at a Planar Interface";
      fontsize = 15, font = :bold)

# ── (a) Glass → air: TIR ─────────────────────────────────────────────────────
ax_a = Axis(fig[1, 1];
            title   = "(a)  glass → air  (n₁=1.5, n₂=1.0)",
            xlabel  = "θᵢ (°)",
            ylabel  = "Reflectance R",
            limits  = (0, 90, -0.05, 1.10),
            titlesize = 12)
vlines!(ax_a, [rad2deg(θ_c)]; color = (:gray, 0.5), linewidth = 1.5, linestyle = :dash)
text!(ax_a, rad2deg(θ_c) + 0.5, 0.55; text = "θ_c", color = :gray, fontsize = 10)
lines!(ax_a, θ_deg, sw_ga.Rs; color = :dodgerblue, linewidth = 2.5, label = "Rₛ")
lines!(ax_a, θ_deg, sw_ga.Rp; color = :crimson,    linewidth = 2.5, label = "Rₚ")
hlines!(ax_a, [1.0]; color = (:black, 0.3), linewidth = 1, linestyle = :dot)
axislegend(ax_a; position = :lt, labelsize = 11)

# ── (b) Air → glass: Brewster angle ──────────────────────────────────────────
ax_b = Axis(fig[1, 2];
            title   = "(b)  air → glass  (n₁=1.0, n₂=1.5)",
            xlabel  = "θᵢ (°)",
            ylabel  = "Reflectance R",
            limits  = (0, 90, -0.05, 1.10),
            titlesize = 12)
vlines!(ax_b, [rad2deg(θ_B)]; color = (:gray, 0.5), linewidth = 1.5, linestyle = :dash)
text!(ax_b, rad2deg(θ_B) + 0.5, 0.55; text = "θ_B", color = :gray, fontsize = 10)
lines!(ax_b, θ_deg, sw_ag.Rs; color = :dodgerblue, linewidth = 2.5, label = "Rₛ")
lines!(ax_b, θ_deg, sw_ag.Rp; color = :crimson,    linewidth = 2.5, label = "Rₚ")
axislegend(ax_b; position = :lt, labelsize = 11)

# ── (c) Snell's law ───────────────────────────────────────────────────────────
ax_c = Axis(fig[1, 3];
            title   = "(c)  Snell's law  θₜ vs θᵢ",
            xlabel  = "θᵢ (°)",
            ylabel  = "θₜ (°)",
            limits  = (0, 90, 0, 90),
            titlesize = 12)
lines!(ax_c, θ_deg, θ_deg;      color = :black,    linewidth = 1.0, linestyle = :dash,
       label = "θᵢ (vacuum)")
lines!(ax_c, θ_deg, θt_water;   color = :dodgerblue, linewidth = 2, label = "n₂=1.33 (water)")
lines!(ax_c, θ_deg, θt_glass;   color = :crimson,    linewidth = 2, label = "n₂=1.50 (glass)")
lines!(ax_c, θ_deg, θt_dense;   color = :darkorange, linewidth = 2, label = "n₂=2.00 (dense)")
axislegend(ax_c; position = :lt, labelsize = 10)

# ── (d) Energy conservation check  (R+T=1) ───────────────────────────────────
ax_d = Axis(fig[2, 1];
            title   = "(d)  Energy conservation  R + T = 1",
            xlabel  = "θᵢ (°)",
            ylabel  = "R + T",
            limits  = (0, 90, 0.90, 1.10),
            titlesize = 12)
lines!(ax_d, θ_deg, sw_ag.Rs .+ sw_ag.Ts; color = :dodgerblue, linewidth = 2, label = "s-pol")
lines!(ax_d, θ_deg, sw_ag.Rp .+ sw_ag.Tp; color = :crimson,    linewidth = 2, label = "p-pol")
hlines!(ax_d, [1.0]; color = :black, linewidth = 1, linestyle = :dash)
axislegend(ax_d; position = :lb, labelsize = 11)

# ── (e) Phase shifts ─────────────────────────────────────────────────────────
ax_e = Axis(fig[2, 2];
            title   = "(e)  Phase shift on reflection  (air → glass)",
            xlabel  = "θᵢ (°)",
            ylabel  = "δ (°)",
            limits  = (0, 90, -10, 200),
            yticks  = ([0, 90, 180], ["0°", "90°", "180°"]),
            titlesize = 12)
vlines!(ax_e, [rad2deg(θ_B)]; color = (:gray, 0.5), linewidth = 1.5, linestyle = :dash)
lines!(ax_e, θ_deg, δs_ag; color = :dodgerblue, linewidth = 2.5, label = "δₛ")
lines!(ax_e, θ_deg, δp_ag; color = :crimson,    linewidth = 2.5, label = "δₚ")
text!(ax_e, rad2deg(θ_B) + 0.5, 95; text = "θ_B", color = :gray, fontsize = 10)
axislegend(ax_e; position = :lb, labelsize = 11)

# ── (f) Wavefront diagram ─────────────────────────────────────────────────────
ax_f = Axis(fig[2, 3];
            title   = "(f)  Wavefronts at θᵢ = 35°  (air → glass)",
            aspect  = DataAspect(),
            limits  = (-2.0, 2.0, -1.8, 1.8),
            xticksvisible = false, yticksvisible = false,
            xticklabelsvisible = false, yticklabelsvisible = false,
            titlesize = 12)

# Interface line
hlines!(ax_f, [0.0]; color = :black, linewidth = 2)
text!(ax_f, 1.6, 0.07; text = "interface", fontsize = 9, color = :black)

# Media labels
text!(ax_f, -1.8, 1.5; text = "n₁ = 1.0 (air)",   fontsize = 10)
text!(ax_f, -1.8, -1.5; text = "n₂ = 1.5 (glass)", fontsize = 10)

θi_wf = deg2rad(35.0)
f_wf  = fresnel(n_air, n_glass, θi_wf)
θt_wf = f_wf.θt

# Draw incident, reflected, and transmitted wavefronts (as lines perpendicular to k)
function wavefronts!(ax, kx, kz, x0, z0, n_lines, spacing, color, label)
    klen = sqrt(kx^2 + kz^2)
    kx_n = kx / klen;  kz_n = kz / klen
    # perp direction
    px = -kz_n;  pz = kx_n
    for i in 0:(n_lines-1)
        cx = x0 + i * spacing * kx_n
        cz = z0 + i * spacing * kz_n
        len = 0.5
        x1 = cx - len*px;  z1 = cz - len*pz
        x2 = cx + len*px;  z2 = cz + len*pz
        lines!(ax, [x1, x2], [z1, z2]; color, linewidth = 1.5,
               label = i == 0 ? label : "")
    end
    # Arrow for direction
    ax_x = x0 + 3 * spacing * kx_n
    ax_z = z0 + 3 * spacing * kz_n
    arrows!(ax, [ax_x - 0.15*kx_n], [ax_z - 0.15*kz_n],
            [0.3*kx_n], [0.3*kz_n]; arrowsize = 10, color)
end

λ = 0.45
n_wf = 5
# Incident: comes from upper-left, k = (sin θi, -cos θi) in (x, z)
wavefronts!(ax_f,  sin(θi_wf), -cos(θi_wf), -1.2, 1.4, n_wf, λ, :dodgerblue, "incident")
# Reflected: (sin θi, +cos θi)
wavefronts!(ax_f,  sin(θi_wf),  cos(θi_wf),  0.0, 0.0, n_wf, λ*0.8, :crimson,   "reflected")
# Transmitted: n1/n2 wavelength compression
λt = λ * n_air / n_glass
wavefronts!(ax_f,  sin(θt_wf), -cos(θt_wf),  0.0, 0.0, n_wf, λt, :seagreen,  "transmitted")

axislegend(ax_f; position = :rb, labelsize = 9)

# ── Save ──────────────────────────────────────────────────────────────────────

mkpath(joinpath(@__DIR__, "../../output/plots"))
outpath = joinpath(@__DIR__, "../../output/plots/14_fresnel.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
