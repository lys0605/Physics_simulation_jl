# Electromagnetic utilities (Gaussian-like units: 1/4πε₀ = 1, μ₀/4π = 1, c = 1)

using LinearAlgebra

# ── Fresnel equations ────────────────────────────────────────────────────────

"""
    fresnel(n1, n2, θi) → NamedTuple

Fresnel reflection/transmission coefficients for a planar interface.
Returns (rs, rp, ts, tp, Rs, Rp, Ts, Tp, θt).
TIR (total internal reflection) when θi ≥ θ_c: returns R=1, T=0.
"""
function fresnel(n1::Real, n2::Real, θi::Real)
    sinθt² = (n1 / n2)^2 * sin(θi)^2
    if sinθt² >= 1.0
        θt = π / 2
        return (; rs=-1.0, rp=1.0, ts=0.0, tp=0.0,
                  Rs=1.0, Rp=1.0, Ts=0.0, Tp=0.0, θt)
    end
    cosθi = cos(θi)
    cosθt = sqrt(1.0 - sinθt²)
    rs = (n1*cosθi - n2*cosθt) / (n1*cosθi + n2*cosθt)
    rp = (n2*cosθi - n1*cosθt) / (n2*cosθi + n1*cosθt)
    ts = 2n1*cosθi / (n1*cosθi + n2*cosθt)
    tp = 2n1*cosθi / (n2*cosθi + n1*cosθt)
    Rs = rs^2
    Rp = rp^2
    Ts = (n2*cosθt) / (n1*cosθi) * ts^2
    Tp = (n2*cosθt) / (n1*cosθi) * tp^2
    θt = asin(sqrt(sinθt²))
    (; rs, rp, ts, tp, Rs, Rp, Ts, Tp, θt)
end

brewster_angle(n1::Real, n2::Real) = atan(n2 / n1)

function critical_angle(n1::Real, n2::Real)
    n1 <= n2 && return NaN
    asin(n2 / n1)
end

# ── Electric multipoles (2-D, static) ────────────────────────────────────────

"""
    potential_2d(charges, positions, x, y; reg=0.03) → Float64

Electrostatic potential V = Σᵢ qᵢ / |r - rᵢ| (2-D field map, sources treated
as 3-D point charges evaluated in the xy-plane so V ∝ 1/r).
`reg` is a softening length to avoid the divergence at the source.
"""
function potential_2d(charges, positions, x::Real, y::Real; reg::Real=0.03)
    V = 0.0
    for (q, (px, py)) in zip(charges, positions)
        r = sqrt((x - px)^2 + (y - py)^2 + reg^2)
        V += q / r
    end
    V
end

"""
    efield_2d(charges, positions, x, y; reg=0.05) → (Ex, Ey)

Electric field components from a set of point charges at positions given as
a vector of (px, py) tuples.
"""
function efield_2d(charges, positions, x::Real, y::Real; reg::Real=0.05)
    Ex = Ey = 0.0
    for (q, (px, py)) in zip(charges, positions)
        dx = x - px;  dy = y - py
        r3 = ((dx^2 + dy^2 + reg^2)^1.5)
        Ex += q * dx / r3
        Ey += q * dy / r3
    end
    Ex, Ey
end

# ── Magnetic dipole (analytical, in xz-plane) ────────────────────────────────

"""
    bfield_dipole_xz(mz, x, z) → (Bx, Bz)

Magnetic field of a point dipole m = mz ẑ at the origin, evaluated in the xz-plane.
B = (μ₀/4π)[3(m·r̂)r̂ - m]/r³  with  μ₀/4π = 1.
"""
function bfield_dipole_xz(mz::Real, x::Real, z::Real; reg::Real=0.05)
    r2 = x^2 + z^2 + reg^2
    r  = sqrt(r2)
    r5 = r2^2.5
    m_dot_r = mz * z     # m = mz ẑ, r = (x, z)
    Bx = 3 * m_dot_r * x / r5
    Bz = (3 * m_dot_r * z - mz * r2) / r5
    Bx, Bz
end

# ── Biot-Savart — circular current loop in z=0 plane ─────────────────────────

"""
    bfield_loop_xz_v2(R, I, x, z; Nφ=600) → (Bx, Bz)

Magnetic field of a circular loop of radius R carrying current I in the z=0
plane, evaluated at a point (x, 0, z). Uses explicit Biot-Savart integration
with full 3-D cross product; By is discarded by azimuthal symmetry at y=0.
"""
function bfield_loop_xz_v2(R::Real, I::Real, x::Real, z::Real; Nφ::Int=600)
    Bx = By = Bz = 0.0
    dφ = 2π / Nφ
    for n in 0:(Nφ-1)
        φ  = (n + 0.5) * dφ
        lx = R * cos(φ);  ly = R * sin(φ)
        dlx = -R * sin(φ) * dφ;  dly = R * cos(φ) * dφ
        rx = x - lx;  ry = -ly;  rz = z
        r3 = (rx^2 + ry^2 + rz^2)^1.5
        # dl × r  (dlz = 0)
        Bx += I * (dly * rz) / r3
        By += I * (-dlx * rz) / r3
        Bz += I * (dlx * ry - dly * rx) / r3
    end
    Bx, Bz
end

# ── Hertzian dipole (oscillating electric dipole) ────────────────────────────
#
# An oscillating dipole p = p₀ cos(ωt) ẑ at the origin.
# Use retarded fields in SI-like form with c=1, k=ω.
#
# E_r = 2cosθ [cos(kr−ωt)/r³ − k sin(kr−ωt)/r²]
# E_θ = sinθ  [cos(kr−ωt)/r³ − k sin(kr−ωt)/r²  − k² cos(kr−ωt)/r] (×−1)
# B_φ = sinθ  [−k sin(kr−ωt)/r²  − k² cos(kr−ωt)/r]  (×−1 for usual sign conv.)
#
# In Cartesian xz coordinates (y=0 plane, dipole along z):

"""
    hertzian_fields(x, z, t; k=2π, p0=1.0, reg=0.05) → (Ex, Ez)

Retarded electric field in the xz-plane from a Hertzian dipole p(t) = p0 cos(ωt) ẑ.
"""
function hertzian_fields(x::Real, z::Real, t::Real;
                          k::Real=2π, p0::Real=1.0, reg::Real=0.05)
    r  = sqrt(x^2 + z^2 + reg^2)
    cosθ = z / r
    sinθ = x / r        # positive x half
    kr   = k * r
    ω    = k            # c = 1

    # Phase factors
    c1 = cos(kr - ω*t)
    s1 = sin(kr - ω*t)

    # Radial and polar components (p0 factor absorbed)
    Er =  2cosθ * p0 * (c1/r^3 - k*s1/r^2)
    Eθ = -sinθ  * p0 * (c1/r^3 - k*s1/r^2 - k^2*c1/r)

    # Convert (r̂, θ̂) → (x̂, ẑ)
    # r̂ = (sinθ, cosθ),  θ̂ = (cosθ, −sinθ)  in (x, z)
    Ex = Er * sinθ + Eθ * cosθ
    Ez = Er * cosθ - Eθ * sinθ
    Ex, Ez
end

"""
    hertzian_power_pattern(θs; k=2π) → Vector{Float64}

Radiation power pattern dP/dΩ ∝ sin²θ for a Hertzian dipole.
"""
hertzian_power_pattern(θs; k::Real=2π) = sin.(θs).^2
