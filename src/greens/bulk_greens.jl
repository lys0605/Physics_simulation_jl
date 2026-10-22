# Bulk Green's function utilities for the QWZ model (periodic BC in x and y).
#
# The single-particle retarded Green's function at crystal momentum k is
#
#   G(k, ω) = [(ω + iη) I − H(k)]⁻¹
#
# where η > 0 is a positive infinitesimal (Lorentzian broadening).
# The physical η → 0⁺ limit gives the spectral representation.
#
# Includes: retarded GF, spectral function A(k,ω), DOS, Chern number (FHS).

using LinearAlgebra

# ─────────────────────────── Green's function ─────────────────────────────────

"""
    retarded_green(model, kx, ky, ω; η=0.02) → 2×2 ComplexF64 matrix

Retarded Green's function:  G(k,ω) = [(ω+iη)I − H(k)]⁻¹.
"""
function retarded_green(model::QWZModel, kx::Real, ky::Real, ω::Real; η::Real=0.02)
    H = bloch_hamiltonian(model, kx, ky)
    inv(complex(ω, η) * I - H)
end

"""
    spectral_function(model, kx, ky, ω; η=0.02) → Float64

Momentum-resolved spectral function (summed over bands):

  A(k, ω) = −(1/π) Im Tr G(k, ω)
           = (1/π) Σₙ η / [(ω − Eₙ(k))² + η²]

The Lorentzian peak at each band Eₙ(k) has FWHM 2η.
"""
function spectral_function(model::QWZModel, kx::Real, ky::Real, ω::Real; η::Real=0.02)
    G = retarded_green(model, kx, ky, ω; η)
    -imag(tr(G)) / π
end

# Fast closed-form version using the 2-band eigenvalues.
function spectral_function_fast(model::QWZModel, kx::Real, ky::Real, ω::Real; η::Real=0.02)
    E₋, E₊ = band_energies(model, kx, ky)
    (η / ((ω - E₋)^2 + η^2) + η / ((ω - E₊)^2 + η^2)) / π
end

# ───────────────────────────── Density of States ──────────────────────────────

"""
    dos(model, ωs; Nk=100, η=0.02) → Vector{Float64}

Density of states  ρ(ω) = (1/N_k²) Σ_k A(k,ω),
integrated over the full Brillouin zone on an Nk×Nk uniform mesh.
"""
function dos(model::QWZModel, ωs; Nk::Int=100, η::Real=0.02)
    kxs = range(0, 2π; length=Nk+1)[1:Nk]
    kys = range(0, 2π; length=Nk+1)[1:Nk]
    Nω  = length(ωs)
    ρ   = zeros(Float64, Nω)
    for kx in kxs, ky in kys
        E₋, E₊ = band_energies(model, kx, ky)
        for (iω, ω) in enumerate(ωs)
            ρ[iω] += η / ((ω-E₋)^2 + η^2) + η / ((ω-E₊)^2 + η^2)
        end
    end
    ρ ./ (π * Nk^2)
end

# ───────────────── High-symmetry path spectral function ──────────────────────

"""
    kpath_points(N_per_seg=60) → (kxs, kys, xticks, xlabels)

High-symmetry path  Γ → X → M → Γ  on the square-lattice Brillouin zone,
with N_per_seg k-points per segment.

  Γ = (0,0),  X = (π,0),  M = (π,π)
"""
function kpath_points(N_per_seg::Int=60)
    seg1_kx = range(0,   π; length=N_per_seg)
    seg1_ky = fill(0.0, N_per_seg)
    seg2_kx = fill(π,   N_per_seg)
    seg2_ky = range(0,   π; length=N_per_seg)
    seg3_kx = range(π,   0; length=N_per_seg+1)[1:end-1]
    seg3_ky = range(π,   0; length=N_per_seg+1)[1:end-1]

    kxs = vcat(seg1_kx, seg2_kx, seg3_kx)
    kys = vcat(seg1_ky, seg2_ky, seg3_ky)

    ticks  = [1, N_per_seg, 2N_per_seg, 3N_per_seg]
    labels = ["Γ", "X", "M", "Γ"]
    kxs, kys, ticks, labels
end

# ─────────────────────── Chern number (Fukui–Hatsugai–Suzuki) ────────────────

"""
    chern_number(model; Nk=60) → Int

Topological Chern number of the lower band computed by the lattice gauge
(Fukui–Hatsugai–Suzuki 2005) formula:

  C = (1/2π) Σ_k Im ln [U_x(k) U_y(k+x̂) U_x†(k+ŷ) U_y†(k)]

where  U_μ(k) = ⟨u(k)|u(k+δμ)⟩/|⟨u(k)|u(k+δμ)⟩|  is the link variable
and |u(k)⟩ is the occupied (lower) Bloch eigenstate.

This expression is gauge-invariant and always yields an exact integer.
"""
function chern_number(model::QWZModel; Nk::Int=60)
    δ  = 2π / Nk
    ks = range(0, 2π; length=Nk+1)[1:Nk]

    # Cache all lower-band eigenvectors on the mesh.
    u = [begin
             H = bloch_hamiltonian(model, kx, ky)
             vecs = eigvecs(H)        # columns sorted by eigenvalue (ascending)
             vecs[:, 1]               # lower band
         end
         for kx in ks, ky in ks]

    C = 0.0
    for i in 1:Nk, j in 1:Nk
        ip = mod1(i+1, Nk);  jp = mod1(j+1, Nk)
        Ux  = _link(u[i,  j],  u[ip, j ])
        Uy  = _link(u[ip, j],  u[ip, jp])
        UxT = _link(u[i,  jp], u[ip, jp])
        UyT = _link(u[i,  j],  u[i,  jp])
        F   = Ux * Uy * conj(UxT) * conj(UyT)
        C  += angle(F)
    end
    round(Int, C / (2π))
end

_link(a, b) = dot(a, b) / abs(dot(a, b))
