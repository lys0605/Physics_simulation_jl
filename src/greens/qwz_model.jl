# Qi-Wu-Zhang (QWZ) model — a minimal 2D Chern insulator on a square lattice.
#
# The Bloch Hamiltonian is the 2×2 matrix
#
#   H(k) = d_x σ_x + d_y σ_y + d_z σ_z
#
# with   d_x = sin k_x,  d_y = sin k_y,  d_z = m + cos k_x + cos k_y.
#
# The real-space tight-binding decomposition is
#
#   H = Σ_j  H₀|j><j|  +  T|j+1><j|  +  T†|j><j+1|   (in x)
#           + T_y|j+ŷ><j| + T_y†|j><j+ŷ|              (in y)
#
# with on-site H₀ = m σ_z and hopping T = (σ_z − i σ_x)/2 in x (and
# T_y = (σ_z − i σ_y)/2 in y).  A single unit cell has 2 orbitals (spin-1/2).
#
# Topological phase diagram (Chern number C):
#   m < −2  →  C = 0  (trivial)
#   −2 < m < 0  →  C = +1
#   0 < m < 2   →  C = −1
#   m > 2   →  C = 0  (trivial)

using LinearAlgebra

struct QWZModel
    m::Float64   # on-site mass; controls topological phase
end

# ─────────────────────────── Pauli matrices ───────────────────────────────────

const σ₀ = ComplexF64[1 0; 0 1]
const σ₁ = ComplexF64[0 1; 1 0]
const σ₂ = ComplexF64[0 -im; im 0]
const σ₃ = ComplexF64[1 0; 0 -1]

# ─────────────────────────── Bloch Hamiltonian ────────────────────────────────

"""
    bloch_hamiltonian(model, kx, ky) → 2×2 ComplexF64 matrix

Bulk Bloch Hamiltonian H(k) of the QWZ model.
"""
function bloch_hamiltonian(model::QWZModel, kx::Real, ky::Real)
    dz = model.m + cos(kx) + cos(ky)
    dx = sin(kx)
    dy = sin(ky)
    ComplexF64[dz         dx - im*dy;
               dx + im*dy    -dz    ]
end

"""
    band_energies(model, kx, ky) → (E₋, E₊)

Closed-form eigenvalues of H(k): E± = ±|d(k)|.
"""
function band_energies(model::QWZModel, kx::Real, ky::Real)
    dz = model.m + cos(kx) + cos(ky)
    dx = sin(kx)
    dy = sin(ky)
    E  = sqrt(dx^2 + dy^2 + dz^2)
    (-E, E)
end

# ──────────────── Per-layer quantities for ribbon / surface GF ────────────────

"""
    layer_onsite(model, ky) → 2×2 matrix

On-site + y-direction hopping summed into a single layer Hamiltonian for fixed k_y.
This is the diagonal block of the ribbon Hamiltonian:

  H₀(k_y) = (m + cos k_y) σ_z + sin k_y σ_y
"""
function layer_onsite(model::QWZModel, ky::Real)
    ComplexF64[(model.m + cos(ky))   -im*sin(ky);
               im*sin(ky)            -(model.m + cos(ky))]
end

"""
    interlayer_hopping() → 2×2 matrix T

Hopping matrix from layer j to layer j+1 in the +x direction:

  T = (σ_z − i σ_x) / 2  =  [0.5  −0.5i; −0.5i  −0.5]

The reverse hopping is T† (Hermitian conjugate).  Verification:
  T e^{ik_x} + T† e^{−ik_x} = cos(k_x) σ_z + sin(k_x) σ_x   ✓
"""
function interlayer_hopping()
    ComplexF64[0.5   -0.5im;
               -0.5im  -0.5]
end

# ──────────────────── Ribbon Hamiltonian (finite slab) ────────────────────────

"""
    ribbon_hamiltonian(model, ky, Nx) → 2Nx × 2Nx Hermitian matrix

Build the strip Hamiltonian for Nx layers in the x-direction at fixed k_y.
Open (hard-wall) boundary conditions in x; periodic BC in y is encoded in
the k_y dependence of the on-site block.

Block structure (each block is 2×2):

  ⌈ H₀   T    0  ⋯ ⌉
  | T†  H₀   T  ⋯ |
  | 0   T†  H₀  ⋯ |
  ⌊ ⋮              ⌋
"""
function ribbon_hamiltonian(model::QWZModel, ky::Real, Nx::Int)
    H  = zeros(ComplexF64, 2Nx, 2Nx)
    H0 = layer_onsite(model, ky)
    T  = interlayer_hopping()
    for j in 1:Nx
        r = (2j-1):2j
        H[r, r] = H0
        if j < Nx
            s = (2j+1):2(j+1)
            H[r, s] = T        # hop right: ⟨j|H|j+1⟩
            H[s, r] = T'       # hop left:  ⟨j+1|H|j⟩
        end
    end
    H
end
