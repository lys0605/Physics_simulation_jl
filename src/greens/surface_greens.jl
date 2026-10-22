# Surface Green's function for a semi-infinite QWZ slab.
#
# Geometry: periodic BC in y (k_y is a good quantum number); the system is
# semi-infinite in x — layers j = 1, 2, 3, … extending into the bulk.
# The surface is at j = 1.
#
# The retarded surface GF satisfies the Dyson self-consistency equation
#
#   G_s(k_y, ω) = [(ω+iη)I − H₀(k_y) − T G_s T†]⁻¹        (*)
#
# where H₀ = layer_onsite, T = interlayer_hopping().
#
# Two solvers are provided:
#
#   surface_green_naive    — fixed-point iteration of (*); simple but O(1/η)
#   surface_green          — López Sancho decimation (1984); O(log(1/η)) steps
#
# The surface spectral function
#
#   A_s(k_y, ω) = −(1/π) Im Tr G_s(k_y, ω)
#
# shows the surface-projected bands, including any in-gap edge modes.

using LinearAlgebra

# ──────────────────────── Simple fixed-point iteration ───────────────────────

"""
    surface_green_naive(model, ky, ω; η, tol, max_iter) → 2×2 ComplexF64

Fixed-point iteration of the Dyson equation (*).  Convergence requires
‖ΔG‖ < tol.  For small η this needs O(1/η²) iterations — use only for
large η or as a reference.
"""
function surface_green_naive(model::QWZModel, ky::Real, ω::Real;
                              η::Real       = 0.05,
                              tol::Real     = 1e-8,
                              max_iter::Int = 5_000)
    z  = complex(ω, η)
    H0 = layer_onsite(model, ky)
    T  = interlayer_hopping()
    Gs = inv(z * I - H0)           # zeroth-order: isolated layer
    for _ in 1:max_iter
        Gs_new = inv(z * I - H0 - T * Gs * T')
        norm(Gs_new - Gs) < tol && return Gs_new
        Gs = Gs_new
    end
    @warn "naive iteration did not converge (η=$η, ω=$ω)"
    Gs
end

# ──────────────────── López Sancho decimation algorithm ──────────────────────

"""
    surface_green(model, ky, ω; η, tol, max_iter) → 2×2 ComplexF64

Retarded surface Green's function via the López Sancho (1984) renormalisation
scheme.  The algorithm doubles the effective slab depth at each step, giving
O(log₂(1/η)) iterations instead of O(1/η²).

Algorithm (initialise with ε_s = ε = H₀, α = T, β = T†):

  for each iteration:
    g    = (z I − ε)⁻¹
    ε_s ← ε_s + α g β
    ε   ← ε   + α g β + β g α
    α   ← α g α
    β   ← β g β
  until ‖α‖ < tol

  G_s = (z I − ε_s)⁻¹
"""
function surface_green(model::QWZModel, ky::Real, ω::Real;
                        η::Real       = 0.02,
                        tol::Real     = 1e-8,
                        max_iter::Int = 200)
    z  = complex(ω, η)
    H0 = layer_onsite(model, ky)
    T  = interlayer_hopping()

    εs = copy(H0)           # surface effective Hamiltonian
    ε  = copy(H0)           # bulk effective Hamiltonian
    α  = copy(T)            # forward (into-bulk) coupling
    β  = T'                 # backward coupling

    for _ in 1:max_iter
        g  = inv(z * I - ε)
        εs = εs + α  * g * β
        ε  = ε  + α  * g * β + β * g * α
        α  = α * g * α
        β  = β * g * β
        norm(α) < tol && break
    end
    inv(z * I - εs)
end

# ──────────────────────── Spectral function helpers ──────────────────────────

"""
    surface_spectral(model, ky, ω; kwargs...) → Float64

Surface spectral function A_s(k_y, ω) = −(1/π) Im Tr G_s(k_y, ω).
"""
function surface_spectral(model::QWZModel, ky::Real, ω::Real; kwargs...)
    -imag(tr(surface_green(model, ky, ω; kwargs...))) / π
end

"""
    surface_spectral_map(model, kys, ωs; η, kwargs...) → Matrix{Float64}

Compute A_s on a (k_y × ω) grid.  Returns an Nky × Nω matrix.
"""
function surface_spectral_map(model::QWZModel, kys, ωs; η::Real=0.02, kwargs...)
    A = zeros(Float64, length(kys), length(ωs))
    for (j, ky) in enumerate(kys), (i, ω) in enumerate(ωs)
        A[j, i] = surface_spectral(model, ky, ω; η, kwargs...)
    end
    A
end

# ───────────────── Convergence diagnostic ────────────────────────────────────

"""
    lopezSancho_convergence(model, ky, ω; η) → (iters, α_norms)

Return the ‖α‖ per iteration for a single (k_y, ω) point.
Useful for understanding how many iterations are actually needed.
"""
function lopezSancho_convergence(model::QWZModel, ky::Real, ω::Real;
                                  η::Real=0.02, tol::Real=1e-10, max_iter::Int=200)
    z  = complex(ω, η)
    H0 = layer_onsite(model, ky)
    T  = interlayer_hopping()
    εs = copy(H0);  ε = copy(H0)
    α  = copy(T);   β = T'
    norms = Float64[]
    for _ in 1:max_iter
        push!(norms, norm(α))
        norm(α) < tol && break
        g  = inv(z * I - ε)
        εs = εs + α * g * β
        ε  = ε  + α * g * β + β * g * α
        α  = α * g * α
        β  = β * g * β
    end
    norms
end
