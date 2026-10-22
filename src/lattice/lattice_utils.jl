# Bravais lattices — primitive / conventional vectors and reciprocal lattices.
#
# A Bravais lattice is the set { Σᵢ nᵢ aᵢ : nᵢ ∈ ℤ } where {aᵢ} are linearly
# independent primitive vectors.  A "conventional" cell may use a larger basis
# (e.g. the cubic edge for BCC/FCC) to reflect the symmetry, at the cost of
# putting more than one lattice site inside the cell.

using LinearAlgebra

# ── 2-D Bravais lattices ─────────────────────────────────────────────────────

"""
    bravais_2d(name; a=1.0, ratio=1.5, angle_deg=70.0) → (a1, a2, a1c, a2c, label)

Primitive vectors `a1, a2` and conventional vectors `a1c, a2c` for one of the
five 2-D Bravais lattices.  For all except `:centered_rectangular` the
conventional cell equals the primitive cell.
"""
function bravais_2d(name::Symbol; a::Real=1.0, ratio::Real=1.5,
                     angle_deg::Real=70.0)
    if name === :oblique
        θ = deg2rad(angle_deg)
        a1 = [a, 0.0];  a2 = [a*ratio*cos(θ), a*ratio*sin(θ)]
        a1, a2, a1, a2, "Oblique"
    elseif name === :rectangular
        a1 = [a, 0.0];  a2 = [0.0, a*ratio]
        a1, a2, a1, a2, "Rectangular"
    elseif name === :centered_rectangular
        # primitive: (a, 0) and (a/2, b/2) with b = a*ratio
        b = a * ratio
        a1 = [a, 0.0];  a2 = [a/2, b/2]
        # conventional: (a, 0) and (0, b)
        a1c = [a, 0.0];  a2c = [0.0, b]
        a1, a2, a1c, a2c, "Centered rectangular"
    elseif name === :square
        a1 = [a, 0.0];  a2 = [0.0, a]
        a1, a2, a1, a2, "Square"
    elseif name === :hexagonal
        a1 = [a, 0.0];  a2 = [a/2, a*sqrt(3)/2]
        a1, a2, a1, a2, "Hexagonal (triangular)"
    else
        error("unknown 2-D lattice: $name")
    end
end

"""
    reciprocal_2d(a1, a2) → (b1, b2)

Reciprocal basis defined by bᵢ · aⱼ = 2π δᵢⱼ.  In 2-D, identify the lattice
plane with a horizontal slice of 3-D so that 'cross products' reduce to a
scalar (signed area).
"""
function reciprocal_2d(a1::AbstractVector, a2::AbstractVector)
    crossz = a1[1]*a2[2] - a1[2]*a2[1]
    b1 = 2π * [ a2[2], -a2[1]] / crossz
    b2 = 2π * [-a1[2],  a1[1]] / crossz
    b1, b2
end

"""
    enumerate_2d(a1, a2, n) → Vector{Tuple{Float64,Float64}}

All lattice points Σ nᵢ aᵢ with |nᵢ| ≤ n.
"""
function enumerate_2d(a1, a2, n::Int)
    pts = NTuple{2,Float64}[]
    for i in -n:n, j in -n:n
        p = i .* a1 .+ j .* a2
        push!(pts, (p[1], p[2]))
    end
    pts
end

# ── 3-D Bravais lattices (cubic family + hexagonal) ──────────────────────────

"""
    bravais_3d(name; a=1.0, c=nothing) → NamedTuple

Returns (; prim, conv, basis_conv, label) where
  prim  = (a1, a2, a3) primitive vectors
  conv  = (A1, A2, A3) conventional vectors (cubic / hexagonal cell)
  basis_conv = basis sites inside conventional cell, in fractional coords
"""
function bravais_3d(name::Symbol; a::Real=1.0, c::Union{Nothing,Real}=nothing)
    if name === :sc
        a1 = [a,0,0];  a2 = [0,a,0];  a3 = [0,0,a]
        (; prim = (a1, a2, a3),
           conv = (a1, a2, a3),
           basis_conv = [(0.0, 0.0, 0.0)],
           label = "Simple cubic")
    elseif name === :bcc
        a1 = a/2 .* [-1,  1,  1]
        a2 = a/2 .* [ 1, -1,  1]
        a3 = a/2 .* [ 1,  1, -1]
        A1 = [a,0,0];  A2 = [0,a,0];  A3 = [0,0,a]
        (; prim = (a1, a2, a3),
           conv = (A1, A2, A3),
           basis_conv = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
           label = "Body-centered cubic")
    elseif name === :fcc
        a1 = a/2 .* [0, 1, 1]
        a2 = a/2 .* [1, 0, 1]
        a3 = a/2 .* [1, 1, 0]
        A1 = [a,0,0];  A2 = [0,a,0];  A3 = [0,0,a]
        (; prim = (a1, a2, a3),
           conv = (A1, A2, A3),
           basis_conv = [(0.0, 0.0, 0.0),
                          (0.5, 0.5, 0.0),
                          (0.5, 0.0, 0.5),
                          (0.0, 0.5, 0.5)],
           label = "Face-centered cubic")
    elseif name === :hex
        cc = isnothing(c) ? a * sqrt(8/3) : c    # ideal c/a for hcp
        a1 = [a, 0.0, 0.0]
        a2 = [a/2, a*sqrt(3)/2, 0.0]
        a3 = [0.0, 0.0, cc]
        (; prim = (a1, a2, a3),
           conv = (a1, a2, a3),
           basis_conv = [(0.0, 0.0, 0.0)],
           label = "Hexagonal")
    else
        error("unknown 3-D lattice: $name")
    end
end

"""
    reciprocal_3d(a1, a2, a3) → (b1, b2, b3)

bᵢ = 2π (aⱼ × aₖ) / V  with cyclic (i,j,k) and V = a₁·(a₂×a₃).
"""
function reciprocal_3d(a1, a2, a3)
    V  = dot(a1, cross(a2, a3))
    b1 = 2π * cross(a2, a3) / V
    b2 = 2π * cross(a3, a1) / V
    b3 = 2π * cross(a1, a2) / V
    b1, b2, b3
end

"""
    enumerate_3d(a1, a2, a3, n) → Vector{NTuple{3,Float64}}

Lattice points with index in [-n, n].  Origin included.
"""
function enumerate_3d(a1, a2, a3, n::Int)
    pts = NTuple{3,Float64}[]
    for i in -n:n, j in -n:n, k in -n:n
        p = i .* a1 .+ j .* a2 .+ k .* a3
        push!(pts, (p[1], p[2], p[3]))
    end
    pts
end
