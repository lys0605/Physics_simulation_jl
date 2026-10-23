# Bravais Lattices, Voronoi Cells, Wigner-Seitz Cells, and Brillouin Zones

> **Interactive companion:** open the Physics Atlas scene
> `cmp.brillouin-zones` (repository lys0605/physics-atlas) to
> rotate any of the 14 Bravais lattices, switch between primitive and
> conventional cells, and watch the Brillouin zone update live as you
> change parameters (`a`, `c`, `α`, …).  Useful while reading this note —
> change to a rhombohedral lattice and sweep `α` from 60° to 109.47° to
> watch the BZ morph from a truncated octahedron (FCC) through a cube (SC)
> to a rhombic dodecahedron (BCC).

This note ties together four ideas that crystallographers, solid-state
physicists, and computational geometers all use under different names:

| Concept            | Domain               | Definition (informal)                              |
|--------------------|----------------------|----------------------------------------------------|
| Voronoi cell       | computational geometry | points closer to a given seed than to any other  |
| Wigner-Seitz cell  | real-space crystallography | Voronoi cell of a Bravais lattice site       |
| Brillouin zone     | momentum-space crystallography | Wigner-Seitz cell of the reciprocal lattice |
| Half-space intersection | convex geometry | the unifying construction                          |

We work out (i) what each object is, (ii) why they coincide, (iii) the
algorithms used to compute them, and (iv) why the resulting cells matter
physically.

---

## 1. Bravais lattices: primitive and conventional cells

A **Bravais lattice** in `d` dimensions is the set

$$
\Lambda \;=\;\bigl\{\, n_1\mathbf{a}_1 + n_2\mathbf{a}_2 + \cdots + n_d\mathbf{a}_d \;:\; n_i \in \mathbb{Z}\,\bigr\}
$$

where `a₁, …, a_d` are linearly independent **primitive vectors**.
The parallelepiped they span (volume `V_p = |\det[a_1\,a_2\,\cdots\,a_d]|`)
is the **primitive cell**.  By construction it contains exactly one lattice
point and tiles space under translation by `Λ`.

A **conventional cell** is any unit cell — primitive or not — chosen to
reveal the point-group symmetry of the lattice.  Examples:

* For body-centered cubic (BCC), `a₁ = (a/2)(-x̂+ŷ+ẑ)` (etc.) gives a primitive
  rhombohedron of volume `a³/2`.  The cubic cell of edge `a` is the
  conventional choice, holds 2 lattice sites, and makes the cubic symmetry
  obvious.
* For face-centered cubic (FCC), the primitive volume is `a³/4`; the
  conventional cubic cell holds 4 sites.
* For all 2-D Bravais lattices except *centered rectangular*, the
  primitive cell already exhibits the full point symmetry; only the
  centered rectangular needs a conventional cell with 2 sites.

The **reciprocal lattice** `\Lambda^*` has primitive vectors `b_i` defined by

$$
\mathbf{b}_i \cdot \mathbf{a}_j \;=\; 2\pi\,\delta_{ij}.
$$

Explicitly, in three dimensions,

$$
\mathbf{b}_1 = \frac{2\pi\,\mathbf{a}_2\times\mathbf{a}_3}{V_p}, \qquad
\mathbf{b}_2 = \frac{2\pi\,\mathbf{a}_3\times\mathbf{a}_1}{V_p}, \qquad
\mathbf{b}_3 = \frac{2\pi\,\mathbf{a}_1\times\mathbf{a}_2}{V_p}.
$$

Cubic-family duality (memorize these):

* simple cubic (SC) → reciprocal SC (edge `2π/a`)
* BCC → FCC reciprocal (edge `4π/a`)
* FCC → BCC reciprocal (edge `4π/a`)

The **reciprocal of the conventional cubic cell** is itself a SC reciprocal
of edge `2π/a` — it gives a smaller cube than the true primitive BZ for BCC
and FCC.  This conventional reciprocal is sometimes seen in band-folding
arguments but is *not* the first Brillouin zone.

---

## 2. Voronoi cells

For a finite point set `S \subset \mathbb{R}^d`, the **Voronoi cell** of
`p \in S` is

$$
V(p) \;=\; \bigl\{\, x \in \mathbb{R}^d \;:\; \lvert x - p\rvert \le \lvert x - q\rvert \;\;\forall\, q \in S,\, q\neq p\,\bigr\}.
$$

The defining inequality unpacks as a half-space:

$$
\lvert x - p\rvert^2 \le \lvert x - q\rvert^2
\;\Longleftrightarrow\;
(q - p)\cdot x \;\le\; \tfrac{1}{2}\bigl(\lvert q\rvert^2 - \lvert p\rvert^2\bigr).
$$

So `V(p) = \bigcap_{q\neq p} H_{p,q}` is an intersection of `|S|-1`
half-spaces — automatically convex and (if bounded) polyhedral.  Each
boundary segment lies on the perpendicular bisector of an edge `pq`.

**Useful properties:**

* The Voronoi diagram tiles `\mathbb{R}^d` without gaps or overlap.
* Cells of nearest neighbors share a face; their **dual** graph is the
  *Delaunay triangulation*.
* Bounded cells have at most `O(n^{\lceil d/2\rceil})` vertices (Upper
  Bound Theorem).

---

## 3. Wigner-Seitz cells

The **Wigner-Seitz (WS) cell** is the Voronoi cell of a lattice point in
its Bravais lattice.  Setting `p = 0` and using only the *nearest few
shells* of lattice neighbors `R_i \in \Lambda \setminus \{0\}`,

$$
\text{WS}\;=\;\bigl\{\, x \;:\; \mathbf{R}_i\cdot x \le \tfrac{1}{2}|\mathbf{R}_i|^2 \;\;\forall\, i\,\bigr\}.
$$

Because the lattice is translationally invariant, every WS cell is a
translate of this central one.  The result is the **unique primitive cell
with the full point symmetry of the lattice**: it inherits all the
rotations, reflections, and inversion that fix a lattice site.

Common 2-D and 3-D shapes:

| Lattice                              | WS shape                          |
|--------------------------------------|-----------------------------------|
| 2-D square                           | square                            |
| 2-D triangular (= hexagonal)         | regular hexagon                   |
| 2-D centered rectangular             | elongated hexagon                 |
| 3-D simple cubic (SC)                | cube                              |
| 3-D body-centered cubic (BCC)        | truncated octahedron (14 faces)   |
| 3-D face-centered cubic (FCC)        | rhombic dodecahedron (12 faces)   |

---

## 4. Brillouin zones

The **first Brillouin zone (BZ)** is the WS cell of the *reciprocal*
Bravais lattice.  Substituting `R_i → G_i` (reciprocal lattice vectors)
into the WS definition,

$$
\text{BZ}\;=\;\bigl\{\, \mathbf{k} \;:\; \mathbf{G}\cdot \mathbf{k} \le \tfrac{1}{2}|\mathbf{G}|^2 \;\;\forall\, \mathbf{G}\in\Lambda^*\setminus\{0\}\,\bigr\}.
$$

The half-space boundaries are the **Bragg planes**: surfaces in `k`-space
where the Laue diffraction condition

$$
2\,\mathbf{k}\cdot\mathbf{G} = |\mathbf{G}|^2
$$

is satisfied.  This is why BZ boundaries are the natural place for
crystal momentum to become degenerate (band gaps open at Bragg planes
in perturbation theory).

By reciprocal-lattice duality:

| Real-space lattice | Reciprocal-space lattice | Brillouin zone |
|--------------------|--------------------------|----------------|
| SC                 | SC                       | cube           |
| BCC                | FCC                      | rhombic dodecahedron |
| FCC                | BCC                      | truncated octahedron |
| 2-D square         | 2-D square               | square         |
| 2-D triangular     | 2-D triangular (rot. 30°)| hexagon        |

---

## 5. The unifying view

All four objects — Voronoi, WS, BZ, and any half-space intersection — are
**convex polytopes built from perpendicular bisectors**.  The only
distinctions are:

1. *Voronoi*: bisect with respect to an arbitrary finite point set.
2. *WS*: bisect with respect to lattice neighbors only.
3. *BZ*: bisect in reciprocal space.

So one algorithm — *intersect a bounding box with the half-spaces from
each bisector* — covers all three uses.  In particular, the WS algorithm
and the BZ algorithm are **the same code** applied to two different
lattices.

---

## 6. Algorithm — 2-D (Sutherland-Hodgman clipping)

For an origin-centered WS cell in 2-D:

```
poly ← large CCW bounding rectangle
for each neighbor R = (Rx, Ry) ≠ 0:
    normal ← (Rx, Ry)
    offset ← (Rx² + Ry²) / 2
    poly ← clip_polygon(poly, normal, offset)   # keep n·x ≤ offset
return poly
```

`clip_polygon` walks the polygon's edges and, for each edge `(p₁, p₂)`,
keeps the inside half:

```
for each edge (p1, p2):
    d1 = n · p1 - offset
    d2 = n · p2 - offset
    if both inside (d1, d2 ≤ 0):   keep p2
    if leaving (d1 ≤ 0, d2 > 0):    keep intersection
    if entering (d1 > 0, d2 ≤ 0):   keep intersection, then p2
    if both outside:                drop
```

Time complexity: each clipping operation is `O(|poly|)`, and the polygon
size stays `O(\text{shells used})`.  Two or three nearest-neighbor shells
suffice in every 2-D Bravais lattice.

---

## 7. Algorithm — 3-D (vertex enumeration from half-spaces)

In 3-D the polygon-clipping idea generalizes to *polyhedron clipping*,
but a simpler approach is **vertex enumeration**:

```
half-spaces:   n_i · x ≤ d_i,  i = 1…N         (n_i = R_i, d_i = |R_i|²/2)

vertices ← []
for each triple (i, j, k):
    solve  [n_i; n_j; n_k] v = (d_i, d_j, d_k)
    if singular (det = 0):          continue
    if any other constraint fails:  continue       (i.e. v lies outside polytope)
    add v to vertices  (with the set of planes it lies on)

for each plane f:
    collect vertices that lie on plane f
    sort them counterclockwise about the face centroid (using a 2-D
        basis (e1, e2) inside the plane)

faces[f] = sorted vertex indices on plane f
```

Cost: `O(N³)` triples × `O(N)` validation per triple ≈ `O(N⁴)`.  Since `N`
(the number of neighbor lattice points considered) is at most a few dozen
for nearest- and next-nearest-shell candidates, this remains fast even
without optimization.

Counterclockwise sorting on a face uses the local frame

```
n̂  = R_f / |R_f|
e₁ = n̂ × ẑ  (or × x̂ if n̂ ≈ ẑ),  normalized
e₂ = n̂ × e₁
angle(v) = atan2( (v - centroid)·e₂, (v - centroid)·e₁ )
```

then sort vertices by `angle`.

**Choosing how many neighbors to include.**  Safe rule of thumb: take
all lattice points within `R_max = 2 · |a_min|` where `a_min` is the
shortest primitive vector.  This always covers the nearest shell that
contributes a face, and the algorithm itself discards spurious vertices.

---

## 8. Physical significance

Several theorems and constructions live entirely on the BZ:

* **Bloch's theorem**.  Single-particle eigenstates of a periodic potential
  are labeled by a crystal momentum `\mathbf{k} \in \text{BZ}`:

  $$
  \psi_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}}\, u_{n\mathbf{k}}(\mathbf{r}),
  \qquad u_{n\mathbf{k}}(\mathbf{r}+\mathbf{R}) = u_{n\mathbf{k}}(\mathbf{r}).
  $$

  Two momenta `\mathbf{k}, \mathbf{k}+\mathbf{G}` (`G \in \Lambda^*`) label
  the same eigenstate, so all distinct states live within one BZ.

* **Bragg planes and gaps.**  At a BZ boundary, two plane waves
  `\mathbf{k}` and `\mathbf{k}-\mathbf{G}` are degenerate at the
  free-electron level.  Any non-zero Fourier component
  `V_{\mathbf{G}}` of the periodic potential lifts the degeneracy and
  opens a band gap.  Energy gaps therefore live on the BZ skeleton.

* **Topological invariants.**  Chern numbers and other invariants are
  integrals over the BZ (a compact torus once opposite faces are
  identified):

  $$
  C_n = \frac{1}{2\pi}\!\int_{\text{BZ}}\!\Omega_n(\mathbf{k})\, d^2k,
  $$

  where `Ω_n` is the Berry curvature of band `n`.

* **Phonons and electrons share the BZ.**  Crystal momentum
  conservation in scattering processes is conservation **modulo `G`** —
  Umklapp processes — which is exactly the BZ-periodicity statement.

---

## 9. Connections to the implementation here

The simulation files mirror each conceptual layer:

| File                                          | Concept                  |
|-----------------------------------------------|--------------------------|
| `src/lattice/lattice_utils.jl`                | Bravais lattices, primitive ↔ conventional, reciprocal |
| `src/lattice/voronoi.jl`                      | half-space intersection — 2-D polygon clip + 3-D vertex enumeration |
| `scripts/lattice/lattices_2d.jl`              | five 2-D Bravais lattices with both cell choices |
| moved to Physics Atlas (`cmp.brillouin-zones`) | SC, BCC, FCC and the other 11 lattices — conventional vs primitive |
| `scripts/lattice/voronoi.jl`                  | general Voronoi diagram (random) → WS cell for a lattice |
| `scripts/lattice/brillouin.jl`                | WS (real) vs BZ (reciprocal), 2-D and 3-D |

Computational notes:

* All 2-D cells use Sutherland-Hodgman clipping (`clip_polygon`,
  `voronoi_cell_2d`, `ws_cell_2d`).  Robust against degenerate cases via
  a small tolerance.
* All 3-D cells use vertex enumeration (`ws_cell_3d`).  Faces are sorted
  CCW so they can be plotted as closed polygons.
* The reciprocal lattice is computed via the `2π (aⱼ × aₖ) / V` formula
  (`reciprocal_3d`); for 2-D, the matrix inverse reduces to a scalar
  cross product (`reciprocal_2d`).
* Computing the BZ is **literally the WS routine** called on reciprocal
  lattice neighbors — emphasizing that real-space and reciprocal-space
  cells share one algorithm.
