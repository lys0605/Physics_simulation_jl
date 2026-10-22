# Related Numerical Methods — When MC Is *Not* the Right Tool

Monte Carlo is one of several pillars of computational many-body physics.
This note sketches the others, with their **principle**, **scope**, and
**limits**, so you can choose the right tool for a problem and so we have
a roadmap for future directories in this repository.

## 1. Exact Diagonalisation (ED)

**Principle.** Write the Hamiltonian as a (sparse) matrix in a chosen
basis and compute eigenvalues / eigenvectors directly.

* Full diagonalisation (LAPACK `eig`): all eigenvalues, `O(D^3)` time,
  `O(D^2)` memory.  Feasible up to `D \sim 10^4`.
* Iterative methods (Lanczos, Arnoldi, Davidson): a few extreme
  eigenvalues by repeated matrix-vector products.  Feasible up to
  `D \sim 10^9` with care.

Symmetries (translation, point group, particle number, spin) **block-
diagonalise** the matrix and reduce `D` by orders of magnitude.

**Lanczos in one box.**

```
v_1 ← random unit vector
β_0 ← 0,  v_0 ← 0
for k = 1, 2, …, K:
    w   = H v_k − β_{k-1} v_{k-1}
    α_k = ⟨v_k, w⟩
    w   = w − α_k v_k
    β_k = ||w||
    v_{k+1} = w / β_k
return tridiagonal T_K with α_k on diagonal and β_k on sub-/super-diagonal
```

Diagonalise the small `T_K` for approximate eigenvalues of `H`.

**Scope.** Spin-½ on `\le 30` sites; fermions on `\sim 18-22` sites; up to
hundreds of orbitals with strong symmetries.  Gold standard for spectra
and dynamics on small systems.

**Limits.** Hilbert-space dimension grows exponentially.  No path forward
beyond ~30 spin-½ sites without compression.

**Where it shines.** Validation of QMC / DMRG / tensor-network codes on
small systems.  Real-time dynamics with no sign problem.  Computing
arbitrary correlators.

## 2. Density-Matrix Renormalisation Group (DMRG) and Matrix Product States

**Principle.** Approximate the wavefunction as a **matrix product state**:

$$
|\Psi\rangle = \sum_{\{\sigma\}} \text{Tr}\!\left[\,A^{\sigma_1}_{1} A^{\sigma_2}_{2} \cdots A^{\sigma_L}_{L}\right] |\sigma_1 \sigma_2 \cdots \sigma_L\rangle,
$$

with each `A^{\sigma}_i` a `\chi \times \chi` matrix.  The **bond dimension**
`\chi` controls accuracy — for gapped 1-D systems, `\chi` need only be
*constant* in `L` (entanglement area law).

**DMRG sweep** (White, 1992):

```
1. Initialize random MPS at bond dimension χ.
2. Sweep left↔right; at each bond, build the effective Hamiltonian
   H_eff on a two-site block; diagonalise; SVD to recover MPS form,
   truncating to χ singular values.
3. Continue until energy converges.
```

DMRG = variational optimisation over the MPS manifold.

**Scope.** 1-D and quasi-1-D systems (ladders, narrow cylinders): the
state of the art.  Gapped 2-D systems on cylinders up to width ~10.  Spin
chains, fermions, bosons, gauge theories.

**Limits.** Entanglement grows linearly with subsystem perimeter (area
law); in `\ge 2`-D the required bond dimension grows exponentially with
the smaller transverse dimension.

**Where it shines.** Ground states and low-lying spectra of 1-D systems
with no sign problem.  DMRG complements QMC: where QMC has a sign
problem, DMRG often does not — and vice versa for high-dimensional
sign-free models.

**Related extensions.**

* **TEBD** (time-evolving block decimation) — real-time evolution of MPS.
* **TDVP** — time-dependent variational principle, geometric integration
  on the MPS manifold.
* **Finite-temperature DMRG** — purification + imaginary-time TEBD.

## 3. Projected Entangled-Pair States (PEPS) — 2-D tensor networks

Generalisation of MPS to 2-D.  Each site carries a tensor with `d`
physical and four virtual indices of dimension `D`.  Norms and
expectation values are evaluated by **boundary MPS** or **CTM**
algorithms — themselves approximate.

**Scope.** 2-D ground states with area-law entanglement, including
frustrated antiferromagnets where QMC fails.  Hubbard-like fermion
models (fPEPS).

**Limits.** Contraction is `\#P`-hard; bond dimensions are kept small
(`D \le 10` typical) and accuracy is harder to certify than DMRG.

## 4. MERA and HOTRG

* **MERA** (Multi-scale Entanglement Renormalisation Ansatz): hierarchical
  tensor network respecting scale invariance; captures critical (log-
  violating) entanglement.
* **HOTRG** (Higher-Order Tensor Renormalisation Group): coarse-grain
  classical or quantum partition functions in any `D`.

Niche but conceptually important; relevant for critical theories and
holographic models.

## 5. Hartree-Fock, Density Functional Theory, Coupled Cluster

**Hartree-Fock** — minimise the energy of a single Slater determinant
self-consistently:

$$
F\,\phi_i = \epsilon_i\,\phi_i, \qquad F = h + \sum_j (J_j - K_j).
$$

Captures the *mean field*; misses correlation.

**DFT** — replace many-body wavefunction with the density:

$$
E[\rho] = T_s[\rho] + \int v_{\rm ext}\rho\,d^3r + E_H[\rho] + E_{xc}[\rho],
$$

solving Kohn-Sham equations `(-\tfrac12\nabla^2 + v_{\rm eff})\phi_i = \epsilon_i\phi_i`.
The exchange-correlation functional `E_{xc}` is approximate (LDA, GGA,
hybrids).  Computational workhorse for materials and chemistry.

**Coupled Cluster** — exponential ansatz `|\Psi\rangle = e^{\hat T}|\Phi_0\rangle`
with `\hat T = \hat T_1 + \hat T_2 + \cdots`.  CCSD(T) is the "gold standard"
for small molecules; scales as `N^7`.

**Scope vs MC.**  HF and DFT are *single-particle* methods — wrong for
strongly correlated systems but extremely fast.  CC is post-HF, captures
dynamical correlation, breaks for static correlation.  QMC is more
accurate than DFT and parallel to CC, often used as a benchmark.

## 6. Configuration Interaction (CI) and selected CI

CI: write `|\Psi\rangle = \sum_I c_I |D_I\rangle` over Slater determinants.
Full CI in a small basis is the exact solution but factorial cost.
**Selected CI** (CIPSI, ASCI, HCI) keeps only the determinants with
largest perturbative weight — bridge between CI and QMC.

## 7. Functional Renormalisation Group (fRG)

Flow equations for `n`-point vertex functions as a UV cutoff is lowered
from `\Lambda \to 0`.  One-loop fRG captures competing orders (SDW, SC,
CDW) in correlated electron systems without sampling.

**Wetterich equation.**

$$
\partial_\Lambda \Gamma_\Lambda \;=\; \tfrac{1}{2}\,\text{STr}\!\left[(\Gamma_\Lambda^{(2)} + R_\Lambda)^{-1}\,\partial_\Lambda R_\Lambda\right].
$$

Truncations introduce uncontrolled errors; results often qualitative,
useful for phase-diagram reconnaissance.

## 8. Dynamical Mean-Field Theory (DMFT)

Map a lattice model with local interactions to an effective **quantum
impurity** coupled self-consistently to a bath:

$$
G_{\rm imp}(i\omega_n) \;=\; G_{\rm loc}(i\omega_n),
$$

where `G_{\rm loc}` is the lattice local Green's function with self-energy
`\Sigma_{\rm imp}`.  Solve the impurity by ED, CT-HYB QMC, NRG, or DMRG.
Exact in `d=\infty`; the leading correlated-electron method for finite
`d`.

**Extensions.** Cluster DMFT (CDMFT, DCA) for short-range non-local
correlations; DMFT+DFT for *ab initio* correlated materials.

## 9. Numerical Renormalisation Group (NRG, Wilson)

Log-discretise the bath of an impurity model and iteratively diagonalise.
Exponentially accurate at low energies; the original method for the
Kondo problem.  Today often combined with DMFT.

## 10. Quantum embedding (DMET, RDMET)

**Density-matrix embedding theory.** Replace a fragment's environment by
a small *bath* whose 1-RDM matches the environment's projected 1-RDM.
Cheaper than DMFT, no `i\omega_n` axis.  Active development.

## 11. Quantum-computing methods (looking ahead)

* **VQE** — variational quantum eigensolver: parametrised quantum
  circuit + classical optimiser.
* **QPE** — quantum phase estimation: exact eigenvalues in polynomial
  depth on fault-tolerant hardware.
* **QITE** — quantum imaginary time evolution.

These are not classical numerical methods, but their classical
*simulations* are becoming routine and overlap with VMC / tensor-network
families.

## 12. How to choose — a cheat sheet

| Problem type | First-line method | Backup |
|---|---|---|
| 1-D spin / fermion ground state | DMRG | ED, QMC |
| 1-D dynamics (small `t`) | TEBD / TDVP | ED (small `L`) |
| 2-D sign-free spin model, large `L` | SSE | DMRG on cylinder |
| 2-D frustrated spin model | DMRG / PEPS | tensor networks |
| Continuum bosons (`^4`He) | PIMC | DMC |
| Small molecules (chemistry) | CCSD(T), DMC, AFQMC | selected CI |
| Solids, weakly correlated | DFT | post-DFT (GW), DMC |
| Solids, strongly correlated | DFT + DMFT (CT-HYB) | DMRG impurity solver |
| Single impurity, low-`T` thermodynamics | NRG | CT-HYB |
| Real-time dynamics, weak coupling | tDMRG, TDVP | exact diagonalisation |
| Lattice gauge theory at `\mu=0` | HMC | tensor networks |
| Lattice gauge theory at `\mu\neq 0` | (no good method yet) | complex Langevin, thimbles |

## 13. Where this directory could expand

Each method below could become its own simulation module — same
notes-plus-scripts style as this directory:

* `docs/ed/` — ED + Lanczos with symmetry adaptation.  Pair with the
  existing `src/greens/qwz_model.jl` to do small-cluster ED comparisons.
* `docs/dmrg/` — MPS, DMRG, TEBD.  Implementation in Julia via
  `ITensors.jl`.
* `docs/dmft/` — DMFT loop with an ED or CT-HYB solver.
* `docs/dft/` — basic plane-wave or Gaussian DFT (could leverage
  existing electronic structure libraries).
* `docs/tensor_networks/` — PEPS contraction, CTM, MERA.

Build out a directory only when you have a *specific* problem that the
method solves — premature abstraction here is as expensive as anywhere
else.
