# Quantum Monte Carlo

Classical MC samples Boltzmann weights `e^{-\beta H(x)}` directly.  Quantum
problems contain *operators* and *non-commuting* terms — the Boltzmann
weight `\langle x|e^{-\beta\hat H}|x\rangle` is not in general non-negative, and the
trace `Z = \text{Tr}\,e^{-\beta\hat H}` cannot be sampled as written.

QMC methods get around this by trading the quantum problem for a higher-
dimensional classical sampling problem.  Three master strategies:

| Family | Idea | Examples |
|---|---|---|
| **Wavefunction-based** | Sample `|\Psi|^2` directly | VMC, DMC |
| **Path-integral** | Map to `d+1`-D classical | PIMC, world-line QMC |
| **Series expansion** | Sample a stochastic expansion of `Z` | SSE |
| **Auxiliary-field** | Decouple interactions with HS transform | DQMC, AFQMC |
| **Continuous-time** | Sample diagrammatic expansion in time | CT-AUX, CT-INT, CT-HYB |

This note sketches each, with the central equations.

## 1. Variational Monte Carlo (VMC)

Choose a parametrised trial wavefunction `|\Psi_T(\alpha)\rangle`.  The
variational energy is

$$
E_T(\alpha) \;=\; \frac{\langle\Psi_T|\hat H|\Psi_T\rangle}{\langle\Psi_T|\Psi_T\rangle}
   \;=\; \frac{\sum_x |\Psi_T(x)|^2\, E_L(x)}{\sum_x |\Psi_T(x)|^2},
\qquad
E_L(x) = \frac{(\hat H\Psi_T)(x)}{\Psi_T(x)}.
$$

`E_L(x)` is the **local energy**.  Sample `x` from `p(x) = |\Psi_T(x)|^2/Z`
via Metropolis; estimate `\langle E_L\rangle`.  Optimise `\alpha` with a stochastic
gradient method.

* If `\Psi_T` is the exact ground state, `E_L(x)` is `x`-independent and
  variance vanishes — the **zero-variance principle** is the most useful
  diagnostic.
* Modern parameterisations: Slater-Jastrow, backflow, neural-network
  wavefunctions (RBM, transformer, CNN — NetKet ecosystem).
* Optimisation: stochastic reconfiguration (Sorella) is essentially
  natural gradient on the Fubini-Study metric.

Limitations: bias of the trial state; the energy is variational, but
expectation values of `\hat A \neq \hat H` are not.

**Companion script:** [`scripts/monte_carlo/vmc_harmonic.jl`](../../scripts/monte_carlo/vmc_harmonic.jl) implements VMC for the 1-D harmonic oscillator and demonstrates the zero-variance principle plus stochastic-gradient optimisation of `\alpha`.

## 2. Diffusion Monte Carlo (DMC)

Project the trial state onto the ground state with imaginary time:

$$
|\Phi_0\rangle \;\propto\; \lim_{\tau\to\infty} e^{-\tau(\hat H - E_T)} |\Psi_T\rangle.
$$

In coordinate basis, `\Psi(\mathbf{R}, \tau)` evolves under a diffusion +
branching equation:

$$
-\partial_\tau \Psi = (\hat H - E_T)\Psi,
\qquad
\hat H = -\tfrac{1}{2}\nabla^2 + V(\mathbf{R}).
$$

The kinetic part is a diffusion; the potential part is a branching/death
process.  Implement with **walkers** `\{\mathbf{R}_i\}`: each step,

```
diffuse:   R_i → R_i + sqrt(δτ) · ξ_i        (Gaussian step)
weight:    w_i → w_i · exp(-δτ · (V(R_i) - E_T))
branch:    duplicate / kill walkers per weight
```

For fermions the wavefunction has nodes — apply the **fixed-node
approximation** to keep `\Psi_T \Psi` non-negative.  The resulting energy
is a variational upper bound *given the nodes of `\Psi_T`*.

DMC is the gold standard for ground-state energies of small molecules and
quantum solids; uncertainty in the nodal surface is the only systematic
left.

## 3. Path-Integral Monte Carlo (PIMC)

Use Trotter's product:

$$
e^{-\beta\hat H} = \bigl(e^{-\Delta\tau\,\hat K} e^{-\Delta\tau\,\hat V}\bigr)^M + O(\Delta\tau^2),
\qquad \Delta\tau = \beta/M.
$$

Insert resolutions of identity between factors.  For distinguishable
particles in a continuous space,

$$
Z \;=\; \int d\mathbf{R}_1\cdots d\mathbf{R}_M\; \prod_{k=1}^M
\,\rho_0(\mathbf{R}_k,\mathbf{R}_{k+1};\Delta\tau)\, e^{-\Delta\tau\,V(\mathbf{R}_k)},
$$

with `\rho_0` the free-particle propagator (a Gaussian).  This is a
classical statistical-mechanics problem in `d+1` dimensions: each quantum
particle becomes a closed *ring polymer* of `M` beads.  Sample with
Metropolis or bisection moves on the beads.

For **bosons**, sum over permutations of imaginary-time end-points
(Ceperley's worm algorithm).  For **fermions**, alternating signs of
permutations cause the sign problem (see §7).

PIMC has been the workhorse for superfluid helium, hydrogen at warm dense
conditions, and many other continuous-space quantum systems.

## 4. Stochastic Series Expansion (SSE)

Expand `Z` as a Taylor series in `\beta`:

$$
Z \;=\; \text{Tr}\sum_{n=0}^{\infty}\frac{(-\beta)^n}{n!}\,\hat H^n
   \;=\; \sum_{n=0}^{\infty}\sum_{S_n}\frac{(-\beta)^n}{n!}\,\langle\alpha|\hat H_{a_n}\cdots\hat H_{a_1}|\alpha\rangle.
$$

Split `\hat H = \sum_b \hat H_b` into local terms.  For sign-problem-free
spin Hamiltonians (e.g., AFM Heisenberg on a bipartite lattice after a
sublattice rotation), every term is non-negative.  Sample the operator
string `S_n` and the boundary state `|\alpha\rangle` simultaneously.

Updates:

* **Diagonal**: insert/remove diagonal operators at random positions
  with the appropriate Metropolis ratio.
* **Loop / operator-loop**: flip a connected set of vertices, an *exact*
  cluster move that exhibits no critical slowing for SU(2) spins.

Implementation references: ALPS, QSWalk, Sandvik's lecture notes.  SSE is
the de facto choice for SU(2)-symmetric spin systems.

## 5. World-line QMC

Closely related to PIMC, but on a lattice.  Discretise imaginary time;
each lattice site has a vertical world-line carrying the local quantum
number through `\tau \in [0,\beta]`.  Off-diagonal operators of `\hat H`
appear as "kinks" or vertices connecting world-lines.  Move kinks via
local or directed-loop updates (Syljuåsen-Sandvik).

Equivalent to SSE in many regimes; choice is largely cultural.

## 6. Auxiliary-field QMC (DQMC, AFQMC)

For lattice fermion problems with quartic interactions
`U\,n_\uparrow n_\downarrow`, the **Hubbard-Stratonovich** transformation
decouples the interaction with an auxiliary field `s`:

$$
e^{-\Delta\tau\, U n_\uparrow n_\downarrow} \;=\; \sum_{s = \pm 1}\,
e^{\,\lambda s\,(n_\uparrow - n_\downarrow)\,+\,\tfrac{\Delta\tau U}{2}(n_\uparrow + n_\downarrow)},
\qquad
\cosh\lambda = e^{\Delta\tau U/2}.
$$

After HS the action is quadratic in fermions, which can be integrated
out:

$$
Z \;=\; \sum_{\{s\}}\, \det\!\bigl[I + \hat B_M(\{s\})\cdots \hat B_1(\{s\})\bigr]_\uparrow\,\det[\cdots]_\downarrow.
$$

Sample `\{s\}` with Metropolis on the determinant ratio.  The half-filled
repulsive Hubbard model has both determinants equal by particle-hole
symmetry — the product is non-negative, no sign problem.  Away from
half-filling, signs appear.

AFQMC for *ab initio* electronic structure (Zhang et al.) uses a
phaseless approximation to control the phase problem and matches
gold-standard CCSD(T) for many molecules.

## 7. The sign problem

For fermions, frustrated spins, or real-time evolution, the weight
`W(\mathcal{C})` can be negative or complex.  Write

$$
\langle A\rangle \;=\; \frac{\langle A\,\sigma\rangle_{|W|}}{\langle\sigma\rangle_{|W|}},\qquad \sigma = W/|W|.
$$

The denominator `\langle\sigma\rangle_{|W|}` shrinks **exponentially** with system
size and inverse temperature:

$$
\langle\sigma\rangle \sim e^{-\beta N \Delta f},
$$

so the relative error of `\langle A\rangle` grows exponentially.  This is **NP-hard**
in general (Troyer-Wiese, 2005) — no universal cure exists.

Mitigations:

* **Sign-problem-free models**: bipartite half-filled Hubbard, AFM
  Heisenberg, Kondo at half-filling.
* **Fixed-node / constrained-path**: project out negative-weight
  contributions at the cost of a controllable bias.
* **Diagrammatic resummation / DiagMC**: rearrange the perturbative
  series so partial cancellations happen before sampling.
* **Reweighting** small-`\beta` results; reformulating with Majorana
  fermions; basis rotations.
* **Phaseless AFQMC** for *ab initio* problems.
* **Lefschetz thimbles / complex Langevin** for finite-density gauge
  theories — still actively researched.

## 8. Continuous-time QMC (CT-QMC)

For impurity problems (quantum impurity in a bath of electrons —
foundational for DMFT) discretise nothing.  Expand `Z` in the
hybridisation `\hat V` between impurity and bath, sample the
*infinite-dimensional* configuration space of pairs of
imaginary times where `\hat V` acts.  Three flavours:

* **CT-INT** — expansion in interaction `U`.
* **CT-AUX** — auxiliary-field decoupling, no `\Delta\tau`.
* **CT-HYB** — expansion in hybridisation; standard for DMFT with
  strong correlation.

CT-QMC plus DMFT (TRIQS, ALPS) is the modern workhorse for correlated
electron systems with local interactions.

## 9. Connections back to classical MC

Every QMC flavour reduces to a *classical* sampling problem:

| QMC | Classical problem | Configuration space |
|---|---|---|
| VMC | sample `|\Psi(x)|^2` | original Hilbert basis |
| PIMC | sample bead positions | `d+1`-D ring polymers |
| SSE | sample operator strings | space-time vertex graphs |
| DQMC/AFQMC | sample HS fields | `L^d \times M` Ising-like fields |
| CT-QMC | sample event sequences | imaginary-time point processes |

This is why all the *technique* of classical MC ([`02_numerical.md`](02_numerical.md))
— autocorrelation, jackknife, parallel tempering, validation —
transfers directly.

## 10. When QMC is the right tool

* Ground state of large lattice models with positive weights (Heisenberg,
  bilinear-biquadratic, dimer models): SSE.
* Bosonic continuum systems: PIMC (helium, ultracold atoms).
* Fermionic small molecules / solids: DMC (with fixed nodes) or AFQMC
  (with phaseless).
* Strongly correlated electrons via DMFT: CT-HYB.
* Critical points in higher-D quantum spin models: SSE with finite-size
  scaling.

QMC is **not** the right tool for:

* Dynamical (real-time) correlators except via analytic continuation,
  which is ill-posed (use MaxEnt / SAC / Padé).
* Frustrated fermion systems away from sign-free points → ED, DMRG,
  tensor networks may be preferred.
* Ground states in 1-D — DMRG is dramatically more efficient.

See [`05_related_methods.md`](05_related_methods.md) for the alternatives.
