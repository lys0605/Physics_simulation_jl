# A Learning Path for Monte Carlo Simulation in Physics

A step-by-step progression from "compute π by random sampling" to AFQMC
for correlated electrons.  Each stage lists **(a) what to build**,
**(b) what to verify**, and **(c) what to read next**.  Plan one to two
weeks per stage if working part-time.

## Stage 0 — Prerequisites

You should already be comfortable with:

* Linear algebra: eigendecomposition, sparse matrices.
* Statistical mechanics: partition function, canonical ensemble, free
  energy, critical exponents.
* Programming: a language with fast loops (Julia, Python+Numba, C++,
  Fortran), basic plotting.
* Probability: variance, central limit theorem, expectation.

If any of these are weak, treat Stage 1 as a chance to firm them up
before moving on.

## Stage 1 — Sampling basics  *(one week)*

### Build

1. **Estimate π** by sampling the unit square and counting how many points
   land in the inscribed circle.  Plot the error vs `\sqrt{N}` to see
   `O(1/\sqrt{N})` convergence.
2. **Sample a Gaussian** by (a) Box-Muller and (b) Metropolis with a
   uniform proposal.  Compare histograms.
3. **Compute** `\int_0^1 \sin^2(\pi x)\,dx` with uniform sampling and
   importance sampling using `p(x) = \sin(\pi x)`.  Show variance reduction.
4. **1-D random walk** — verify diffusive scaling `\langle x^2\rangle \sim t`.

### Verify

* Reproduce the analytical answer for each integral.
* Plot statistical errors vs `1/\sqrt{N}` (loglog) — slope should be `-0.5`.

### Read

* MacKay, *Information Theory, Inference, and Learning Algorithms*, Ch. 29-30.
* Newman & Barkema, *Monte Carlo Methods in Statistical Physics*, Ch. 1-3.

## Stage 2 — The 2-D Ising model  *(one to two weeks)*

The Ising model is the hydrogen atom of computational statistical
mechanics — everything generalises from it.

### Build

1. **Metropolis Ising** on an `L×L` square lattice with periodic boundary
   conditions, `H = -J\sum_{\langle ij\rangle}\sigma_i\sigma_j`.
2. Implement observables: `E`, `M`, `M^2`, `M^4`.
3. Run a temperature sweep from `T = 1.0` to `T = 4.0` at `L = 8, 16, 32`.
4. Implement **autocorrelation analysis** and **jackknife error bars**.

### Verify

* `T_c \approx 2.269\,J/k_B` (Onsager).  Locate by:
  * peak of `\chi(T)` and `C_V(T)`.
  * crossing of Binder cumulants `U_4(T, L)` for different `L`.
* Critical exponents: `\beta = 1/8, \gamma = 7/4, \nu = 1` via finite-size
  scaling.

### Extend

* Add **Wolff cluster** updates; compare autocorrelation time vs Metropolis
  near `T_c`.  Critical slowing should drop from `\tau \sim L^2` to `\tau \sim L^{0.2}`.
* Try the **Wang-Landau** algorithm to compute `g(E)` directly.

### Read

* Sandvik lecture notes, *Computational Studies of Quantum Spin Systems*,
  §3-4 (the Ising part).
* Wolff, Phys. Rev. Lett. **62**, 361 (1989).
* Ferrenberg & Swendsen, Phys. Rev. Lett. **63**, 1195 (1989) — reweighting.

## Stage 3 — Beyond Ising  *(one to two weeks)*

Diversify the model zoo.

### Build (pick two)

* **Potts model** (`q`-state) — first-order transition for `q > 4` in 2-D.
* **XY model** — continuous spins, Kosterlitz-Thouless transition.
* **Hard-sphere fluid** — event-driven MC or NPT sampling, equation of state.
* **Lennard-Jones liquid** — Wales-Doye basin-hopping for ground states.

### Verify

* For XY: locate `T_{\rm KT}` via the helicity modulus jump
  `\rho_s(T_{\rm KT}^-) = (2/\pi) T_{\rm KT}`.
* For Potts: check latent heat across the transition.

### Read

* Krauth, *Statistical Mechanics: Algorithms and Computations*.
* Frenkel & Smit, *Understanding Molecular Simulation* (for continuum).

## Stage 4 — Variational Monte Carlo  *(two weeks)*

Bridge from classical to quantum.

### Build

1. **Harmonic oscillator VMC**: trial `\Psi_T(x) = e^{-\alpha x^2}`.  Sample
   `x \sim |\Psi_T|^2`, compute `\langle E_L(x)\rangle`, optimise `\alpha`.
2. **Helium atom VMC**: `\Psi_T = e^{-Z(r_1+r_2)}\,e^{r_{12}/(2(1+\alpha r_{12}))}`
   (Slater-Jastrow).  Aim for `E \approx -2.890` Ha (exact: `-2.903`).
3. **1-D Heisenberg chain** with Jastrow factor; compare to exact ED for
   `L = 10`.

### Verify

* **Zero-variance principle**: as `\alpha` approaches the optimum, the
  variance of `E_L` shrinks dramatically.  Plot it.

### Read

* Foulkes et al., Rev. Mod. Phys. **73**, 33 (2001) — *Quantum Monte
  Carlo simulations of solids*.
* Becca & Sorella, *Quantum Monte Carlo Approaches for Correlated
  Systems*, Ch. 3-4.

## Stage 5 — Imaginary-time projection  *(two weeks)*

### Build

1. **Diffusion Monte Carlo** for the 1-D harmonic oscillator (no nodes).
   Confirm `E_0 = \hbar\omega/2`.
2. **DMC for `H_2`**: fixed-node with single-determinant nodes.
3. **PIMC for a single quantum particle** in a 1-D well at `T > 0`.  Plot
   the imaginary-time correlation `\langle x(\tau)x(0)\rangle`.

### Verify

* `\Delta\tau \to 0` extrapolation for DMC energies.
* PIMC: `M \to \infty` extrapolation of Trotter error.

### Read

* Foulkes et al., as above.
* Ceperley, Rev. Mod. Phys. **67**, 279 (1995) — *Path integrals in the
  theory of condensed helium*.

## Stage 6 — Sign-free lattice QMC  *(two to three weeks)*

### Build

1. **SSE for the 1-D AFM Heisenberg chain**, with diagonal + operator-loop
   updates.  Compare ground-state energy with Bethe ansatz: `E_0/N \to -\ln 2 + 1/4 \approx -0.443` (units `J=1`).
2. **SSE for the 2-D AFM Heisenberg model** on `L = 8, 16, 32`.
3. **Stochastic Series Expansion thermalisation** at multiple `\beta`,
   extract `\chi(T)`, magnetic susceptibility, structure factor.

### Verify

* 2-D Néel order: spontaneous staggered magnetisation `m_s \approx 0.307`
  at `T = 0` for nearest-neighbour Heisenberg.

### Read

* Sandvik, *Stochastic series expansion methods*, arXiv:1909.10591.
* Syljuåsen & Sandvik, Phys. Rev. E **66**, 046701 (2002).

## Stage 7 — Sign problem and constrained methods  *(two to three weeks)*

### Build

1. **DQMC for half-filled square-lattice Hubbard** (`U=4t`, `L=8`).
   Sign-free: every Metropolis ratio is positive.
2. Move *away from half filling* and watch `\langle\sigma\rangle` decay.  Plot
   `-\ln\langle\sigma\rangle` vs `\beta N`.
3. **Constrained-path AFQMC** with a trial wavefunction; reproduce the
   half-filled energy and try a doped point.

### Verify

* Half-filled `U=4t` square-lattice ground-state energy: `E_0/N \approx -0.860 \pm 0.005` (extrapolated).
* Antiferromagnetic structure factor `S(\pi,\pi)` grows with `L`.

### Read

* Assaad & Evertz, *Computational Many-Particle Physics*, Lect. Notes Phys. 739, Ch. 10.
* Zhang & Krakauer, Phys. Rev. Lett. **90**, 136401 (2003) — phaseless AFQMC.

## Stage 8 — Continuous-time and embedding methods  *(three weeks)*

### Build

1. **CT-HYB single-impurity Anderson model** (one orbital, semicircular
   bath).  Extract spectral function via MaxEnt.
2. **DMFT loop**: alternate CT-HYB impurity solver with self-consistency
   on the bath.  Reproduce the Mott transition in the half-filled Hubbard
   model on the Bethe lattice.

### Verify

* Mott transition at `U_c \approx 5\!-\!6\,t` in single-site DMFT.

### Read

* Gull et al., Rev. Mod. Phys. **83**, 349 (2011) — *Continuous-time Monte
  Carlo methods for quantum impurity models*.
* Georges et al., Rev. Mod. Phys. **68**, 13 (1996) — DMFT.

## Stage 9 — Frontiers  *(open-ended)*

By this point you can read the current QMC literature.  Possible directions:

* **Neural-network quantum states** (Carleo & Troyer) — VMC with RBM /
  transformer / CNN ansatz, plus stochastic reconfiguration.
* **Diagrammatic Monte Carlo** for unbiased fermionic field theory.
* **Tensor-network Monte Carlo hybrids**.
* **Quantum-circuit Monte Carlo** for simulating noisy quantum devices.
* **Real-time / Keldysh** QMC and dynamical sign-problem mitigation.
* **Gauge-theory QMC** with Lefschetz thimbles or complex Langevin.

## Milestones to know you have arrived

You are ready to do original work when, on a model you have not seen before:

1. you can derive whether MC sign-free updates exist;
2. you can choose between VMC / DMC / SSE / DQMC / CT-HYB and justify it;
3. you can budget the simulation (system size, sweeps, parallelisation);
4. you can quantify *every* systematic and statistical error you report.

## Companion code

The course of study above benefits from a working Metropolis Ising as a
"hello world".  Once that exists, every subsequent stage modifies one
piece at a time: lattice → continuous spins → quantum particle in
imaginary time → operator string → auxiliary field.  Code reuse is
substantial; the conceptual leaps are in the *measure*, not the loop.

The Julia scripts in `scripts/monte_carlo/` already cover Stages 1, 2, and
the first half of Stage 4 — see [`07_scripts_walkthrough.md`](07_scripts_walkthrough.md)
for a per-script guide with exercises.  The remaining stages (SSE, DQMC,
CT-HYB) are not yet in the repo; treat them as "build it yourself" using
the equations in [`03_quantum_monte_carlo.md`](03_quantum_monte_carlo.md).
