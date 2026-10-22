# Companion Scripts — Walkthrough

A guided tour of the eleven Julia scripts in `scripts/monte_carlo/`.  Each
entry ties the script to the underlying equations in the other notes,
explains what the plot should show, and suggests modifications so the
script becomes a *workbench*, not a fixed demo.

The order matches [`04_learning_path.md`](04_learning_path.md): start with
`pi_estimation.jl` and work down.  Every script is self-contained — you
can read just one entry and run it.

---

## 1. `pi_estimation.jl` — hit-or-miss MC and `1/\sqrt N` convergence

**Concept.**  The simplest possible MC algorithm: sample uniformly,
average an indicator.  Demonstrates the universal convergence law of MC
integration.

**Algorithm.**

```
for i = 1 … N:
    x, y ← uniform[0,1]²
    if x² + y² ≤ 1:  hit += 1
estimate π̂ = 4 · hit / N
```

**Why it works.**  The probability of a hit is `\int_0^1\!\int_0^1 \mathbb{1}[x^2+y^2 \le 1]\,dx\,dy = \pi/4`.
By the central limit theorem,

$$
|\hat\pi_N - \pi| \;\sim\; \frac{\sigma}{\sqrt{N}}, \qquad \sigma^2 = \pi(4-\pi)/4.
$$

**Expected output.**

* Panel (a): the unit square, points coloured by inside/outside; the
  estimated π is printed in the title.
* Panel (b): log-log of `|\hat\pi - \pi|` (median over 200 trials) vs `N`;
  slope `= -1/2`.  Compare with the dashed `1/\sqrt N` reference.

**Try modifying.**

* Replace `rand()` with a low-discrepancy sequence (Sobol, Halton) and
  see super-`\sqrt N` convergence — *quasi-Monte Carlo*.
* Use importance sampling with `p(x,y) = (\text{const}) \cdot (1 - x^2 - y^2)`
  inside the circle; the variance should drop.
* Try a higher-dimensional ball; observe that the relative volume in a
  unit hyper-cube shrinks toward zero — the *curse of dimensionality* and
  why importance sampling matters.

**See also.** [`01_theory.md`](01_theory.md) §1-2.

---

## 2. `random_walk.jl` — trajectories and diffusive scaling

**Concept.**  Random walks are the building block of every MCMC chain.
The scaling `\langle x^2(t)\rangle \propto t` is the discrete analogue of
Einstein's diffusion law `\langle x^2\rangle = 2Dt`.

**Why it scales as `t`.**  For independent ±1 steps in 1-D,

$$
\langle x_t^2\rangle = \langle (\textstyle\sum_{i=1}^t \xi_i)^2\rangle = \sum_i\langle\xi_i^2\rangle = t.
$$

**Expected output.**

* Panel (a): six 2-D walkers, each of 4000 steps, coloured by walker.
* Panel (b): `\langle x^2\rangle_t` from 5000 1-D walks, overlaid with `y=t`.
  They should overlap within statistical noise.

**Try modifying.**

* Bias the steps: `\xi \in \{+1, -1\}` with `p \neq 1/2`.  Now
  `\langle x_t\rangle \propto t` (drift) and `\text{Var}(x_t) \propto t`.
* Self-avoiding walk: only accept a move if the new site has not been
  visited.  Numerical exponent `\nu \approx 0.75` in 2-D.
* Lévy flight: draw the step length from a power-law `P(\ell)\sim \ell^{-1-\mu}`
  with `\mu < 2`; `\langle x^2\rangle` becomes super-diffusive.

**See also.** [`01_theory.md`](01_theory.md) §3 (Markov chains).

---

## 3. `ising_metropolis.jl` — the canonical Ising thermodynamics

**Concept.**  Single-spin-flip Metropolis on `H = -J\sum_{\langle ij\rangle}\sigma_i\sigma_j`,
sweeping `T` to locate the critical point.

**Detailed-balance ratio.**

$$
a(\sigma \to \sigma') \;=\; \min\!\left(1,\; e^{-\beta\,\Delta E}\right), \qquad \Delta E = 2 J \sigma_i\!\sum_{j \in \partial i}\!\sigma_j.
$$

Only the five values `\Delta E \in \{-8, -4, 0, 4, 8\}` occur (in units `J=1`),
so `e^{-\beta\Delta E}` is pre-tabulated in a dictionary.

**Expected output.**

* (a) `E/N(T)` — smooth decrease with a kink at `T_c`.
* (b) `\langle|m|\rangle(T)` — order parameter vanishes for `T > T_c`.
* (c) `\chi(T)` — divergent peak at `T_c`; on a finite lattice the peak
  is rounded and shifts slightly above `T_c \approx 2.269`.
* (d) `C_V(T)` — logarithmic peak at `T_c`.
* (e) three snapshots, one ordered, one near-critical, one disordered.

**Try modifying.**

* External field: add `-h\sum_i \sigma_i`.  `\Delta E \to \Delta E + 2 h \sigma_i`.
* Anti-ferromagnetic coupling `J < 0` on the bipartite lattice.  Onsager's
  `T_c` is unchanged; the order parameter is staggered magnetisation.
* `L = 8, 16, 24, 32` → finite-size scaling collapse:
  `\chi_L(T) = L^{\gamma/\nu}\, f(L^{1/\nu}(T-T_c))`.
* Replace Metropolis with **heat-bath** (Gibbs):
  `\sigma_i \leftarrow +1` with prob `\frac{e^{\beta h_i}}{e^{\beta h_i} + e^{-\beta h_i}}`
  where `h_i = J \sum_j \sigma_j`.  Compare autocorrelation.

**See also.** [`01_theory.md`](01_theory.md) §4-5, §7.

---

## 4. `ising_wolff.jl` — cluster updates and critical slowing down

**Concept.**  Local Metropolis suffers near `T_c`: `\tau_{\rm int} \sim L^z`
with `z \approx 2.17`.  Wolff's cluster algorithm flips a Fortuin-Kasteleyn
cluster grown with bond probability

$$
p \;=\; 1 - e^{-2\beta J}.
$$

The dynamical exponent collapses to `z \approx 0.2`.

**Expected output.**

* Left panel: `\tau_{\rm int}(\langle m^2\rangle)` vs `T` for Metropolis vs Wolff on the
  same lattice.  Both curves peak at `T_c`; the Metropolis peak is one to
  two orders of magnitude higher.
* Right panel: speed-up ratio `\tau_{\rm Metro}/\tau_{\rm Wolff}` peaked at `T_c`.

**Try modifying.**

* Sweep `L` and check the dynamical exponent by fitting
  `\tau_{\rm int}(L)` at `T = T_c`.
* Implement Swendsen-Wang (flip *every* cluster, not just one).  Compare.
* Apply Wolff to the 2-D Potts model.  How does `p` change?  (Answer:
  `p = 1 - e^{-\beta J}` since the alignment of an FK bond now requires
  matching one of `q` colours.)

**See also.** [`01_theory.md`](01_theory.md) §6 (cluster algorithms).

---

## 5. `autocorrelation.jl` — `\tau_{\rm int}` and block-binning errors

**Concept.**  Three central tools of MC error analysis:

* normalised autocorrelation `\rho(t) = (\langle x_i x_{i+t}\rangle - \langle x\rangle^2)/\text{Var}(x)`
* integrated time `\tau_{\rm int} = \tfrac12 + \sum_{t \ge 1}\rho(t)` (Madras-Sokal truncation)
* binning analysis: `\sigma_b(\bar x)` flattens once `b \gtrsim 2 \tau_{\rm int}`.

**Naïve error pitfall.**  Reporting `\sigma_x/\sqrt N` claims `\tau_{\rm int} = 0`.
The correct standard error is the plateau of the block-binning curve, or
equivalently `\sigma_x\sqrt{2\tau_{\rm int}/N}`.

**Expected output.**

* Row 1: magnetisation time series at `T = 2.0`, `T_c`, `T = 2.8`.
  At `T_c` the series visibly drifts, with long excursions.
* Row 2: `\rho(t)` decays to zero quickly off-critical; slowly at `T_c`.
  Each panel reports the estimated `\tau_{\rm int}` in its title.
* Row 3: `\sigma(b)` vs `b` on log-log.  The horizontal dashed line is
  the naïve `\sigma/\sqrt N`; the plateau (much higher than the naïve
  line at `T_c`) is the correct standard error.

**Try modifying.**

* Implement **jackknife** for `C_V = \beta^2 N \text{Var}(E/N)` — naive
  error propagation gets this wrong because `C_V` is a non-linear
  functional of the data.
* Plot `2\tau_{\rm int}(L)` vs `L` at `T_c` to extract the dynamical exponent.
* Apply the same diagnostics to the Wolff chain from script 4; the
  binning plateaus should be much closer to the naïve line.

**See also.** [`02_numerical.md`](02_numerical.md) §3-4.

---

## 6. `vmc_harmonic.jl` — variational MC, zero-variance principle, SGD

**Concept.**  Sample `|\Psi_\alpha(x)|^2` by Metropolis and estimate

$$
E(\alpha) = \frac{\sum_x |\Psi_\alpha(x)|^2\,E_L(x)}{\sum_x |\Psi_\alpha(x)|^2}, \qquad
E_L(x) = \frac{(\hat H \Psi_\alpha)(x)}{\Psi_\alpha(x)}.
$$

For the harmonic oscillator with `\Psi_\alpha = e^{-\alpha x^2}`,

$$
E_L(x) = \alpha + \left(\tfrac12 - 2\alpha^2\right)x^2.
$$

At `\alpha = 1/2` the `x^2` coefficient vanishes ⇒ `E_L \equiv 1/2` and
`\text{Var}(E_L) = 0`.  This is the **zero-variance principle**: the more
accurate the trial wavefunction, the lower the statistical noise.

**Stochastic gradient.**  With `O = \partial_\alpha \ln\Psi_\alpha^2 = -2x^2`,

$$
\frac{\partial \langle E\rangle}{\partial \alpha} \;=\; 2\bigl(\langle E_L\, O\rangle - \langle E_L\rangle\langle O\rangle\bigr).
$$

Update `\alpha \leftarrow \alpha - \eta \cdot \widehat{\partial\langle E\rangle/\partial\alpha}`.

**Expected output.**

* (a) `\langle E_L\rangle(\alpha)` with a minimum at `\alpha = 1/2` (value `1/2`).
* (b) `\text{Var}(E_L)(\alpha)` collapsing to zero at `\alpha = 1/2`.
* (c) histogram of samples at `\alpha = 1/2` overlaying `|\Psi_{1/2}|^2 \propto e^{-x^2}`.
* (d) SGD trajectory from `\alpha_0 = 0.30` converging to `\approx 1/2`.

**Try modifying.**

* Add a quartic term `\lambda x^4` to `\hat H` and recompute `E_L`.  Now
  the optimum `\alpha^*` is shifted; find it.
* Implement **stochastic reconfiguration** (Sorella):
  `\delta\alpha = -\eta\, S^{-1}\, \partial E/\partial\alpha` with
  `S = \langle O^2\rangle - \langle O\rangle^2`.  For a single parameter this is just
  rescaling; for multi-parameter trial functions it is the natural
  gradient on the Fubini-Study metric.
* Use a richer trial function `\Psi_T(x) = e^{-\alpha x^2}\,(1 + \gamma x^2)`
  and optimise jointly over `(\alpha, \gamma)`.

**See also.** [`03_quantum_monte_carlo.md`](03_quantum_monte_carlo.md) §1.

---

## 7. `percolation.jl` — geometric phase transition

**Concept.**  Site percolation: occupy each lattice site with probability
`p`.  The largest connected cluster's relative size is an order
parameter.  For 2-D square site percolation,

$$
p_c \;\approx\; 0.59274621\ldots\quad\text{(not exactly known)}
$$

Below `p_c` clusters are small; above `p_c` a spanning cluster occupies
a finite fraction of the system.  Near `p_c`, the cluster-size distribution
follows `P(s) \sim s^{-\tau}` with `\tau \approx 187/91`.

**Expected output.**

* Top panel: largest-cluster fraction vs `p` for `L = 32, 64, 128`.  The
  three curves cross near `p_c` and the transition sharpens with `L`.
* Bottom row: visualisations at `p = 0.5, p_c, 0.7`.  The biggest cluster
  is highlighted in red.

**Try modifying.**

* Switch to **bond percolation**: occupy each lattice *bond* with
  probability `p`.  `p_c^{bond}(\text{square}) = 1/2` exactly (Kesten 1980).
* Cluster size *distribution*: at `p_c`, fit `\log N(s)` vs `\log s` to
  extract `\tau`.
* Different lattices: triangular (site `p_c = 1/2`), honeycomb (site
  `p_c \approx 0.6962`).
* Continuum percolation: discs of random centres, with `p` replaced by
  density `\eta`.

**See also.** [`06_cross_field_examples.md`](06_cross_field_examples.md)
"Statistical geometry and dynamical systems".

---

## 8. `schelling.jl` — micromotives and macrobehaviour

**Concept.**  Two-colour residential grid; each agent moves to an empty
cell if the same-colour fraction among its (≤ 8) Moore neighbours is
below tolerance `\tau`.  For `\tau \gtrsim 0.3`, near-perfect segregation
emerges — Thomas Schelling's 1971 result that *individual mild preferences
produce extreme aggregate sorting*.

**Why it matters.**  The update rule satisfies *no* detailed balance with
respect to any natural energy.  Schelling is a *non-equilibrium* MC and a
canonical example of how MC tooling generalises far beyond physics.

**Expected output.**

* Animated grid (`schelling.mp4`): random initial state evolves to
  large monochrome domains.
* Side panel: a "segregation index" — the mean same-colour fraction
  across occupied cells — climbs from `\approx 0.5` to `\gtrsim 0.85`.

**Try modifying.**

* Vary `\tau \in [0.1, 0.6]`; below some critical tolerance no
  segregation forms.  Plot final segregation vs `\tau` (this is a real
  *phase diagram* of a non-equilibrium model).
* Three or more colours; how do the patterns differ?
* Add tolerance heterogeneity: each agent has a private `\tau_i`.
* Add a cost-of-moving threshold so only agents with same-colour
  fraction *below* `\tau - \epsilon` actually move.  This corresponds to
  Tiebout's "voting with feet" model.

**See also.** [`06_cross_field_examples.md`](06_cross_field_examples.md)
"Social / behavioural sciences".

---

## 9. `wright_fisher.jl` — genetic drift and fixation

**Concept.**  Population of `N` individuals at a biallelic locus.  Each
generation samples `N` offspring from the previous allele pool:

$$
p_{t+1} \;\sim\; \frac{1}{N}\,\text{Binomial}(N, p_t').
$$

With selection, the "effective" frequency before sampling is biased:

$$
p_t' \;=\; \frac{p_t(1+s)}{1 + s\,p_t}.
$$

For neutral evolution (`s = 0`), the **fixation probability** equals the
initial frequency `p_0`.  For weak selection,

$$
P_{\rm fix}(p_0) \;=\; \frac{1 - e^{-2 s\,p_0}}{1 - e^{-2 s}} \qquad(\text{Kimura, 1962}).
$$

**Expected output.**

* (a) thirty neutral trajectories from `p_0 = 0.5`: random meander,
  eventual absorption at 0 or 1.
* (b) thirty selected trajectories from `p_0 = 0.05` with `s = 0.04`:
  most still go extinct, a few sweep to fixation.
* (c) fixation probability vs `p_0`: dots are simulation; the dashed
  diagonal is the neutral theory `P = p_0`; the curve is Kimura with `s=0.02`.

**Try modifying.**

* Mutation: at rate `\mu` per generation, flip alleles.  Stationary
  distribution becomes Beta with parameters depending on `N\mu`.
* Moran process (continuous-time alternative): at each step, one random
  individual dies and a random other reproduces.  Slower, smoother.
* Time-varying `N(t)`: bottleneck (sudden drop in `N`) accelerates drift.
* Two loci with linkage: track haplotype frequencies and compute linkage
  disequilibrium `D = p_{AB} - p_A p_B`.

**See also.** [`06_cross_field_examples.md`](06_cross_field_examples.md)
"Biology and evolution".

---

## 10. `forest_fire.jl` — self-organised criticality

**Concept.**  Drossel-Schwabl cellular automaton with three states
(empty, tree, burning) and two parameters: tree growth `g` and
lightning `f`.  In the **slow-driving limit** `f/g \to 0`, the model
self-organises to a critical state with power-law cluster sizes,

$$
N(s) \;\propto\; s^{-\tau}, \qquad \tau \approx 2 \;\text{in 2-D.}
$$

This is the cleanest *self-organised critical* (SOC) model: no parameter
tuning, criticality emerges from the dynamics.

**Expected output.**

* Animation (`forest_fire.mp4`): black (empty), green (tree), red (burning).
  Fires sweep through forested regions, regrowth fills back in.
* Side panel: log-log histogram of tree-cluster sizes, with a `s^{-2}`
  reference line.  The fit is approximate at this `L`; the true SOC
  scaling appears at much larger lattices.

**Try modifying.**

* Vary `f` and `g` over orders of magnitude.  Confirm critical behaviour
  emerges as `f/g \to 0`.
* Add wind: bias the burning-spread direction.  How does the cluster
  shape change?
* Sandpile (Bak-Tang-Wiesenfeld) is an even cleaner SOC model — try
  implementing it for comparison.

**See also.** [`06_cross_field_examples.md`](06_cross_field_examples.md)
"Statistical geometry and dynamical systems".

---

## 11. `gillespie.jl` — exact stochastic simulation algorithm

**Concept.**  For a reaction network with propensities `a_i(\mathbf{n})`,
the exact next-event distribution is

$$
\Delta t \sim \text{Exp}\!\left(\textstyle\sum_i a_i\right), \qquad
\Pr(\text{reaction }k\mid\text{event}) = \frac{a_k}{\sum_i a_i}.
$$

Demonstrated on stochastic Lotka-Volterra:

```
R ──────►  R + R     (rate α R)        prey birth
R + F ──►  F + F     (rate β R F)      predation + predator birth
F ──────►  ∅         (rate γ F)        predator death
```

Compared with the deterministic ODE, the stochastic version can drive
either species to extinction; oscillations have random period and
amplitude.

**Expected output.**

* (a) populations vs `t` — characteristic oscillations.
* (b) phase portrait `F(R)` — a noisy closed orbit drifting outward.
* (c) ensemble of thirty independent runs showing population variability
  and occasional extinction.

**Try modifying.**

* `\tau`-leaping: advance by a small `\Delta\tau` and draw Poisson counts
  for each reaction.  Compare runtime and accuracy.
* Add migration: small influx of rabbits at constant rate.  Eliminates
  extinction; introduces a stationary distribution.
* Replace LV with the SIR epidemic
  `S+I \to 2I` (rate `\beta SI/N`),
  `I \to R` (rate `\gamma I`).
  Compute the final epidemic size distribution and the basic
  reproduction number `R_0 = \beta/\gamma`.
* Spatial Gillespie: place agents on a grid; reactions are local.  Compare
  the well-mixed and spatial cases.

**See also.** [`06_cross_field_examples.md`](06_cross_field_examples.md)
"Chemistry and reaction kinetics", "Biology and evolution".

---

## Suggested progression

If using these scripts as a course, the natural order is:

1. `pi_estimation.jl`  →  the central convergence law
2. `random_walk.jl`    →  trajectories and diffusion
3. `ising_metropolis.jl` →  the canonical lattice model
4. `autocorrelation.jl`  →  errors done right
5. `ising_wolff.jl`    →  cluster acceleration
6. `vmc_harmonic.jl`   →  bridge to quantum
7. `percolation.jl`    →  geometric criticality
8. `gillespie.jl`      →  master equations and exact stochastic dynamics
9. `wright_fisher.jl`  →  noisy population dynamics
10. `schelling.jl`     →  non-equilibrium social dynamics
11. `forest_fire.jl`   →  self-organised criticality

Each step adds one new concept on top of the previous.  Once 1-6 are
internalised, items 7-11 can be read in any order — they show that the
*same* sampling machinery works far outside physics.

## How to extend

The scripts deliberately favour clarity over performance.  When you
modify one and find yourself wanting:

* a faster RNG → switch to `Xoshiro256++` (the Julia default);
* parallel chains → use `Threads.@threads` on the outer loop;
* hot inner loops → add `@inbounds`, replace `Dict` lookups with arrays;
* larger systems → reduce `Float64 → Int8` for spin storage, bit-pack;
* automatic differentiation of `E_L` in VMC → use `ForwardDiff.jl`;
* full QMC framework → look at `Carlo.jl`, `NetKet`, or `ALPS`.

For more demanding projects (DMRG, AFQMC, CT-HYB), see
[`05_related_methods.md`](05_related_methods.md).
