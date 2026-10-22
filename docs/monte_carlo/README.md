# Monte Carlo Methods for Physics

A self-contained set of notes on Monte Carlo (MC) and Quantum Monte Carlo (QMC)
methods, oriented toward someone learning to *simulate* physics with these
tools.  The five companion notes are layered:

| File | Scope |
|---|---|
| [`01_theory.md`](01_theory.md) | Classical statistical MC: importance sampling, Markov chains, Metropolis-Hastings, detailed balance, cluster algorithms, critical slowing down |
| [`02_numerical.md`](02_numerical.md) | Pragmatic numerical work: RNGs, equilibration, autocorrelation, binning/jackknife/bootstrap errors, parallel tempering, Wang-Landau, profiling |
| [`03_quantum_monte_carlo.md`](03_quantum_monte_carlo.md) | Quantum MC families: VMC, DMC, PIMC, SSE, world-line, auxiliary-field, CT-QMC; the sign problem and mitigations |
| [`04_learning_path.md`](04_learning_path.md) | Step-by-step learning outline from MC pi-estimation to AFQMC, with milestones, exercises, and reading per stage |
| [`05_related_methods.md`](05_related_methods.md) | Brief survey of complementary methods (ED, DMRG/MPS, PEPS, DMFT, DFT, HF/CC, fRG) for future expansion |
| [`06_cross_field_examples.md`](06_cross_field_examples.md) | Famous MC models outside physics — sociology, biology, economics, chemistry, ML — with central equations and links to companion scripts |
| [`07_scripts_walkthrough.md`](07_scripts_walkthrough.md) | Per-script guide: what each Julia file demonstrates, key equations, expected output, exercises, and pointers back into notes 01-06 |

## Companion scripts

Runnable Julia demos live in [`scripts/monte_carlo/`](../../scripts/monte_carlo/),
indexed in the top-level [`README.md`](../../README.md).  Each is short and
self-contained:

| Script | Concept |
|---|---|
| `pi_estimation.jl` | hit-or-miss MC, `1/\sqrt N` convergence |
| `random_walk.jl`   | trajectories + `\langle x^2\rangle \sim t` |
| `ising_metropolis.jl` | the canonical Ising thermodynamics sweep |
| `ising_wolff.jl`   | cluster algorithm vs critical slowdown |
| `autocorrelation.jl` | `\tau_{\rm int}`, block-binning error plateau |
| `vmc_harmonic.jl`  | VMC, zero-variance principle, gradient descent |
| `percolation.jl`   | geometric phase transition at `p_c \approx 0.5927` |
| `schelling.jl`     | residential segregation (sociology, animated) |
| `wright_fisher.jl` | genetic drift, fixation probability (biology) |
| `forest_fire.jl`   | self-organised criticality (ecology, animated) |
| `gillespie.jl`     | exact SSA on stochastic Lotka-Volterra (chemistry/biology) |

## How to use

* If you want the **principle and equations** of a method, start with
  [`01_theory.md`](01_theory.md) and [`03_quantum_monte_carlo.md`](03_quantum_monte_carlo.md).
* If you want to **write and debug working simulations**, read
  [`02_numerical.md`](02_numerical.md) alongside any algorithm note.
* If you are deciding **what to learn next**, see
  [`04_learning_path.md`](04_learning_path.md) for the recommended order.
* If you wonder **what other methods exist** and when to reach for them,
  see [`05_related_methods.md`](05_related_methods.md).

## Why Monte Carlo

In statistical and quantum many-body physics one almost always needs to
evaluate sums or integrals of the form

$$
\langle A\rangle \;=\; \frac{1}{Z}\sum_{x\in\mathcal{X}}\,A(x)\,w(x), \qquad
Z = \sum_{x\in\mathcal{X}} w(x),
$$

where the configuration space `\mathcal{X}` is exponentially large.
Exact enumeration is impossible past `~30` sites or `~10` particles.
Importance-sampled MC reduces the work to a *constant prefactor* times
the autocorrelation time of the chain, independent of dimension — at the
cost of an `O(1/\sqrt{N})` statistical error.  This trade-off is the
single reason MC dominates so much of computational many-body physics.
