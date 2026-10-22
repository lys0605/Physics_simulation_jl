# Classical Monte Carlo — Theory and Principles

## 1. The problem

In equilibrium statistical mechanics the canonical expectation of an
observable `A` at inverse temperature `β = 1/(k_B T)` is

$$
\langle A\rangle \;=\; \frac{1}{Z}\sum_{x\in\mathcal{X}} A(x)\, e^{-\beta H(x)}, \qquad
Z = \sum_{x\in\mathcal{X}} e^{-\beta H(x)}.
$$

For a lattice of `N` spin-½ degrees of freedom the configuration space has
size `|\mathcal{X}| = 2^N`.  Even for `N = 100`, `2^N \approx 10^{30}` —
direct enumeration is hopeless.

**Idea:** approximate the average by sampling `M` configurations `\{x_i\}`
drawn from the **Boltzmann distribution** `\pi(x) = e^{-\beta H(x)}/Z`, then

$$
\langle A\rangle \;\approx\; \frac{1}{M}\sum_{i=1}^{M} A(x_i).
$$

The central limit theorem gives an error `\sigma/\sqrt{M_{\text{eff}}}` where
`\sigma^2 = \langle A^2\rangle - \langle A\rangle^2` and `M_{\text{eff}}` is the number
of *independent* samples.  Crucially, the convergence rate `1/\sqrt{M}`
does **not** depend on the dimension of `x` — this is what makes MC
viable in 10²–10⁹ degrees of freedom.

## 2. Importance sampling

Naive uniform sampling fails because almost all of `\mathcal{X}` has
`e^{-\beta H} \approx 0`.  Importance sampling instead generates `x_i`
with frequency proportional to `\pi(x)`.

But we do not know `Z`!  The trick: build a stochastic process whose
**stationary distribution** is `\pi(x)` — without ever computing `Z`.
A Markov chain whose ratio `\pi(x')/\pi(x) = e^{-\beta\Delta H}` only
involves an *energy difference* sidesteps the partition function entirely.

## 3. Markov chains

A Markov chain on `\mathcal{X}` is defined by transition probabilities
`P(x\to x')` satisfying `\sum_{x'} P(x\to x') = 1`.  A distribution `\pi`
is **stationary** if

$$
\sum_x \pi(x) P(x\to x') = \pi(x') \qquad \text{(stationarity).}
$$

For a chain to converge to `\pi`, two conditions suffice:

1. **Ergodicity / irreducibility.**  Every configuration is reachable
   from every other in a finite number of steps.
2. **Aperiodicity.**  The chain does not return to a state only at fixed
   integer multiples of some period `T > 1`.

A *sufficient* (and convenient) replacement for stationarity is the
**detailed balance** condition

$$
\boxed{\;\pi(x)\,P(x\to x') \;=\; \pi(x')\,P(x'\to x)\;}
$$

(Detailed balance ⇒ stationarity, but not the reverse.)  Almost every MC
algorithm in physics is built to satisfy detailed balance because it is
local and easy to verify.

## 4. Metropolis-Hastings

Decompose `P(x\to x') = q(x\to x')\,a(x\to x')` into a **proposal**
distribution `q` and an **acceptance** probability `a`.  Detailed balance
fixes the only choice consistent with arbitrary `q`:

$$
a(x\to x') \;=\; \min\!\left(1,\; \frac{\pi(x')\, q(x'\to x)}{\pi(x)\, q(x\to x')}\right).
$$

For symmetric proposals `q(x\to x') = q(x'\to x)`:

$$
a(x\to x') \;=\; \min\!\bigl(1,\; e^{-\beta\,[H(x') - H(x)]}\bigr).
$$

This is the **Metropolis criterion** — uphill moves in energy are accepted
with the Boltzmann factor `e^{-\beta\Delta H}` and downhill moves always.

**Companion script:** [`scripts/monte_carlo/ising_metropolis.jl`](../../scripts/monte_carlo/ising_metropolis.jl) implements this acceptance rule for the 2-D Ising model.

### Single-spin-flip Metropolis (Ising example)

```
for sweep = 1 … N_sweeps:
    for each site i (random or sequential):
        ΔE = energy change from flipping σ_i
        if ΔE ≤ 0 or rand() < exp(-β ΔE):
            σ_i ← −σ_i
            E  += ΔE
    record observables
```

One *sweep* = `N_{\text{sites}}` proposed flips so that on average every
site is touched once.

## 5. Gibbs / heat-bath sampling

A different update sets a single variable to its **conditional**
distribution given the others:

$$
\sigma_i \,\sim\, p(\sigma_i \mid \{\sigma_{j\neq i}\}) \;=\; \frac{e^{-\beta H(\sigma_i, \cdots)}}{\sum_{\sigma_i'}\, e^{-\beta H(\sigma_i', \cdots)}}.
$$

For Ising this reduces to a single sigmoid; for Potts models or continuous
spins it remains an exact local update.  No rejection — but the
configuration must permit the conditional to be sampled cheaply.

## 6. Cluster algorithms

Single-spin-flip Metropolis suffers near a phase transition because
correlated regions of size `\xi` flip incoherently.  Critical slowing
down: autocorrelation time scales as `\tau \sim \xi^z` with dynamical
exponent `z \approx 2` for local algorithms.

**Swendsen-Wang** (Potts/Ising): build clusters via Fortuin-Kasteleyn
bond percolation,

$$
p_{\text{bond}} = 1 - e^{-2\beta J}\quad\text{(if spins agree, else } p_{\text{bond}} = 0\text{)},
$$

then flip each cluster with probability ½.  Detailed balance follows from
the FK joint measure on (spins, bonds).

**Wolff:** like Swendsen-Wang but grows a *single* cluster from a random
seed and flips it deterministically.  Dynamical exponent collapses to
`z \approx 0.2` (Ising d=2).

Cluster moves are the right idea wherever the model has a percolation
representation: Ising, Potts, O(N) sigma models, quantum spins via SSE
loop updates (see [`03_quantum_monte_carlo.md`](03_quantum_monte_carlo.md)).

**Companion script:** [`scripts/monte_carlo/ising_wolff.jl`](../../scripts/monte_carlo/ising_wolff.jl) implements Wolff cluster updates and compares `\tau_{\rm int}` against single-spin Metropolis near `T_c`.

## 7. Observables and fluctuation relations

Several observables and their fluctuations follow directly from `Z`:

$$
E = -\partial_\beta \ln Z, \qquad
C_V = \frac{\partial \langle E\rangle}{\partial T} = \frac{\langle E^2\rangle - \langle E\rangle^2}{k_B T^2}.
$$

For magnetic systems with magnetisation `M`:

$$
\chi \;=\; \frac{\beta}{N}\!\left(\langle M^2\rangle - \langle M\rangle^2\right), \qquad
U_4 \;=\; 1 - \frac{\langle M^4\rangle}{3\,\langle M^2\rangle^2}\quad(\text{Binder cumulant}).
$$

The Binder cumulant `U_4(T,L)` is `L`-independent at the critical point —
crossing the curves for different system sizes pins down `T_c`.

## 8. Ergodicity in practice: equilibration, autocorrelation

The Markov chain converges geometrically: `\|P^t(x_0,\cdot) - \pi\|_{\rm TV} \sim e^{-t/\tau_{\rm mix}}`.
Two operationally important times:

* **Equilibration / burn-in** `\tau_{\rm eq}` — number of sweeps to forget
  the initial condition.  Estimate by running from two extreme
  configurations (`T=0` ordered, `T=\infty` random) and waiting for
  observables to agree.
* **Integrated autocorrelation time** `\tau_{\rm int}` — controls the
  variance of the time-averaged estimator,

  $$
  \tau_{\rm int}^A \;=\; \frac{1}{2}\,+\,\sum_{t=1}^{\infty}\,\frac{C^A(t)}{C^A(0)},\qquad
  C^A(t) = \langle A_{i+t}A_i\rangle - \langle A\rangle^2.
  $$

  The effective number of independent samples is `M_{\rm eff} = M/(2\tau_{\rm int})`.

## 9. Sources of error

| Error | Cause | Detection |
|---|---|---|
| Statistical | finite sample variance | binning / jackknife |
| Equilibration bias | starting too far from equilibrium | run two initial conditions |
| Ergodicity failure | metastable basins (1st-order transitions, glassy systems) | replica exchange tests |
| Critical slowing | local updates near `T_c` | autocorrelation explodes |
| Auto-correlated errors | naive `\sigma/\sqrt{N}` | always multiply by `\sqrt{2\tau_{\rm int}}` |

The pragmatic side of these issues is the subject of
[`02_numerical.md`](02_numerical.md).

**Companion scripts:**
[`pi_estimation.jl`](../../scripts/monte_carlo/pi_estimation.jl) — `1/\sqrt N` convergence;
[`random_walk.jl`](../../scripts/monte_carlo/random_walk.jl) — Markov-chain trajectories.

## 10. Beyond canonical sampling

Sometimes the canonical ensemble is the wrong target:

* **Multi-canonical / Wang-Landau** sample the density of states `g(E)`
  directly, giving access to all temperatures from one run.
* **Replica exchange (parallel tempering)** swaps configurations between
  chains at different `β` to traverse free-energy barriers.
* **Umbrella sampling** with reweighting estimates free-energy profiles
  along reaction coordinates.

These are introduced in [`02_numerical.md`](02_numerical.md).

## 11. Summary in one paragraph

Classical Monte Carlo works because any *ergodic* Markov chain whose
*detailed-balance* equation enforces the Boltzmann weight will, in the
long run, sample configurations with frequency `\pi(x) \propto e^{-\beta H(x)}`.
Estimators of equilibrium observables converge as `1/\sqrt{M_{\rm eff}}`,
*independent* of configuration-space dimension.  The art is (i) designing
local updates that satisfy detailed balance, (ii) accelerating them with
cluster or non-local moves when correlations grow, and (iii) honestly
estimating the autocorrelation time so error bars are correct.
