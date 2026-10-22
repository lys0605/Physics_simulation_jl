# Monte Carlo — Numerical and Implementation Notes

A working MC code is half algorithm, half book-keeping.  This note covers
the book-keeping: random number generation, equilibration tests,
autocorrelation, error bars, accelerated samplers, and validation.

## 1. Random number generators

A modern simulation needs a generator that is:

* **Long period.**  At least `2^{64}`; ideally `2^{128}` or more.
* **Statistically clean.**  Passes BigCrush.
* **Fast.**  Inner loops can generate `10^9` deviates.
* **Splittable / streamable.**  So parallel chains do not share output.

Acceptable choices today:

| Generator | Period | Notes |
|---|---|---|
| Xoshiro256++ / Xoshiro256** | `2^{256}-1` | Julia default since 1.7 |
| PCG64 (DXSM) | `2^{128}` | numpy default since 1.17 |
| Mersenne Twister MT19937 | `2^{19937}-1` | OK but old, fails some tests |
| Philox4×64 | counter-based | trivially splittable across MPI ranks |

**Bad choices** still found in code: `rand()` from libc, linear congruential
generators with `< 2^{32}` period, `srand(time(NULL))` for seeding.

**Seeding.**  For reproducibility, log the seed.  For parallel chains, use
`Random.seed!(state)` with distinct seeds *per chain* (counter-based RNGs
make this trivial; otherwise use `Random.Xoshiro` instances or
`StableRNGs.jl`).

## 2. Equilibration / thermalisation

Discard the first `\tau_{\rm eq}` sweeps to remove dependence on the initial
condition.  Two complementary checks:

1. **Hot vs cold start.**  Run twice — one chain from a random
   configuration (`T=\infty`), one from an ordered (`T=0`) configuration.
   Equilibration is reached when their running averages overlap within
   the error.
2. **Block averages.**  Split the chain into `k` consecutive blocks.  If
   the first block differs systematically from the rest, increase
   `\tau_{\rm eq}`.

Conservative rule of thumb: discard `5\,\tau_{\rm int}` (after estimating
`\tau_{\rm int}` from a long pilot run).

## 3. Autocorrelation time

Define the time autocorrelation of observable `A` at lag `t`:

$$
\rho_A(t) \;=\; \frac{\langle A_i A_{i+t}\rangle - \langle A\rangle^2}{\langle A^2\rangle - \langle A\rangle^2}.
$$

The **integrated autocorrelation time** is

$$
\tau_{\rm int}^A \;=\; \tfrac{1}{2} \,+\, \sum_{t=1}^{\infty} \rho_A(t).
$$

Practical estimator (Madras-Sokal): truncate the sum at the first `t^*`
where `t^* \ge c\,\tau_{\rm int}(t^*)` with `c \in [4, 10]`.  Otherwise
the estimator's variance diverges.

Effective number of independent samples:

$$
M_{\rm eff} \;=\; \frac{M}{2\,\tau_{\rm int}}.
$$

**Anything not multiplied by `\sqrt{2\,\tau_{\rm int}}` will under-estimate
errors.**  Common bug: report `\sigma/\sqrt{M}` directly — you are claiming
`\tau_{\rm int} = 0`.

**Companion script:** [`scripts/monte_carlo/autocorrelation.jl`](../../scripts/monte_carlo/autocorrelation.jl) plots `\rho(t)` and the block-binning curve at three temperatures, showing the gap between the naive and corrected errors. The `tau_int` and `block_curve` helpers used there live in [`src/monte_carlo/mc_utils.jl`](../../src/monte_carlo/mc_utils.jl).

## 4. Error estimation

Three robust methods, in order of complexity:

### Binning analysis

Group consecutive samples into blocks of size `b`, compute the block
means, and estimate the variance of the block means:

$$
\sigma_b^2(\bar A) \;=\; \frac{1}{N_b(N_b-1)}\sum_{k=1}^{N_b}\bigl(\bar A_k - \bar A\bigr)^2.
$$

Increase `b` until the estimate plateaus — that plateau value is the
true variance.  The corresponding `b^*` satisfies `2\tau_{\rm int} \approx b^*`.

### Jackknife

Remove one block at a time, recompute the estimator on the remaining
`N_b - 1` blocks, get `\hat\theta_{(k)}`.  Then

$$
\hat\theta_{\rm jack} = \frac{1}{N_b}\sum_k \hat\theta_{(k)},\qquad
\sigma_{\rm jack}^2 = \frac{N_b-1}{N_b}\sum_k (\hat\theta_{(k)} - \hat\theta_{\rm jack})^2.
$$

Works for **derived quantities** (`C_V = \beta^2 \mathrm{Var}(E)`,
Binder cumulants, ratios) where naive error propagation is wrong.

### Bootstrap

Resample blocks with replacement `R` times, compute `\hat\theta^{(r)}` on
each, report the standard deviation across `r=1\ldots R`.  Slightly more
flexible than jackknife (handles bias correction, confidence intervals),
slower (`R \sim 10^3` typical).

## 5. Histogram reweighting

A single run at `\beta_0` carries information about nearby temperatures.
Multi-histogram reweighting (Ferrenberg-Swendsen) combines runs at several
`\beta_k`:

$$
\langle A\rangle_\beta \;=\; \frac{\sum_k\,\sum_i\, w_i^{(k)}(\beta)\, A_i^{(k)}}{\sum_k\,\sum_i\, w_i^{(k)}(\beta)},
\qquad
w_i^{(k)}(\beta) = \frac{e^{-\beta E_i^{(k)}}}{\sum_j n_j\, e^{f_j - \beta_j E_i^{(k)}}},
$$

with free energies `f_j` solved self-consistently.  Used heavily in
finite-size scaling near `T_c`.

## 6. Parallel tempering (replica exchange)

Run `R` chains at temperatures `T_1 < T_2 < \cdots < T_R`.  Periodically
attempt to swap configurations between neighbouring chains with
probability

$$
P_{\rm swap}(i, i+1) \;=\; \min\!\left(1,\, e^{(\beta_i - \beta_{i+1})(E_i - E_{i+1})}\right).
$$

The hot chains traverse free-energy barriers, the cold chains accumulate
statistics at the temperature of interest.  Tune `\{T_k\}` so the average
swap acceptance is `\sim 20\!-\!40\%` between every adjacent pair.

## 7. Wang-Landau sampling

Estimate the density of states `g(E)` directly by performing a random
walk in energy with probability `\propto 1/g(E)` and refining `g(E)` on
the fly:

```
init g(E) ← 1 for all E,  histogram H(E) ← 0
log-modification factor f ← 1
while f > f_min:
    propose move, accept with prob min(1, g(E)/g(E'))
    update g(E_new) ← g(E_new) · exp(f),  H(E_new) ← H(E_new) + 1
    if H(E) is "flat" (e.g. 80% of mean):
        f ← f/2,  reset H(E) ← 0
return g(E)
```

Given `g(E)`, compute partition functions and observables at *any*
temperature in post-processing.  Particularly useful for first-order
transitions where canonical samplers tunnel poorly.

## 8. Validation strategies

Always test new code against:

* **Analytical limits.**  `T\to 0` ground state, `T\to\infty` mean energy,
  exact solutions for 1-D Ising or 2-D Ising on small lattices, Onsager
  `T_c = 2J/\ln(1+\sqrt 2)`.
* **Detailed balance.**  Run on a small system, histogram visited states,
  compare with `e^{-\beta E}/Z`.
* **Symmetry breaking.**  After symmetrising estimators (e.g.
  `\langle |M|\rangle` rather than `\langle M\rangle`), confirm the symmetric value at
  high `T` is zero within errors.
* **Reproducibility.**  Same seed ⇒ same output, bit-for-bit.

## 9. Implementation patterns

* **Local energy updates.**  Compute `\Delta E` from neighbours of the
  flipped site; never recompute total energy.  Maintain a running total
  for observables.
* **Lookup tables.**  For Ising, `e^{-\beta\Delta E}` has only a handful of
  possible values — pre-tabulate.
* **Random site selection vs sequential sweeps.**  Both satisfy detailed
  balance.  Sequential is slightly faster (better cache, no RNG call per
  site choice).
* **Cluster algorithms.**  Use a union-find data structure or a stack-based
  flood fill; both run in amortised `O(N)` per sweep.
* **Vectorisation.**  Checkerboard updates: in models with bipartite
  couplings, all "black" sites can be updated in parallel.  Friendly to
  SIMD and GPUs.
* **Measurement frequency.**  Measure every `O(\tau_{\rm int})` sweeps,
  *not* every sweep — measurement is typically the most expensive step.
* **Memory layout.**  Store spins as `Int8` (or even bit-packed) to fit
  large lattices in cache.

## 10. Profiling and scaling tests

Before optimising, measure where time is spent.  Three traps:

1. **RNG dominance** at `N\!\sim\!10^4` — use a counter-based RNG and a
   tabulated `e^{-\beta\Delta E}`.
2. **Cache misses** at `N\!\gtrsim\!10^6` — switch to bit-packed spins or
   row-major iteration.
3. **Reallocation in the inner loop** — preallocate temporaries; in
   Julia, avoid `global` and ensure type stability.

A good scaling test: plot wall-clock per sweep against `N`.  For Metropolis
Ising this should be linear; deviations point to memory traffic issues.

## 11. Reproducibility checklist

Before claiming a result:

- [ ] seed logged
- [ ] equilibration discarded (with justification)
- [ ] `\tau_{\rm int}` reported for each observable
- [ ] error bars computed via binning or jackknife, not `\sigma/\sqrt{N}`
- [ ] system-size scan (at least three `L`) for any thermodynamic claim
- [ ] code passes a small-system comparison against exact diagonalisation
- [ ] one run reproduced bit-for-bit from the same seed
