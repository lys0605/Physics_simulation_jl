# Famous Monte Carlo Models Across Fields

Monte Carlo methods spread quickly out of physics.  The same Markov-chain
machinery describes opinion formation, biological evolution, market
fluctuations, traffic flow, and machine learning.  This note collects the
canonical "named" stochastic models in each field, with the central
equation, and points to a companion Julia script (if one exists in
`scripts/monte_carlo/`).

The list is curated, not exhaustive.  Each row is a model worth knowing
because it appears in multiple textbooks, has been validated against data,
or is a recurring teaching example.

---

## Social / behavioural sciences

| Model | Year | What it shows | Update rule (sketch) |
|---|---|---|---|
| **Schelling segregation** | 1971 | Sharp segregation emerges from mild local preference for like-coloured neighbours | An unhappy agent (fewer than `τ N_neighbours` similar neighbours) moves to a random empty cell |
| **Voter model** | 1973 | Consensus or coexistence depending on dimension and graph | Site `i` adopts neighbour `j`'s opinion at rate 1; `\sum_j (s_j - s_i) / \text{deg}(i)` |
| **Sznajd model** | 2000 | Two-agent persuasion → opinion clustering | If `s_i = s_{i+1}`, neighbours of the pair adopt that opinion |
| **Axelrod culture model** | 1997 | Multidimensional traits, "homophily" → cultural domains | Site copies a feature of a similar neighbour with prob ∝ overlap |
| **Hegselmann-Krause** | 2002 | Bounded-confidence opinion dynamics | `o_i ← \text{mean}\{o_j : |o_j - o_i| < \epsilon\}` |

Companion script: [`schelling.jl`](../../scripts/monte_carlo/schelling.jl).

These are the bread and butter of *computational sociology* and
*sociophysics*.  Once translated to a Metropolis-style update, the same
analysis machinery (autocorrelation, finite-size scaling, phase
transitions) applies.

---

## Biology and evolution

| Model | Year | Domain | Central quantity |
|---|---|---|---|
| **Wright-Fisher** | 1930s | Population genetics | Allele frequency `p_{t+1} \sim \text{Binom}(N, p_t)/N` |
| **Moran process** | 1958 | Population genetics, continuous-time | Birth-death updates, fixation probability `(1-r^{-i})/(1-r^{-N})` for selection coefficient `r` |
| **Kimura's neutral theory** | 1968 | Molecular evolution | Drift dominates over selection for `\|s\| \ll 1/N` |
| **Galton-Watson branching** | 1874 | Family extinction, virus spread, polymerase chain reaction | Offspring `Z_{t+1} = \sum_i X_i^{(t)}` |
| **SIR / SIS epidemics** | 1927 | Disease spread on graphs | Stochastic transitions S→I (rate `β`), I→R (rate `γ`) |
| **Lotka-Volterra (stochastic)** | 1925 | Predator-prey | Gillespie SSA on rabbit-fox reactions |
| **Eigen's quasispecies** | 1971 | Molecular evolution under mutation | Eigenproblem for replication+mutation matrix |
| **Neutral evolution / Markov mutation** | various | Sequence evolution | Substitution rate matrices (Jukes-Cantor, HKY, GTR) |

Companion scripts: [`wright_fisher.jl`](../../scripts/monte_carlo/wright_fisher.jl),
[`gillespie.jl`](../../scripts/monte_carlo/gillespie.jl).

The link to MC: every one of these can be cast as a Markov chain or
master equation, and stationary properties become MC estimators.

---

## Economics and finance ("econophysics")

| Model | Year | What it shows | Mechanism |
|---|---|---|---|
| **Kinetic exchange (Chakraborti-Chakrabarti)** | 2000 | Wealth distribution → Boltzmann-like tail | Pairwise wealth-conservation exchange with random split `\epsilon` |
| **Yard-sale / Bouchaud-Mézard** | 2000 | Wealth condensation | Multiplicative noise + redistribution |
| **Minority game** | 1997 | Self-organising market participants | Agents pick the side with fewer recent winners |
| **Cont-Bouchaud herd model** | 2000 | Fat-tailed asset returns | Random clustering of traders → percolation-like volatility |
| **Geometric Brownian motion** | 1900 (Bachelier) | Option pricing | `dS/S = \mu\,dt + \sigma\,dW` |
| **Heston stochastic volatility** | 1993 | Smile in implied volatility | Mean-reverting variance + price coupled |
| **Agent-based macro (Kirman ants)** | 1993 | Herding-induced bimodality | Two-state agents flip per neighbour count |

The connection to MC: option pricing via random sampling of terminal
payoffs (Boyle, 1977) is **the** original use of MC in finance.  Wealth
distribution toy models are direct analogues of the Maxwell-Boltzmann
distribution.

---

## Chemistry and reaction kinetics

| Model | Year | Domain | Algorithm |
|---|---|---|---|
| **Gillespie SSA (exact)** | 1976 | Well-stirred chemical kinetics | Sample (time-to-next-event, which-event) per master equation |
| **τ-leaping** | 2001 | Approximate, faster SSA | Advance by `\tau`, draw Poisson counts of reactions |
| **Kinetic MC for surface growth** | 1980s | Crystal / film growth | Lattice cells with deposition / hopping / desorption rates |
| **Reactive MC / FFS (forward flux sampling)** | 2000s | Rare events in reaction networks | Sequential sampling of trajectories crossing interfaces |

Companion script: [`gillespie.jl`](../../scripts/monte_carlo/gillespie.jl).

---

## Computer science and machine learning

| Model | Year | What for | Key idea |
|---|---|---|---|
| **Simulated annealing** | 1983 | Combinatorial optimisation | Metropolis on `\exp(-E/T)` with `T` cooled slowly |
| **Boltzmann machine / RBM** | 1985 | Generative model | Gibbs sampling on hidden + visible binary units |
| **Markov-chain Monte Carlo for Bayesian inference** | 1990s | Posterior sampling | Same Metropolis-Hastings, different `\pi` |
| **HMC / NUTS** | 1987 | Continuous parameter sampling | Hamiltonian dynamics + accept/reject |
| **TD-Gammon, MCTS** | 1990s, 2006 | Reinforcement learning, game AI | Tree search with random rollouts |
| **Stochastic gradient Langevin** | 2011 | Bayesian deep learning | `\theta ← \theta - \eta\nabla L + \sqrt{2\eta T}\xi` |

The relationship to physics MC is **literal** — most ML samplers were
invented by physicists, and concepts like detailed balance, mixing time,
and free energy carry across.

---

## Statistical geometry and dynamical systems

| Model | Year | What it shows | Companion script |
|---|---|---|---|
| **Site / bond percolation** | 1957 | Geometric phase transition; cluster scaling | [`percolation.jl`](../../scripts/monte_carlo/percolation.jl) |
| **Forest fire (Drossel-Schwabl)** | 1992 | Self-organised criticality, power-law clusters | [`forest_fire.jl`](../../scripts/monte_carlo/forest_fire.jl) |
| **Bak-Sneppen evolution** | 1993 | SOC, punctuated equilibria | minimum fitness replaced + neighbours |
| **Eden growth / DLA** | 1961, 1981 | Fractal aggregation | random walkers stick on contact |
| **Random walk on graphs (PageRank)** | 1998 | Network ranking | Stationary distribution of teleporting walk |

Companion script: [`random_walk.jl`](../../scripts/monte_carlo/random_walk.jl).

---

## Physics (for completeness)

These are covered in the other notes; we list them so readers can map MC
back from cross-field analogues:

| Topic | Method | Companion script |
|---|---|---|
| Monte-Carlo integration | importance sampling | [`pi_estimation.jl`](../../scripts/monte_carlo/pi_estimation.jl) |
| Ising critical point | Metropolis + Wolff | [`ising_metropolis.jl`](../../scripts/monte_carlo/ising_metropolis.jl), [`ising_wolff.jl`](../../scripts/monte_carlo/ising_wolff.jl) |
| Autocorrelation / error analysis | binning, jackknife | [`autocorrelation.jl`](../../scripts/monte_carlo/autocorrelation.jl) |
| Quantum harmonic oscillator | VMC | [`vmc_harmonic.jl`](../../scripts/monte_carlo/vmc_harmonic.jl) |

---

## A unified view

The common pattern across **every** entry in this note:

1. Define a *state space* of configurations.
2. Define a *stochastic update rule* (Markov chain).
3. Identify the *stationary distribution* or *transient observables* of
   interest.
4. Estimate them by averaging over many sampled trajectories.

When the update rule satisfies detailed balance with respect to a target
distribution `π`, you have classical MC.  When it does not — when the
dynamics itself is the object of study (segregation, epidemics, traffic) —
the same tooling estimates *non-equilibrium* averages, return times,
extinction probabilities, and so on.

The Julia scripts in `scripts/monte_carlo/` are intentionally short
(≤ 100 lines each).  Pick one whose physics you do not yet know — try
to predict the outcome, then run it.
