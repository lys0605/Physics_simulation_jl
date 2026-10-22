# Green's Functions for Tight-Binding Models
## Theory, Computation, and Topological Systems

---

## Table of Contents

1. [Why Green's Functions?](#1-why-greens-functions)
2. [Definitions and Representations](#2-definitions-and-representations)
3. [The Spectral Function](#3-the-spectral-function)
4. [Tight-Binding Models on 2D Lattices](#4-tight-binding-models-on-2d-lattices)
5. [The QWZ Chern Insulator](#5-the-qwz-chern-insulator)
6. [Topological Invariants: Chern Number](#6-topological-invariants-chern-number)
7. [Bulk Green's Function (Periodic BC)](#7-bulks-greens-function-periodic-bc)
8. [Ribbon Green's Function (Open BC in x)](#8-ribbon-greens-function-open-bc-in-x)
9. [Surface Green's Function (Semi-infinite BC)](#9-surface-greens-function-semi-infinite-bc)
10. [Hartree-Fock Self-Energy and Interaction-Driven Phase Transitions](#10-hartree-fock-self-energy-and-interaction-driven-phase-transitions)
11. [Simulation Details and Convergence](#11-simulation-details-and-convergence)
12. [Physical Quantities Summary](#12-physical-quantities-summary)

---

## 1. Why Green's Functions?

The single-particle retarded Green's function

$$G^R(\mathbf{k}, \omega) = \langle\langle c_{\mathbf{k}} ; c^\dagger_{\mathbf{k}} \rangle\rangle^R_\omega$$

is arguably the most information-dense object in condensed-matter physics. From it one can extract:

| Quantity | Formula |
|---|---|
| Spectral function | $A(\mathbf{k}, \omega) = -\frac{1}{\pi}\,\mathrm{Im\,Tr}\,G^R(\mathbf{k}, \omega)$ |
| Density of states | $\rho(\omega) = \frac{1}{N_k}\sum_\mathbf{k} A(\mathbf{k}, \omega)$ |
| Quasiparticle lifetime | $\tau^{-1} = -2\,\mathrm{Im}\,\Sigma^R(\mathbf{k}, \omega)$ |
| Linear-response functions | $\chi(\omega) \propto \int d\omega'\, G^R(\omega+\omega')\,G^A(\omega')$ |
| Self-energy / interactions | Dyson: $G = G_0 + G_0 \Sigma G$ |

The Green's function encodes the full single-particle excitation spectrum (via the spectral function) and, through diagrammatic perturbation theory, can be systematically dressed with interactions. For non-interacting systems it reduces to the resolvent of the Hamiltonian — exactly what we compute here.

---

## 2. Definitions and Representations

### 2.1 Retarded and Advanced Green's Functions

In the energy domain the retarded (+) and advanced (–) Green's functions for a Hamiltonian $H$ are

$$G^\pm(\omega) = \lim_{\eta\to 0^+} \bigl[(\omega \pm i\eta)\,I - H\bigr]^{-1}$$

where $\eta > 0$ is an infinitesimal positive broadening. Numerically we keep $\eta$ finite:

$$G(\omega) \equiv G^R(\omega;\eta) = \bigl[(\omega + i\eta)\,I - H\bigr]^{-1}$$

The two are related by $G^A = (G^R)^\dagger$ (Hermitian conjugate).

### 2.2 Lehmann Representation

For a non-interacting system with eigenstates $|n\rangle$, eigenenergies $E_n$:

$$G_{ab}(\omega) = \sum_n \frac{\langle a | n\rangle\langle n | b\rangle}{\omega + i\eta - E_n}$$

This shows that $G$ has poles at the exact eigenvalues. The spectral weight carried by state $|n\rangle$ in the $\{a, b\}$ channel is $\langle a|n\rangle\langle n|b\rangle$.

### 2.3 Dyson Equation

When a perturbation $\Sigma$ (the self-energy) dresses the bare propagator $G_0$:

$$G = G_0 + G_0\,\Sigma\,G \quad\Longleftrightarrow\quad G^{-1} = G_0^{-1} - \Sigma$$

For a semi-infinite geometry, $\Sigma = T\,G_s\,T^\dagger$ is the self-energy due to the semi-infinite bulk, leading to the self-consistent equation solved by the López Sancho algorithm (Section 9).

---

## 3. The Spectral Function

### 3.1 Definition

$$A(\mathbf{k}, \omega) = -\frac{1}{\pi}\,\mathrm{Im\,Tr}\,G^R(\mathbf{k}, \omega;\eta)$$

For a non-interacting system this evaluates to a sum of Lorentzians:

$$A(\mathbf{k}, \omega) = \frac{1}{\pi}\sum_n \frac{\eta}{(\omega - E_n(\mathbf{k}))^2 + \eta^2}$$

In the limit $\eta \to 0^+$ each Lorentzian becomes a delta function: $A(\mathbf{k},\omega) \to \sum_n \delta(\omega - E_n(\mathbf{k}))$.

### 3.2 Sum Rule

The spectral function obeys the exact sum rule

$$\int_{-\infty}^{\infty} A(\mathbf{k}, \omega)\,d\omega = 1 \qquad \text{(per site, per spin)}$$

which is a useful numerical check: for finite $\eta$ and a sufficiently wide $\omega$ window, the discrete sum over the $\omega$ grid should approach 1 after rescaling by $d\omega$.

### 3.3 Density of States

$$\rho(\omega) = \frac{1}{(2\pi)^2}\int_\mathrm{BZ} A(\mathbf{k},\omega)\,d^2k \approx \frac{1}{N_k^2}\sum_{\mathbf{k}} A(\mathbf{k},\omega)$$

The factor $(2\pi)^{-2}$ is absorbed into the uniform mesh weight $N_k^{-2}$.

---

## 4. Tight-Binding Models on 2D Lattices

### 4.1 Real-Space Hamiltonian

A general nearest-neighbour tight-binding model on the square lattice with $m$ orbitals per site has the Hamiltonian

$$H = \sum_{\mathbf{R}} \mathbf{c}^\dagger_{\mathbf{R}} H_0 \mathbf{c}_{\mathbf{R}} + \sum_{\mathbf{R}}\left(\mathbf{c}^\dagger_{\mathbf{R}+\hat{x}} T_x \mathbf{c}_{\mathbf{R}} + \mathbf{c}^\dagger_{\mathbf{R}+\hat{y}} T_y \mathbf{c}_{\mathbf{R}} + \mathrm{h.c.}\right)$$

where $H_0$ is the on-site (intra-cell) Hamiltonian and $T_{x,y}$ are the inter-cell hopping matrices.

### 4.2 Bloch Hamiltonian (Periodic BC)

Under the Fourier transform $\mathbf{c}_{\mathbf{R}} = N^{-1/2}\sum_\mathbf{k} e^{i\mathbf{k}\cdot\mathbf{R}}\mathbf{c}_\mathbf{k}$, the Hamiltonian decouples:

$$H = \sum_\mathbf{k} \mathbf{c}^\dagger_\mathbf{k}\,\mathcal{H}(\mathbf{k})\,\mathbf{c}_\mathbf{k}$$

$$\mathcal{H}(\mathbf{k}) = H_0 + T_x\,e^{-ik_x} + T_x^\dagger\,e^{ik_x} + T_y\,e^{-ik_y} + T_y^\dagger\,e^{ik_y}$$

For the square lattice, $k_x, k_y \in (-\pi, \pi]$ (first Brillouin zone, lattice spacing $a = 1$).

---

## 5. The QWZ Chern Insulator

### 5.1 Model Definition

The Qi-Wu-Zhang (QWZ) model (Qi, Wu, Zhang 2006) is a minimal 2-band Chern insulator on the square lattice. In the Bloch basis $|A\rangle, |B\rangle$ (two sublattices / orbitals):

$$\mathcal{H}(\mathbf{k}) = \sin k_x\,\sigma_x + \sin k_y\,\sigma_y + (m + \cos k_x + \cos k_y)\,\sigma_z$$

where $\sigma_{x,y,z}$ are Pauli matrices and $m$ is a tunable mass parameter.

The band energies are $E_\pm(\mathbf{k}) = \pm |\mathbf{d}(\mathbf{k})|$ with

$$\mathbf{d} = (\sin k_x,\; \sin k_y,\; m + \cos k_x + \cos k_y)$$

$$|\mathbf{d}|^2 = \sin^2 k_x + \sin^2 k_y + (m + \cos k_x + \cos k_y)^2$$

The spectrum is gapped whenever $\mathbf{d} \neq 0$ for all $\mathbf{k}$.  The gap closes at:

| $m$ | Gapless point |
|---|---|
| $m = 0$ | $\mathbf{k} = (0, 0)$ |
| $m = -2$ | $\mathbf{k} = (\pi, 0), (0,\pi)$ |
| $m = 2$ | $\mathbf{k} = (\pi, \pi)$ |

### 5.2 Phase Diagram

$$C = \begin{cases} +1 & -2 < m < 0 \\ -1 & 0 < m < 2 \\ 0 & |m| > 2 \end{cases}$$

For the simulations we use $m = -1$ (topological, $C = +1$) and $m = 3$ (trivial, $C = 0$).

### 5.3 Real-Space Hoppings

Expanding $\mathcal{H}(\mathbf{k})$ in real space, the hopping matrices are:

$$T_x = \frac{\sigma_z - i\sigma_x}{2} = \frac{1}{2}\begin{pmatrix}1 & -i \\ -i & -1\end{pmatrix}$$

$$T_y = \frac{\sigma_z - i\sigma_y}{2} = \frac{1}{2}\begin{pmatrix}1 & -1 \\ 1 & -1\end{pmatrix}$$

The on-site term is $H_0(k_y) = m\,\sigma_z + 2T_y\cos k_y + \ldots$ (layer-dependent when $k_y$ is still a good quantum number in the ribbon geometry, see Section 8).

---

## 6. Topological Invariants: Chern Number

### 6.1 Berry Curvature

The Berry curvature of the occupied (lower) band $|u_-(\mathbf{k})\rangle$ is

$$\Omega(\mathbf{k}) = -2\,\mathrm{Im}\sum_{n \neq -} \frac{\langle u_-|\partial_{k_x}\mathcal{H}|u_n\rangle\langle u_n|\partial_{k_y}\mathcal{H}|u_-\rangle}{(E_n - E_-)^2}$$

For the 2-band QWZ model this simplifies to

$$\Omega(\mathbf{k}) = \frac{1}{2}\hat{\mathbf{d}}\cdot\left(\frac{\partial\hat{\mathbf{d}}}{\partial k_x} \times \frac{\partial\hat{\mathbf{d}}}{\partial k_y}\right), \qquad \hat{\mathbf{d}} = \mathbf{d}/|\mathbf{d}|$$

### 6.2 Chern Number

$$C = \frac{1}{2\pi}\int_\mathrm{BZ} \Omega(\mathbf{k})\,d^2k \in \mathbb{Z}$$

$C$ counts how many times $\hat{\mathbf{d}}$ wraps around the unit sphere as $\mathbf{k}$ covers the BZ.

### 6.3 Fukui–Hatsugai–Suzuki (FHS) Lattice Formula

For a discretized BZ with $N_k \times N_k$ plaquettes, the FHS formula (2005) gives an exact integer:

$$C = \frac{1}{2\pi}\sum_{\mathbf{k}} \mathrm{Im}\ln\bigl[U_x(\mathbf{k})\,U_y(\mathbf{k}+\hat{x})\,U_x^*(\mathbf{k}+\hat{y})\,U_y^*(\mathbf{k})\bigr]$$

where the link variable (Wilson line) is

$$U_\mu(\mathbf{k}) = \frac{\langle u(\mathbf{k})|u(\mathbf{k}+\delta\mu)\rangle}{|\langle u(\mathbf{k})|u(\mathbf{k}+\delta\mu)\rangle|}$$

**Why use FHS?** The continuum Berry-phase integral requires a smooth gauge (fixing the phase of $|u(\mathbf{k})\rangle$ globally), which is impossible in a non-trivial topology. The FHS construction is gauge-invariant by design and always returns an exact integer, even for coarse meshes.

**Convergence:** For the QWZ model, $N_k = 30$ gives $C = \pm 1$ or $0$ exactly.  The eigenvectors $|u(\mathbf{k})\rangle$ can be chosen as the output of `eigvecs(H(k))` with no phase fixing required.

---

## 7. Bulk Green's Function (Periodic BC)

### 7.1 Setup

With periodic BC in both $x$ and $y$, $\mathbf{k}$ is a conserved quantum number and

$$G(\mathbf{k}, \omega) = \bigl[(\omega + i\eta)\,I - \mathcal{H}(\mathbf{k})\bigr]^{-1}$$

This is a $2\times 2$ matrix inversion (2 orbitals). The spectral function plotted along a high-symmetry path $\Gamma \to X \to M \to \Gamma$ visualises both bulk bands simultaneously and shows the gap structure.

### 7.2 High-Symmetry Points

On the square BZ (lattice spacing $a = 1$):

$$\Gamma = (0,0), \quad X = (\pi, 0), \quad M = (\pi, \pi)$$

The path $\Gamma \to X \to M \to \Gamma$ traverses all distinct symmetry points and passes near the gap-closing momenta at $m = 0, \pm 2$.

### 7.3 Spectral Function Along the Path

For the $n$-th segment of the path with $N_\mathrm{seg}$ points, one evaluates $A(k_i, \omega)$ for all $(\omega_j, k_i)$ pairs on a grid — $O(N_k N_\omega)$ matrix inversions, each $O(d^3)$ for an $d\times d$ Hamiltonian. Here $d = 2$, so cost is negligible.

The closed-form version for 2 bands is especially fast:

$$A(\mathbf{k},\omega) = \frac{\eta}{\pi}\left[\frac{1}{(\omega-E_+)^2+\eta^2} + \frac{1}{(\omega-E_-)^2+\eta^2}\right]$$

---

## 8. Ribbon Green's Function (Open BC in x)

### 8.1 Ribbon Hamiltonian

Periodic BC in $y$ ($k_y$ conserved); open (hard-wall) BC at $x = 1$ and $x = N_x$. With $N_x$ unit cells in $x$, each with 2 orbitals, the Hamiltonian at each $k_y$ is a $2N_x \times 2N_x$ block-tridiagonal matrix:

$$H^\mathrm{ribbon}(k_y) = \begin{pmatrix} H_0(k_y) & T^\dagger & & \\ T & H_0(k_y) & T^\dagger & \\ & \ddots & \ddots & \ddots \\ & & T & H_0(k_y) \end{pmatrix}$$

where $H_0(k_y) = (m + \cos k_y)\sigma_z - \sin k_y\,\sigma_y$ is the on-site term in the $k_y$-Fourier basis and $T = T_x$ is the hopping in $x$.

**Boundary conditions:** The hard-wall BC means no hopping matrix outside the first/last block — the $(1, N_x)$ and $(N_x, 1)$ blocks are absent.

### 8.2 Band Structure

Full diagonalisation of $H^\mathrm{ribbon}(k_y)$ gives $2N_x$ bands $\{E_n(k_y)\}$. In the topological phase ($C = \pm 1$), two edge bands traverse the bulk gap and connect the valence and conduction bands — these are the **topological edge states** guaranteed by the bulk-boundary correspondence:

$$N_\mathrm{edge} = |C|$$

edge modes crossing the gap (per edge, per spin).

### 8.3 Layer-Resolved DOS

The layer-resolved density of states (LDOS) at layer $j$ and energy $\omega$ is

$$\mathrm{LDOS}(j, \omega; k_y) = \frac{1}{\pi}\sum_n \bigl(|V_{2j-1,n}|^2 + |V_{2j,n}|^2\bigr)\frac{\eta}{(\omega - E_n)^2 + \eta^2}$$

where $V_{mn}$ is the $m$-th component of the $n$-th eigenvector.  Bright weight at $j = 1$ or $j = N_x$ confirms edge localisation of in-gap states.

**Why LDOS?** A global DOS cannot distinguish bulk from edge.  The LDOS reveals which layers carry the in-gap spectral weight and confirms the spatial profile of the edge mode.

### 8.4 Finite-Size Effects

The ribbon has width $N_x$. Edge states on opposite edges (at $j=1$ and $j=N_x$) hybridise exponentially:

$$\Delta E \sim e^{-N_x/\xi}$$

where $\xi \sim v_F/\Delta$ is the edge-state localisation length ($\Delta$ = bulk gap, $v_F$ = edge velocity). For $N_x = 40$ and the QWZ gap $\Delta \approx 2$ this splitting is negligible.

---

## 9. Surface Green's Function (Semi-infinite BC)

### 9.1 Geometry and Dyson Equation

Periodic BC in $y$; semi-infinite in $x$ with the surface at layer $j = 1$ extending to $j \to +\infty$. The surface retarded GF $G_s(k_y, \omega)$ satisfies the self-consistent Dyson equation

$$G_s = \bigl[(\omega+i\eta)\,I - H_0(k_y) - \Sigma_s\bigr]^{-1}$$

where the bulk self-energy is

$$\Sigma_s = T\,G_s\,T^\dagger$$

This is a fixed-point (self-consistent) equation for the $2 \times 2$ matrix $G_s$.

### 9.2 Naive Fixed-Point Iteration

Start with $G_s^{(0)} = [zI - H_0]^{-1}$ ($z = \omega + i\eta$) and iterate:

$$G_s^{(n+1)} = \bigl[zI - H_0 - T\,G_s^{(n)}\,T^\dagger\bigr]^{-1}$$

**Convergence rate:** Each iteration effectively adds one more layer to the semi-infinite stack. The error decreases as

$$\|G_s^{(n)} - G_s\| \sim \left|\frac{t}{\eta}\right|^{2n}$$

where $t$ is the hopping amplitude.  For small $\eta \ll t$, this requires $O(t^2/\eta^2)$ iterations — catastrophically slow.

**Practical use:** The naive method works for $\eta \gtrsim 0.1$ but becomes intractable for $\eta \lesssim 0.01$.

### 9.3 López Sancho Decimation Algorithm

The López Sancho (1984) renormalisation-group algorithm doubles the effective slab depth at each step, achieving $O(\log_2(1/\eta))$ convergence.

**Initialise:**

$$\varepsilon_s = \varepsilon = H_0(k_y), \quad \alpha = T, \quad \beta = T^\dagger$$

**Iterate until $\|\alpha\| < \varepsilon_\mathrm{tol}$:**

$$g_n = \bigl[zI - \varepsilon_n\bigr]^{-1}$$

$$\varepsilon_s^{(n+1)} = \varepsilon_s^{(n)} + \alpha_n\,g_n\,\beta_n$$

$$\varepsilon^{(n+1)} = \varepsilon^{(n)} + \alpha_n\,g_n\,\beta_n + \beta_n\,g_n\,\alpha_n$$

$$\alpha_{n+1} = \alpha_n\,g_n\,\alpha_n, \quad \beta_{n+1} = \beta_n\,g_n\,\beta_n$$

**Extract surface GF:**

$$G_s = \bigl[zI - \varepsilon_s\bigr]^{-1}$$

**Why does this work?**  $\alpha_n$ represents the effective hopping across $2^n$ unit cells; once this vanishes, the surface layer is fully decoupled from the bulk. The sequence $\|\alpha_0\| \to \|\alpha_1\| \to \cdots$ decreases geometrically:

$$\|\alpha_n\| \sim \left(\frac{t^2}{\eta^2}\right)^{2^n - 1} \cdot \|\alpha_0\|$$

so the number of iterations required scales as $\log_2\log(t/\eta) \approx \log_2(1/\eta)$ — far fewer than the naive $O(\eta^{-2})$.

**Typical iteration count:**

| $\eta$ | Naive iters | López Sancho iters |
|---|---|---|
| $0.2$ | ~25 | ~4 |
| $0.05$ | ~400 | ~7 |
| $0.01$ | ~10 000 | ~10 |
| $0.001$ | ~$10^6$ | ~13 |

### 9.4 Surface Spectral Function

$$A_s(k_y, \omega) = -\frac{1}{\pi}\,\mathrm{Im\,Tr}\,G_s(k_y, \omega)$$

This is the spectral function projected onto the surface layer. Unlike the ribbon, there is **no finite-size hybridisation** between the two edges — the system is truly semi-infinite. Topological edge modes appear as sharp, perfectly in-gap features even at small $\eta$.

---

## 10. Hartree-Fock Self-Energy and Interaction-Driven Phase Transitions

### 10.1 Two-Body Interaction

We add an inter-orbital density-density (Hubbard) interaction to the QWZ model:

$$H_\text{int} = U \sum_{\mathbf{R}} n_{\mathbf{R},1}\, n_{\mathbf{R},2}$$

where $n_{\mathbf{R},a} = c^\dagger_{\mathbf{R},a} c_{\mathbf{R},a}$ is the number operator for orbital $a \in \{1,2\}$ on site $\mathbf{R}$.  In momentum space, a contact interaction $U$ (constant in $\mathbf{q}$) is

$$H_\text{int} = \frac{U}{N}\sum_{\mathbf{k},\mathbf{k}',\mathbf{q}} c^\dagger_{\mathbf{k}+\mathbf{q},1}\,c^\dagger_{\mathbf{k}'-\mathbf{q},2}\,c_{\mathbf{k}',2}\,c_{\mathbf{k},1}$$

### 10.2 Hartree-Fock Decoupling

Mean-field (Hartree-Fock) decoupling of the quartic interaction retains two classes of contractions:

**Hartree (direct):** $c^\dagger_{\mathbf{k},a}\,c_{\mathbf{k},a} \to \langle n_a \rangle = \rho_{aa}$, giving a diagonal self-energy

$$\Sigma^H_{11} = U\rho_{22}, \qquad \Sigma^H_{22} = U\rho_{11}$$

**Fock (exchange):** $c^\dagger_{\mathbf{k},1}\,c_{\mathbf{k},2} \to \rho_{12}$, giving off-diagonal corrections

$$\Sigma^F_{12} = -U\rho_{12}, \qquad \Sigma^F_{21} = -U\rho_{21}$$

The full $2\times 2$ HF self-energy for a contact interaction is **k-independent**:

$$\Sigma^{\text{HF}} = \begin{pmatrix} U\rho_{22} & -U\rho_{12} \\ -U\rho_{21} & U\rho_{11} \end{pmatrix}$$

### 10.3 Self-Consistent Field (SCF) Equations

The density matrix is determined self-consistently from the lower band of the HF Hamiltonian $\mathcal{H}^{\text{HF}}(\mathbf{k}) = \mathcal{H}_0(\mathbf{k}) + \Sigma^{\text{HF}}$:

$$\rho_{ab} = \frac{1}{N_k^2}\sum_\mathbf{k} \langle a | u_-(\mathbf{k})\rangle\langle u_-(\mathbf{k}) | b\rangle$$

where $|u_-(\mathbf{k})\rangle$ is the lower eigenstate of $\mathcal{H}^{\text{HF}}(\mathbf{k})$.  The SCF loop reads:

1. Start with $\Sigma^{(0)} = 0$.
2. Compute $\rho$ from the lower band of $\mathcal{H}_0 + \Sigma^{(n)}$.
3. Compute $\Sigma^{(n+1)}_\text{new}$ from the HF formula.
4. Mix: $\Sigma^{(n+1)} = \alpha\,\Sigma^{(n+1)}_\text{new} + (1-\alpha)\,\Sigma^{(n)}$ (linear mixing, $\alpha \approx 0.4$).
5. Repeat until $\|\Sigma^{(n+1)} - \Sigma^{(n)}\| < \varepsilon_\text{tol}$.

Linear mixing is necessary for convergence near the phase transition where the fixed-point map has a near-unit eigenvalue.

### 10.4 Vanishing of the Fock Off-Diagonal (Inversion Symmetry)

The QWZ model is invariant under independent sign reversal of each momentum component: $k_x \to -k_x$ maps $d_x \to -d_x$, $d_y \to d_y$, $d_z \to d_z$. Since $d_x$ is odd under this reflection, the lower-band eigenvector contribution to $\rho_{12}$ integrates to zero:

$$\rho_{12} = \frac{1}{N_k^2}\sum_\mathbf{k} \langle 1 | u_-(\mathbf{k})\rangle\langle u_-(\mathbf{k}) | 2\rangle = -\frac{1}{2}\left\langle \frac{d_x + id_y}{|\mathbf{d}|}\right\rangle_\text{BZ} = 0$$

because $d_x/|\mathbf{d}|$ is odd in $k_x$ and integrates to zero. Consequently **the Fock term vanishes exactly** and $\Sigma^{\text{HF}}$ is purely diagonal:

$$\Sigma^{\text{HF}} = \operatorname{diag}(U\rho_{22},\; U\rho_{11})$$

### 10.5 Orbital Polarization and Effective Mass Renormalization

The diagonal self-energy acts as a momentum-independent staggered potential, equivalent to renormalizing the topological mass parameter:

$$m_\text{eff}(U) = m + \frac{U(\rho_{22} - \rho_{11})}{2}$$

(the $U\cdot I$ part shifts the chemical potential only).  For $m = -1$, the bare occupancies at $U=0$ satisfy $\rho_{11} > \rho_{22}$ (orbital 1 is preferentially occupied throughout most of the BZ where $d_z < 0$), so $\rho_{22}-\rho_{11} < 0$ and **$m_\text{eff}$ decreases with $U$**.

This creates a positive feedback: higher $U$ increases the Hartree potential on orbital 2 ($\Sigma_{22} = U\rho_{11}$), which pushes electrons out of orbital 2 and into orbital 1, further increasing $\rho_{11}$ and decreasing $m_\text{eff}$ — an **orbital polarization instability** analogous to the Stoner instability in ferromagnetism.

### 10.6 Interaction-Driven Topological Phase Transition

The topological phase boundary of the QWZ model is at $m_\text{eff} = -2$ (where the bulk gap closes at $\mathbf{k} = (\pi, 0)$ and $(0, \pi)$).  Setting $m_\text{eff} = -2$ gives the critical condition

$$U_c\,[\rho_{11}(U_c) - \rho_{22}(U_c)] = 2$$

which must be solved self-consistently since $\rho_{aa}$ depends on $U$.  Numerically for $m=-1$: $U_c \approx 2.46$.

For $U < U_c$: $m_\text{eff} > -2$, topological phase $C = -1$, in-gap edge modes present.  
For $U > U_c$: $m_\text{eff} < -2$, trivial phase $C = 0$, gap re-opens without edge modes.

### 10.7 HF Green's Function and Spectral Function

The interacting retarded Green's function in HF is

$$G^{\text{HF}}(\mathbf{k}, \omega) = \bigl[(\omega+i\eta)\,I - \mathcal{H}_0(\mathbf{k}) - \Sigma^{\text{HF}}\bigr]^{-1}$$

Since $\Sigma^{\text{HF}}$ is real and static, the Dyson equation gives a renormalized band structure with no imaginary self-energy. The spectral function

$$A^{\text{HF}}(\mathbf{k},\omega) = -\frac{1}{\pi}\,\mathrm{Im\,Tr}\,G^{\text{HF}}(\mathbf{k},\omega)$$

remains a sum of Lorentzians (width $2\eta$) centred at the renormalized eigenvalues.  The quasiparticle weight $Z = 1$ exactly within HF — spectral weight is not transferred to incoherent backgrounds. To obtain $Z < 1$ (finite quasiparticle lifetime), one must go beyond HF to the **second-order Born approximation**:

$$\Sigma^{(2)}(\mathbf{k},\omega) = \frac{U^2}{N^2}\sum_{\mathbf{k}',\mathbf{q}} \frac{f(E_-(\mathbf{k}'+\mathbf{q}))(1-f(E_-(\mathbf{k}')))}{\omega - E_-(\mathbf{k}') - E_-(\mathbf{k}+\mathbf{k}'+\mathbf{q}) + i\eta} + \ldots$$

which has a non-zero imaginary part and gives $Z = [1 - \partial\,\mathrm{Re}\,\Sigma^{(2)}/\partial\omega|_{\omega=0}]^{-1} < 1$.

---

## 11. Simulation Details and Convergence

### 11.1 Broadening Parameter $\eta$

The $\eta$ parameter plays two roles:

1. **Physical broadening** — represents inverse quasiparticle lifetime (disorder, phonons). Real materials have $\eta \sim 1$–$100\ \mathrm{meV}$.
2. **Numerical regularisation** — ensures the matrix $zI - H$ is invertible. As $\eta \to 0$, features sharpen to delta functions.

**Resolution vs cost tradeoff:**  Features narrower than $2\eta$ are unresolved. Edge modes in topological systems are intrinsically sharp (no self-energy), so $\eta = 0.02$–$0.05$ gives good resolution. Bulk van Hove singularities are slightly over-broadened.

**Choice in scripts:**

| Script | $\eta$ | Justification |
|---|---|---|
| `09_bulk_spectral.jl` | $0.04$ | Bulk bands, moderate resolution |
| `10_ribbon_bands.jl` | $0.03$ | LDOS, Nx=40 spacing ≈ 0.25 |
| `11_surface_greens.jl` | $0.025$ | Edge modes need sharper features |
| `12_interacting_hf.jl` | $0.05$ | HF bands broader near phase transition |

### 11.2 $\omega$ Grid

The $\omega$ grid should span the full spectrum with enough points to resolve the sharpest feature (width $\sim 2\eta$):

$$N_\omega \gtrsim \frac{\omega_\mathrm{max} - \omega_\mathrm{min}}{2\eta}$$

For $\omega \in [-3.5, 3.5]$ and $\eta = 0.025$:  $N_\omega \gtrsim 7/0.05 = 140$. We use $N_\omega = 250$, oversampling by $\sim 2\times$.

### 11.3 BZ Mesh for DOS

The DOS integral is approximated by a uniform $N_k \times N_k$ mesh. The error in a smooth integrand goes as $O(N_k^{-2})$ (midpoint rule on a torus). Van Hove singularities (band extrema) slow convergence; they appear as shoulders in $\rho(\omega)$ rather than true singularities at finite $\eta$.

Typical choice: $N_k = 100$–$200$ for the DOS, $N_k = 60$ for the Chern number.

### 11.4 López Sancho Convergence

The tolerance for the loop termination condition $\|\alpha_n\| < \varepsilon_\mathrm{tol}$ controls the accuracy of $G_s$:

$$\|G_s^\mathrm{computed} - G_s^\mathrm{exact}\| \lesssim \frac{\|\alpha_n\|}{\eta^2}$$

For $\eta = 0.025$ and $\varepsilon_\mathrm{tol} = 10^{-8}$, the GF error is $\lesssim 10^{-8}/\eta^2 \approx 1.6\times 10^{-5}$, well below the broadening-limited spectral resolution.

### 11.5 Ribbon Size Convergence

The ribbon dispersion $E_n(k_y)$ converges to the semi-infinite result as $N_x \to \infty$. The finite-size hybridisation gap of edge states scales as $e^{-N_x/\xi}$. For the QWZ model near $m = -1$:

$$\xi \approx \frac{v_F}{\Delta} \approx \frac{1}{|m+2|} \approx 1$$

so $N_x = 40$ gives $e^{-40} \approx 10^{-17}$ — completely negligible.

### 11.6 Spectral Sum Rule Check

A quick sanity check for any spectral computation: integrate $A(\mathbf{k},\omega)$ over $\omega$ and verify it approaches 1 (per site). For finite $\eta$ and $\omega \in [-\omega_\mathrm{max}, \omega_\mathrm{max}]$:

$$\int_{-\omega_\mathrm{max}}^{\omega_\mathrm{max}} A(\mathbf{k},\omega)\,d\omega \approx 1 - \frac{2}{\pi}\arctan\!\left(\frac{\eta}{\omega_\mathrm{max}}\right)$$

The correction is $\approx 2\eta/(\pi\omega_\mathrm{max}) \approx 0.004$ for $\eta = 0.025$, $\omega_\mathrm{max} = 3.5$ — negligible.

---

## 12. Physical Quantities Summary

### What to compute and why

#### Bulk spectral function $A(\mathbf{k}, \omega)$

$$A(\mathbf{k},\omega) = -\frac{1}{\pi}\,\mathrm{Im\,Tr}\,G(\mathbf{k},\omega)$$

**Why:** Directly measured by angle-resolved photoemission (ARPES). Reveals the band structure, gap, and van Hove singularities. In an interacting system, it shows quasiparticle renormalisation and damping.

#### Density of states $\rho(\omega)$

$$\rho(\omega) = \frac{1}{N_k^2}\sum_\mathbf{k} A(\mathbf{k},\omega)$$

**Why:** Measured by tunnelling spectroscopy and specific heat. The gap in $\rho(\omega)$ is a bulk diagnostic for the insulating phase. Van Hove singularities appear as peaks.

#### Chern number $C$

$$C = \frac{1}{2\pi}\int_\mathrm{BZ}\Omega(\mathbf{k})\,d^2k \in \mathbb{Z}$$

**Why:** Topological invariant that predicts the number of chiral edge modes (bulk-boundary correspondence). $C \neq 0$ implies a topologically non-trivial phase that cannot be adiabatically connected to a trivial insulator without closing the gap.

#### Ribbon band structure + LDOS

$$\mathrm{LDOS}(j,\omega) = \frac{1}{N_{k_y}}\sum_{k_y} A(j; k_y, \omega)$$

**Why:** Finite-size probe of edge states. In-gap bands visible in the ribbon spectrum confirm topological character. The LDOS shows spatial localisation to $j = 1$ or $j = N_x$, distinguishing edge from bulk states.

#### Surface spectral function $A_s(k_y, \omega)$

$$A_s(k_y,\omega) = -\frac{1}{\pi}\,\mathrm{Im\,Tr}\,G_s(k_y,\omega)$$

**Why:** The semi-infinite geometry eliminates finite-size effects entirely. Edge modes appear as perfectly sharp (width $\sim \eta$) in-gap features, without any hybridisation gap. The surface GF also forms the starting point for multi-layer perturbation theory (e.g.\ adding disorder, magnetic field, or superconducting pairing at the boundary).

---

## References

1. **Qi, Wu, Zhang** (2006) — "Topological quantization of the spin Hall effect in two-dimensional paramagnetic semiconductors." *Phys. Rev. B* 74, 085308. (QWZ model)

2. **López Sancho, López Sancho, Rubio** (1984) — "Quick iterative scheme for the calculation of transfer matrices: application to Mo (100)." *J. Phys. F: Met. Phys.* 14, 1205. (Decimation algorithm)

3. **Fukui, Hatsugai, Suzuki** (2005) — "Chern numbers in discretized Brillouin zone: efficient method of computing (spin) Hall conductances." *J. Phys. Soc. Jpn.* 74, 1674. (FHS lattice Chern formula)

4. **Hasan and Kane** (2010) — "Colloquium: Topological insulators." *Rev. Mod. Phys.* 82, 3045. (Review of topological phases)

5. **Datta** (1995) — *Electronic Transport in Mesoscopic Systems.* Cambridge University Press. (Recursive GF methods)
