# Physics Simulation

Classical-mechanics simulations written in Julia and rendered straight to
MP4, intended for use as short teaching clips in a classroom setting.

Physics is integrated with
[DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) and
frames are written by [CairoMakie](https://docs.makie.org/)'s `record()`
function, which pipes into `ffmpeg` to produce the video file. The result
is a project where adding a new scene is: write an ODE, write a ~30-line
script, run it, get an MP4.

The **nonlinear dynamics** module uses
[DynamicalBilliards.jl](https://juliadynamics.github.io/DynamicalBilliards.jl/stable/)
for collision geometry; no ODE solver is needed because the particle travels
in a straight line between bounces.

## Scenes included

### Geometry — minimal surfaces

| Script | What it shows | Output |
|---|---|---|
| [`scripts/geometry/13_catenoid.jl`](scripts/geometry/13_catenoid.jl) | Catenoid as a soap-film minimal surface (H = 0): two equal rings stretch apart; the neck radius shrinks self-consistently until the critical separation h/R ≈ 0.6627 is reached, at which point the film snaps into two flat discs | MP4 |

### Moved to Physics Atlas

These scenes now live, rebuilt and checked, in [lys0605/physics-atlas](https://github.com/lys0605/physics-atlas). The old versions are in this repository's history (commits 8397691, 1dce21e and f70930e).

| Old | Now |
|---|---|
| `scripts/02_pendulum.jl`, `src/mechanics/pendulum.jl` | `mec.pendulum-large-angle` |
| `scripts/01_projectile.jl`, `src/mechanics/projectile.jl` | `mec.projectile-drag` |
| `scripts/03_double_pendulum.jl`, `src/mechanics/double_pendulum.jl` | `nld.double-pendulum` |
| `scripts/04_spring_mass.jl`, `src/mechanics/spring_mass.jl` | `mec.driven-oscillator` |
| `scripts/05_cyclic_pursuit.jl`, `src/mechanics/cyclic_pursuit.jl` | `mec.cyclic-pursuit` |
| `scripts/lattice/lattices_3d.jl` | `cmp.brillouin-zones` (conventional and primitive cells of all 14 lattices) |
| `scripts/nonlinear/06–08_billiards_*.jl`, `src/nonlinear/billiards.jl`, `docs/nonlinear/billiards_theory.md` | `nld.billiards` |
| `scripts/em/14_fresnel.jl` | `opt.fresnel` |
| `docs/lattice/interactive_viewer.html` | `cmp.brillouin-zones` |
| `scripts/lattice/lattices_2d.jl` | `cmp.lattices-2d` (all five 2D lattices, cells, zones, lattice lines) |
| `scripts/lattice/brillouin.jl` | `cmp.lattices-2d` (2D) and `cmp.brillouin-zones` (3D) |

### Condensed-matter — Green's functions & topological systems

Built on the [Qi-Wu-Zhang (QWZ)](https://link.aps.org/doi/10.1103/PhysRevB.74.085308) 2-band Chern insulator on the square lattice. Three boundary conditions are explored in order of increasing sophistication.

| Script | Boundary condition | What it shows | Output |
|---|---|---|---|
| [`scripts/greens/09_bulk_spectral.jl`](scripts/greens/09_bulk_spectral.jl) | Periodic in x and y | Spectral function A(k,ω) along Γ→X→M→Γ, density of states ρ(ω), Chern number C | PNG |
| [`scripts/greens/10_ribbon_bands.jl`](scripts/greens/10_ribbon_bands.jl) | Periodic in y, open in x (ribbon, Nx=40) | Band structure E(ky) with topological edge states highlighted; layer-resolved DOS confirming edge localisation | PNG |
| [`scripts/greens/11_surface_greens.jl`](scripts/greens/11_surface_greens.jl) | Periodic in y, semi-infinite in x | Surface spectral function via López Sancho decimation; convergence diagnostics; comparison with naive fixed-point iteration | PNG |
| [`scripts/greens/12_interacting_hf.jl`](scripts/greens/12_interacting_hf.jl) | Periodic (bulk), Hartree-Fock interaction U n₁n₂ | Self-consistent HF renormalization; orbital polarization; interaction-driven topological phase transition at U_c ≈ 2.46 | PNG |

Theory and mathematical background: [`docs/greens/greens_theory.md`](docs/greens/greens_theory.md)

### Electromagnetism & optics

Shared utilities live in `src/em/em_utils.jl` (Fresnel equations, point-charge
fields, Biot-Savart, Hertzian dipole retarded fields).

| Script | What it shows | Output |
|---|---|---|
| [`scripts/em/15_multipoles.jl`](scripts/em/15_multipoles.jl) | Electrostatic potential heat maps and field-line streamplots for monopole, dipole, quadrupole, and octupole charge configurations | PNG |
| [`scripts/em/16_magnetic_fields.jl`](scripts/em/16_magnetic_fields.jl) | Biot-Savart field maps (log-magnitude + streamlines) for a magnetic dipole, a single circular loop, a Helmholtz coil (uniform-field region), and an anti-Helmholtz / Maxwell coil (linear gradient) | PNG |
| [`scripts/em/17_radiation.jl`](scripts/em/17_radiation.jl) | Animated E_z snapshot of a Hertzian dipole radiating into the xz-plane over one full oscillation period; side panel shows the sin²θ polar power pattern | MP4 |

### Crystallography — lattices, Voronoi cells, Wigner-Seitz, Brillouin zones

Shared utilities: `src/lattice/lattice_utils.jl` (Bravais lattices, reciprocal
basis), `src/lattice/voronoi.jl` (2-D Sutherland-Hodgman clipping, 3-D
vertex enumeration for half-space intersection).

| Script | What it shows | Output |
|---|---|---|
| [`scripts/lattice/voronoi.jl`](scripts/lattice/voronoi.jl) | Voronoi tessellation of random points; specialisation to lattices showing the Wigner-Seitz cell as the Voronoi cell of a Bravais lattice; explicit bisector-construction zoom | PNG |

Theory and mathematical background: [`docs/lattice/lattice_theory.md`](docs/lattice/lattice_theory.md)

**Interactive viewer:** moved to Physics Atlas (`cmp.brillouin-zones`).

### Monte Carlo and Quantum Monte Carlo

Notes and companion scripts.  The notes can be read independently as a
syllabus; the scripts are short (≤ 100 lines) and meant to be modified.

**Notes** (`docs/monte_carlo/`):

| Note | Scope |
|---|---|
| [`README.md`](docs/monte_carlo/README.md) | Index and motivation |
| [`01_theory.md`](docs/monte_carlo/01_theory.md) | Classical MC theory — importance sampling, Markov chains, Metropolis-Hastings, detailed balance, cluster algorithms |
| [`02_numerical.md`](docs/monte_carlo/02_numerical.md) | Numerical practice — RNGs, equilibration, autocorrelation, binning/jackknife/bootstrap, parallel tempering, Wang-Landau |
| [`03_quantum_monte_carlo.md`](docs/monte_carlo/03_quantum_monte_carlo.md) | VMC, DMC, PIMC, SSE, DQMC/AFQMC, CT-QMC, and the sign problem |
| [`04_learning_path.md`](docs/monte_carlo/04_learning_path.md) | Step-by-step learning outline from π-estimation to AFQMC |
| [`05_related_methods.md`](docs/monte_carlo/05_related_methods.md) | Brief survey of ED, DMRG/MPS, PEPS, DMFT, DFT, HF/CC, fRG |
| [`06_cross_field_examples.md`](docs/monte_carlo/06_cross_field_examples.md) | Famous MC models in sociology, biology, economics, chemistry, ML, ecology, finance |
| [`07_scripts_walkthrough.md`](docs/monte_carlo/07_scripts_walkthrough.md) | Per-script guide tying each Julia file to the equations in notes 01-06, with exercises |

**Companion scripts** (`scripts/monte_carlo/`) — each is self-contained and
runnable with `julia --project=. scripts/monte_carlo/<name>.jl`:

| Script | What it shows | Output |
|---|---|---|
| [`pi_estimation.jl`](scripts/monte_carlo/pi_estimation.jl) | π by hit-or-miss sampling; demonstrates the 1/√N convergence law | PNG |
| [`random_walk.jl`](scripts/monte_carlo/random_walk.jl) | 2-D random walk trajectories; mean-square displacement ∼ t (diffusive scaling) | PNG |
| [`ising_metropolis.jl`](scripts/monte_carlo/ising_metropolis.jl) | 2-D Ising Metropolis sweep over T: E, ⟨\|m\|⟩, χ, C_V, snapshots at T = 1.8 / T_c / 3.0 | PNG |
| [`ising_wolff.jl`](scripts/monte_carlo/ising_wolff.jl) | Wolff cluster vs Metropolis: dramatic τ_int reduction near T_c | PNG |
| [`autocorrelation.jl`](scripts/monte_carlo/autocorrelation.jl) | Autocorrelation function, integrated time, block-binning error plateau at three temperatures | PNG |
| [`vmc_harmonic.jl`](scripts/monte_carlo/vmc_harmonic.jl) | VMC for the 1-D harmonic oscillator; energy & variance vs α, zero-variance principle, stochastic-gradient optimisation | PNG |
| [`percolation.jl`](scripts/monte_carlo/percolation.jl) | 2-D site percolation; largest cluster fraction vs p for L = 32, 64, 128; snapshots around p_c ≈ 0.5927 | PNG |
| [`schelling.jl`](scripts/monte_carlo/schelling.jl) | Schelling segregation: residential dynamics evolve from random to segregated; global segregation index in side panel | MP4 |
| [`wright_fisher.jl`](scripts/monte_carlo/wright_fisher.jl) | Population-genetics drift; neutral vs selected trajectories; fixation probability vs initial frequency (Kimura) | PNG |
| [`forest_fire.jl`](scripts/monte_carlo/forest_fire.jl) | Drossel-Schwabl forest-fire automaton — self-organised criticality, animated state + cluster-size power law | MP4 |
| [`gillespie.jl`](scripts/monte_carlo/gillespie.jl) | Exact stochastic simulation algorithm on Lotka-Volterra; single trajectory, phase portrait, ensemble | PNG |

Shared utilities in `src/monte_carlo/mc_utils.jl`: autocorrelation, integrated
autocorrelation time, jackknife, block-binning error curves.

## Requirements

- [Julia](https://julialang.org/downloads/) ≥ 1.10 (tested on 1.11)
- [`ffmpeg`](https://ffmpeg.org/) on your `PATH` (for MP4 encoding)

No display server is required — CairoMakie renders headlessly.

## Setup

```bash
git clone <this-repo-url>
cd Physics_simulation
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The first `instantiate` takes several minutes because it precompiles
DifferentialEquations and CairoMakie. Subsequent launches are fast.

## Running a scene

```bash
julia --project=. scripts/geometry/13_catenoid.jl
```

The MP4 appears in `output/videos/13_catenoid.mp4`. Same pattern for the
other scripts.

## Project layout

```
src/greens/      Green's function / condensed-matter models
src/em/          electromagnetic utilities (Fresnel, multipoles, Biot-Savart, dipole)
src/lattice/     Bravais lattices, reciprocal lattice, Voronoi / WS / BZ algorithms
src/monte_carlo/ MC utilities (autocorrelation, jackknife, block-binning errors)
src/plot_utils.jl  shared plotting helpers (see below)
scripts/         one file per scene; each produces one MP4 or PNG
output/videos/   rendered MP4 output (gitignored)
output/plots/    rendered PNG output (gitignored)
Project.toml     pinned Julia environment
```

The split is deliberate: physics lives in `src/<area>/` and knows
nothing about plotting, while `scripts/` files know nothing about how the
model is computed. The same physics file can be reused across several
scenes — for example, `src/greens/qwz_model.jl` serves all four
Green's-function scripts.

## Plotting utilities (`src/plot_utils.jl`)

Common plot types are wrapped in four helpers so scripts stay short.
Include the file and call the function — no boilerplate needed.

```julia
include(joinpath(@__DIR__, "..", "src", "plot_utils.jl"))
```

### `animate_trail(positions, outpath; ...)`

Animates a particle moving through a `Vector{Point2f}` with an
accumulating trail and saves an MP4.

```julia
ts        = range(sol.t[1], sol.t[end]; length=240)
positions = [Point2f(sol(t)[1], sol(t)[2]) for t in ts]
animate_trail(positions, "output/videos/01_projectile.mp4";
              title="Projectile motion", xlabel="x (m)", ylabel="y (m)",
              ylims=(0.0, 50.0))
```

Key keyword arguments: `framerate`, `trail_color`, `particle_color`,
`markersize`, `padding`, `xlims`, `ylims`, `figsize`.

### `savefig(path, fig_or_result; ...)`

Convenience wrapper around CairoMakie's `save` that also creates any missing
parent directories. Accepts a `Figure` or the `(; fig, ax)` named tuple
returned by any `plot_*` function, so you can pass the result directly.

```julia
savefig("output/plots/sin.png", plot_function1d(sin, xs))

# with higher resolution for PNG
fig, ax = plot_heatmap2d(xs, ys, zs)
savefig("output/plots/heatmap.png", fig; px_per_unit=2)

# vector format
savefig("output/plots/curve.pdf", plot_curve2d(xs, ys))
```

Key keyword arguments (forwarded to CairoMakie's `save`):
- `px_per_unit = 2`  raster resolution multiplier (PNG)
- `pt_per_unit = 1`  point size for vector formats (PDF, SVG)

---

All static plot functions return a named tuple `(; fig, ax)`. You can unpack
it in two ways and then keep modifying the plot:

```julia
# positional unpacking
fig, ax = plot_function1d(sin, xs)

# named unpacking (more explicit)
(; fig, ax) = plot_function1d(sin, xs)

# add a second curve to the same axes
lines!(ax, xs, cos.(xs); color=:crimson, label="cos")
axislegend(ax)
save("output/sincos.png", fig)
```

### `plot_heatmap2d(xs, ys, zs; ...)` → `(; fig, ax)`

2-D scalar field rendered as a colour-mapped heat map.
`zs[i, j]` is the value at `(xs[i], ys[j])`.

```julia
xs = range(-2, 2; length=200)
ys = range(-2, 2; length=200)
zs = [exp(-(x^2 + y^2)) for x in xs, y in ys]
fig, ax = plot_heatmap2d(xs, ys, zs; title="Gaussian", clabel="intensity")
contour!(ax, xs, ys, zs; color=:white, levels=5)   # overlay contour lines
save("output/gaussian.png", fig)
```

Key keyword arguments: `colormap`, `clabel`, `figsize`.

### `plot_surface3d(xs, ys, zs; ...)` → `(; fig, ax)`

3-D surface plot (`ax` is an `Axis3`).
`zs[i, j]` is the height at `(xs[i], ys[j])`.

```julia
xs  = range(-π, π; length=100)
ys  = range(-π, π; length=100)
zs  = [sin(x) * cos(y) for x in xs, y in ys]
fig, ax = plot_surface3d(xs, ys, zs; title="sin(x)cos(y)", zlabel="z")
ax.elevation = 0.4   # adjust camera angle after the fact
save("output/surface.png", fig)
```

Key keyword arguments: `colormap`, `zlabel`, `figsize`.

### `plot_vectorfield2d(xs, ys, us, vs; ...)` → `(; fig, ax)`

2-D vector field drawn with arrows.
`us[i, j]` / `vs[i, j]` are the x/y components at `(xs[i], ys[j])`.
Arrows are auto-scaled so the longest fits within one grid cell
(`normalize=true` by default).

```julia
xs = range(-2, 2; length=20)
ys = range(-2, 2; length=20)
us = [-y for x in xs, y in ys]
vs = [ x for x in xs, y in ys]
fig, ax = plot_vectorfield2d(xs, ys, us, vs; title="Circular flow")
save("output/vectorfield.png", fig)
```

Key keyword arguments: `arrowcolor`, `arrowsize`, `normalize`, `figsize`.

### `plot_function1d(xs, ys; ...)` / `plot_function1d(f, xs; ...)` → `(; fig, ax)`

Plot y = f(x) as a line. Pass either pre-computed arrays or a callable.

```julia
xs  = range(0, 2π; length=500)
fig, ax = plot_function1d(sin, xs; title="sin(x)", xlabel="x", ylabel="sin(x)")
save("output/sin.png", fig)

# overlay a second curve
ys = @. exp(-xs/4) * sin(4xs)
fig, ax = plot_function1d(xs, ys; title="Damped oscillation")
hlines!(ax, [0.0]; color=:gray, linestyle=:dash)
save("output/damped.png", fig)
```

Key keyword arguments: `color`, `linewidth`, `figsize`.

### `plot_polar(θs, rs; ...)` / `plot_polar(f, θs; ...)` → `(; fig, ax)`

Plot a polar curve r = f(θ). Cartesian conversion is handled internally;
the axis is always square so circles look circular.

```julia
θs  = range(0, 2π; length=500)
fig, ax = plot_polar(θ -> 1 + cos(θ), θs; title="Cardioid")
save("output/cardioid.png", fig)

# overlay two polar curves on the same axes
fig, ax = plot_polar(θ -> cos(3θ), θs; title="Rose curves")
lines!(ax, sin.(3θs) .* cos.(θs), sin.(3θs) .* sin.(θs); color=:crimson)
save("output/roses.png", fig)
```

Key keyword arguments: `color`, `linewidth`, `figsize`.

### `plot_curve2d(xs, ys; ...)` → `(; fig, ax)`

Plot a parametric plane curve (x(t), y(t)). Unlike `plot_vectorfield2d` this
draws only the path — no vectors at grid points. The axis is square by default.
Pass `color_values` to colour the curve by a scalar (e.g. speed or time).

```julia
# Lissajous figure
ts  = range(0, 2π; length=1000)
fig, ax = plot_curve2d(sin.(3ts), sin.(2ts); title="Lissajous 3:2")
save("output/lissajous.png", fig)

# pendulum phase portrait coloured by speed
fig, ax = plot_curve2d(θs, ωs; title="Phase portrait",
                       xlabel="θ", ylabel="ω",
                       color_values=abs.(ωs), clabel="|ω|",
                       aspect=nothing)
save("output/phase.png", fig)
```

Key keyword arguments: `color`, `color_values`, `colormap`, `clabel`,
`linewidth`, `aspect`, `figsize`.

## How a scene is built

Every script follows the same skeleton (shown with the projectile scene, which has since moved to Physics Atlas):

```julia
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
include(joinpath(@__DIR__, "..", "src", "mechanics", "projectile.jl"))
include(joinpath(@__DIR__, "..", "src", "plot_utils.jl"))

sol = simulate_projectile(v0=30.0, angle_deg=55.0, drag=0.02,
                          tspan=(0.0, 6.0))

ts        = range(sol.t[1], sol.t[end]; length=240)  # 240 frames @ 30 fps = 8 s
positions = [Point2f(sol(t)[1], sol(t)[2]) for t in ts]

animate_trail(positions, "output/videos/01_projectile.mp4";
              title="Projectile motion with quadratic drag",
              xlabel="x (m)", ylabel="y (m)")
```

Two conventions keep the animations smooth and reproducible:

1. **Frame timing is decoupled from ODE steps.** Scripts sample the
   solution with `sol(t)` at evenly-spaced times — framerate is a render
   knob, not a solver knob, so stiff systems still animate smoothly.
2. **`Figure` is built once, then `Observable`s are mutated inside the
   `record()` loop.** Nothing is rebuilt per frame.

## Adding a new scene

1. Create `src/<area>/<name>.jl` exposing a `simulate_<name>` helper
   that returns an `ODESolution`. Keep it pure — no plotting, no I/O.
2. Copy any script in `scripts/` as a template, swap the `include` and
   `simulate_` call, then use a helper from `plot_utils.jl` to render.
3. `julia --project=. scripts/<your_script>.jl`

## Optional: interactive / 3D

Everything here uses CairoMakie (headless, 2D, clean vector output). To
add interactive or 3D rendering:

```bash
julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'
```

Then swap `using CairoMakie` for `using GLMakie` in the script — the
Makie plotting API is identical across backends.

## License

MIT — see [LICENSE](LICENSE) if added.
