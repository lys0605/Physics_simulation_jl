# Physics Simulation

Classical-mechanics simulations written in Julia and rendered straight to
MP4, intended for use as short teaching clips in a classroom setting.

Physics is integrated with
[DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) and
frames are written by [CairoMakie](https://docs.makie.org/)'s `record()`
function, which pipes into `ffmpeg` to produce the video file. The result
is a project where adding a new scene is: write an ODE, write a ~30-line
script, run it, get an MP4.

## Scenes included

| Script | What it shows |
|---|---|
| [`scripts/01_projectile.jl`](scripts/01_projectile.jl) | Ballistic trajectory with quadratic air drag |
| [`scripts/02_pendulum.jl`](scripts/02_pendulum.jl) | Nonlinear pendulum vs small-angle approximation, side by side |
| [`scripts/03_double_pendulum.jl`](scripts/03_double_pendulum.jl) | Chaotic double pendulum with fading trail |
| [`scripts/04_spring_mass.jl`](scripts/04_spring_mass.jl) | Damped harmonic oscillator + live phase-space plot |

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
julia --project=. scripts/01_projectile.jl
```

The MP4 appears in `output/videos/01_projectile.mp4`. Same pattern for the
other scripts.

## Project layout

```
src/mechanics/   reusable physics models — pure ODE right-hand-sides
                 plus simulate_* helpers that return an ODESolution
scripts/         one file per scene; each produces one MP4
output/videos/   rendered MP4 output (gitignored)
Project.toml     pinned Julia environment
```

The split is deliberate: physics lives in `src/mechanics/` and knows
nothing about plotting, while `scripts/` files know nothing about how the
ODE is integrated. The same mechanics file can be reused across multiple
scenes — for example, the pendulum model is used to visualise both the
raw motion and the small-angle comparison.

## How a scene is built

Every script follows the same skeleton:

```julia
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using CairoMakie
include(joinpath(@__DIR__, "..", "src", "mechanics", "projectile.jl"))

sol = simulate_projectile(v0=30.0, angle_deg=55.0, drag=0.02,
                          tspan=(0.0, 6.0))

ts  = range(sol.t[1], sol.t[end]; length=240)  # 240 frames @ 30 fps = 8 s

fig = Figure(size=(960, 540))
ax  = Axis(fig[1, 1]; xlabel="x (m)", ylabel="y (m)")
pt    = Observable(Point2f(0, 0))
trail = Observable(Point2f[])
lines!(ax, trail); scatter!(ax, pt)

record(fig, "output/videos/01_projectile.mp4", eachindex(ts); framerate=30) do i
    u = sol(ts[i])
    pt[] = Point2f(u[1], u[2])
    push!(trail[], pt[]); notify(trail)
end
```

Two conventions keep the animations smooth and reproducible:

1. **Frame timing is decoupled from ODE steps.** The script samples the
   solution with `sol(t)` at evenly-spaced times — framerate is a render
   knob, not a solver knob, so stiff systems still animate smoothly.
2. **`Figure` is built once, then `Observable`s are mutated inside the
   `record()` loop.** Nothing is rebuilt per frame.

## Adding a new scene

1. Create `src/mechanics/<name>.jl` exposing a `simulate_<name>` helper
   that returns an `ODESolution`. Keep it pure — no plotting, no I/O.
2. Copy any script in `scripts/` as a template, change the `include`,
   `simulate_` call, axes, and `record` loop body.
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
