using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, Random, Statistics
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Schelling segregation model (1971) ───────────────────────────────────────
#
# L × L grid populated by two groups (red and blue) and some empty cells.
# Each agent is "happy" if at least a fraction τ of its (≤ 8) neighbours
# share its colour.  Unhappy agents move to a random empty cell.
#
# Even at modest τ ≈ 0.3, segregated clusters emerge — the canonical
# illustration of how micromotives produce macro-patterns.

Random.seed!(2026)

const L = 60
const FRAC_RED   = 0.45
const FRAC_BLUE  = 0.45
const TOLERANCE  = 0.35           # τ — desired same-colour fraction
const N_STEPS    = 240
const FRAMERATE  = 30

# ── Initialise grid ──────────────────────────────────────────────────────────
function init_grid(L, fr, fb)
    grid = zeros(Int8, L, L)        # 0 empty, +1 red, −1 blue
    N    = L * L
    nr   = round(Int, fr * N);  nb = round(Int, fb * N)
    cells = shuffle(collect(1:N))
    for k in 1:nr;   grid[cells[k]]            =  1; end
    for k in 1:nb;   grid[cells[nr + k]]       = -1; end
    grid
end

# ── Single update ────────────────────────────────────────────────────────────
@inline function same_color_fraction(grid, i, j, L)
    color = grid[i, j];   color == 0 && return 1.0
    same = 0;  total = 0
    for di in -1:1, dj in -1:1
        (di == 0 && dj == 0) && continue
        ni, nj = mod1(i+di, L), mod1(j+dj, L)
        c = grid[ni, nj]
        if c != 0
            total += 1
            same  += (c == color)
        end
    end
    total == 0 ? 1.0 : same / total
end

function segregation_step!(grid, L; τ=TOLERANCE)
    empties = [(i, j) for i in 1:L, j in 1:L if grid[i, j] == 0]
    occupants = [(i, j) for i in 1:L, j in 1:L if grid[i, j] != 0]
    shuffle!(occupants)
    moves = 0
    for (i, j) in occupants
        if same_color_fraction(grid, i, j, L) < τ && !isempty(empties)
            target_idx = rand(1:length(empties))
            ei, ej = empties[target_idx]
            grid[ei, ej] = grid[i, j]
            grid[i, j]   = 0
            empties[target_idx] = (i, j)
            moves += 1
        end
    end
    moves
end

function global_segregation(grid, L)
    s = 0.0;  count = 0
    for i in 1:L, j in 1:L
        if grid[i, j] != 0
            s += same_color_fraction(grid, i, j, L)
            count += 1
        end
    end
    s / max(count, 1)
end

# ── Animate ──────────────────────────────────────────────────────────────────
grid = init_grid(L, FRAC_RED, FRAC_BLUE)

set_theme!(theme_light())
fig = Figure(size=(900, 480))
Label(fig[0, :], "Schelling Segregation Model  —  tolerance τ = $TOLERANCE";
      fontsize=14, font=:bold)

ax = Axis(fig[1, 1]; aspect=DataAspect(),
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false,
          title="L = $L,  red = $(round(Int, 100*FRAC_RED))%,  blue = $(round(Int, 100*FRAC_BLUE))%")
ax_seg = Axis(fig[1, 2]; xlabel="step", ylabel="mean same-colour fraction",
              title="global segregation", limits=(0, N_STEPS, 0.4, 1.0))

grid_obs = Observable(copy(grid))
seg_hist = Observable(Float64[global_segregation(grid, L)])
step_hist = Observable(Float64[0.0])

heatmap!(ax, 1:L, 1:L, grid_obs;
         colormap=cgrad([:dodgerblue, :white, :crimson]),
         colorrange=(-1, 1))
lines!(ax_seg, step_hist, seg_hist; color=:purple, linewidth=2)
hlines!(ax_seg, [TOLERANCE]; color=:gray, linestyle=:dash)

mkpath(joinpath(@__DIR__, "../../output/videos"))
out = joinpath(@__DIR__, "../../output/videos/schelling.mp4")

record(fig, out, 1:N_STEPS; framerate=FRAMERATE) do step
    moves = segregation_step!(grid, L)
    grid_obs[] = copy(grid)
    push!(step_hist[], Float64(step));  step_hist[] = step_hist[]
    push!(seg_hist[],  global_segregation(grid, L));  seg_hist[]  = seg_hist[]
end

@info "wrote $out"
