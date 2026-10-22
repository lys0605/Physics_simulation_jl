using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, DynamicalBilliards
include(joinpath(@__DIR__, "../../src/nonlinear/billiards.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Poincaré sections: integrable vs chaotic billiards ────────────────────────
#
# The Poincaré (or boundary) map records where the particle hits the boundary
# and at what angle.  DynamicalBilliards.jl's `boundarymap` returns Birkhoff
# canonical coordinates:
#
#   s  ∈ [0, L)   arc-length position along the perimeter
#   p  = sin α    tangential velocity component  ∈ (−1, 1)
#
# These coordinates are area-preserving (ds ∧ dp = const, Liouville's theorem),
# so a random scatter of initial conditions should fill phase space uniformly in
# a chaotic system, while an integrable system shows a foliation of 1-D curves.
#
# Four billiards side-by-side:
#   (a) Square           — integrable; straight walls → shear flow in (s, p)
#   (b) 100-gon ≈ circle — integrable; conservation of angular momentum →
#                          horizontal invariant lines
#   (c) Bunimovich stadium — completely chaotic / ergodic (Bunimovich 1979)
#   (d) Sinai billiard    — completely chaotic / ergodic (Sinai 1970)

const N_ICS  = 30      # initial conditions per billiard
const N_COLS = 800     # collisions per initial condition

# ── billiards ─────────────────────────────────────────────────────────────────
billiards = [
    ("Square",            billiard_rectangle(1.0, 1.0)),
    ("Circle  (100-gon)", make_circle(1.0; n_sides=100)),
    ("Bunimovich stadium", billiard_stadium()),
    ("Sinai billiard",     billiard_sinai(0.25)),
]

# ── compute Poincaré sections ─────────────────────────────────────────────────
function poincare_section(bd, n_ics, n_cols)
    s_all = Float64[];  p_all = Float64[]
    for _ in 1:n_ics
        part = randominside(bd)
        bm, _ = boundarymap(part, bd, n_cols)
        for pt in bm
            push!(s_all, pt[1]);  push!(p_all, pt[2])
        end
    end
    s_all, p_all
end

sections = [poincare_section(bd, N_ICS, N_COLS) for (_, bd) in billiards]

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1400, 680))

for (col, ((label, bd), (s_vals, p_vals))) in enumerate(zip(billiards, sections))
    _, ints = boundarymap(randominside(bd), bd, 1)
    L = ints[end]

    ax = Axis(fig[1, col];
              title   = label,
              xlabel  = col == 1 ? "arc-length  s" : "s",
              ylabel  = col == 1 ? "sin α" : "",
              limits  = (0.0, L, -1.0, 1.0),
              xticklabelsize = 11, yticklabelsize = 11,
              titlesize = 13)

    scatter!(ax, s_vals, p_vals;
             color      = (:navy, 0.25),
             markersize = 1.5,
             strokewidth = 0)

    # Annotate integrable / chaotic
    tag = col ≤ 2 ? "integrable" : "chaotic"
    clr = col ≤ 2 ? :darkgreen   : :darkred
    text!(ax, L*0.5, 0.88; text=tag, align=(:center,:center),
          color=clr, fontsize=13, font=:bold)
end

Label(fig[0, :], "Poincaré sections: integrable vs chaotic billiards";
      fontsize=16, font=:bold)

outpath = joinpath(@__DIR__, "../../output/plots/08_billiards_poincare.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
