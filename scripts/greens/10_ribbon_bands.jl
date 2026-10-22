using Pkg; Pkg.activate(joinpath(@__DIR__, "../.."))
using CairoMakie, LinearAlgebra
include(joinpath(@__DIR__, "../../src/greens/qwz_model.jl"))
include(joinpath(@__DIR__, "../../src/plot_utils.jl"))

# ── Ribbon geometry: periodic BC in y, open BC in x ───────────────────────────
#
# The strip has Nx unit cells in x (each with 2 orbitals), so the Hamiltonian
# at each k_y is a 2Nx × 2Nx Hermitian matrix.  Diagonalising it gives 2Nx
# bands E_n(k_y).  For a topological system (|C| = 1), in-gap edge bands
# cross the bulk gap and connect the valence and conduction bands.
#
# Three panels per phase (trivial | topological):
#   Left  : E(k_y) band structure — edge states visible as lines in the gap
#   Right : Layer-resolved DOS  LDOS(j, ω) at k_y = π/2
#           bright at j = 1 or Nx → edge localisation

const Nx  = 40           # unit cells in x (80 orbitals per k_y)
const η   = 0.03
const Nky = 300
const kys = range(-π, π; length = Nky)

models = [QWZModel(3.0), QWZModel(-1.0)]
labels = ["Trivial  (m = 3, C = 0)", "Topological  (m = −1, C = +1)"]
clrs   = [:gray70, :dodgerblue]

# ── ribbon spectrum ────────────────────────────────────────────────────────────
function ribbon_spectrum(model, kys, Nx)
    [real.(eigvals(ribbon_hamiltonian(model, ky, Nx))) for ky in kys]
end

# ── layer-resolved spectral weight at a single k_y ────────────────────────────
# Uses the eigendecomposition:  LDOS(j,ω) = Σ_n |⟨j|n⟩|² η/[(ω−Eₙ)²+η²] / π
function layer_ldos_slice(model, ky, ωs, Nx; η=0.03)
    H  = ribbon_hamiltonian(model, ky, Nx)
    Es = real.(eigvals(H))
    Vs = eigvecs(H)                   # columns = eigenvectors
    A  = zeros(Float64, Nx, length(ωs))
    for j in 1:Nx
        for n in 1:2Nx
            w = abs2(Vs[2j-1, n]) + abs2(Vs[2j, n])   # orbital weight on layer j
            for (iω, ω) in enumerate(ωs)
                A[j, iω] += w * η / ((ω - Es[n])^2 + η^2)
            end
        end
    end
    A ./ π
end

@info "Computing ribbon spectra and LDOS…"
spectra = [ribbon_spectrum(m, kys, Nx) for m in models]
ωs_ldos = range(-3.5, 3.5; length = 200)
ky_sel  = π / 2          # k_y at which to show the LDOS
ldos    = [layer_ldos_slice(m, ky_sel, ωs_ldos, Nx; η) for m in models]

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1300, 600))

for (col, (lbl, spec, A)) in enumerate(zip(labels, spectra, ldos))
    # Band structure
    ax_b = Axis(fig[1, 2col-1];
                title      = lbl,
                xlabel     = "k_y",
                ylabel     = col == 1 ? "E" : "",
                limits     = (-π, π, -3.5, 3.5),
                xticks     = ([-π, 0, π], ["-π", "0", "π"]),
                titlesize  = 13)

    # Flatten: for each k-point i, repeat kys[i] for all 2Nx eigenvalues
    all_ks = [kys[i] for (i, E_vec) in enumerate(spec) for _ in E_vec]
    all_Es = reduce(vcat, spec)
    scatter!(ax_b, all_ks, all_Es; color = (:black, 0.15), markersize = 1.5)

    # Highlight in-gap states (|E| < 1.2 for m=-1 gap ≈ 2)
    if col == 2
        ks_edge = [kys[i] for (i, E_vec) in enumerate(spec) for E in E_vec if abs(E) < 1.2]
        Es_edge = [E       for E_vec      in spec            for E in E_vec if abs(E) < 1.2]
        scatter!(ax_b, ks_edge, Es_edge; color = :crimson, markersize = 2.5)
    end
    hlines!(ax_b, [0.0]; color = :gray, linewidth = 0.7, linestyle = :dot)

    # LDOS heatmap
    ax_l = Axis(fig[1, 2col];
                title    = "LDOS(layer, ω)  at  k_y = π/2",
                xlabel   = "layer j",
                ylabel   = col == 1 ? "ω" : "",
                limits   = (1, Nx, -3.5, 3.5),
                titlesize = 12,
                yticklabelsvisible = col == 1)

    heatmap!(ax_l, 1:Nx, ωs_ldos, A;
             colormap = :inferno, colorrange = (0, 0.6))
end

Label(fig[0, :], "Ribbon band structure and layer-resolved DOS — open BC in x, periodic in y";
      fontsize = 14, font = :bold)

outpath = joinpath(@__DIR__, "../../output/plots/10_ribbon_bands.png")
savefig(outpath, fig; px_per_unit = 2)
@info "wrote $outpath"
