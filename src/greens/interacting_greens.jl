# Hartree-Fock Green's function for the QWZ model.
#
# Interaction:  H_int = U Σ_R  n_{R,1} n_{R,2}
#
# For a contact (k-independent) U, the HF self-energy is also k-independent:
#
#   Σ^HF = [[ U ρ₂₂,   −U ρ₁₂ ],
#            [−U ρ₂₁,    U ρ₁₁ ]]
#
# Diagonal (Hartree): orbital a feels the average density of orbital b.
# Off-diagonal (Fock / exchange): mixes the two orbitals.
#
# Self-consistency: ρ_{ab} = (1/N_k²) Σ_k ⟨a|u₋(k)⟩⟨u₋(k)|b⟩
# where |u₋(k)⟩ is the lower eigenstate of H_HF(k) = H₀(k) + Σ^HF.
#
# Note: Σ^HF is real/static → quasiparticle weight Z = 1 exactly.
# For Z < 1 (finite lifetime) one needs the second-order Born approximation.
#
# Depends on _link defined in bulk_greens.jl (include that file first).

using LinearAlgebra

"""
    density_matrix_hf(model, Σ_HF; Nk=60) → 2×2 ComplexF64

Ground-state density matrix at half-filling (lower band occupied):
  ρ_{ab} = (1/N_k²) Σ_k ⟨a|u₋(k)⟩⟨u₋(k)|b⟩
"""
function density_matrix_hf(model::QWZModel, Σ_HF::AbstractMatrix; Nk::Int=60)
    ks = range(0, 2π; length=Nk+1)[1:Nk]
    ρ  = zeros(ComplexF64, 2, 2)
    for kx in ks, ky in ks
        u = eigvecs(bloch_hamiltonian(model, kx, ky) + Σ_HF)[:, 1]
        ρ .+= u * u'
    end
    ρ / Nk^2
end

"""
    hf_selfenergy(U, ρ) → 2×2 ComplexF64

Hartree-Fock self-energy for H_int = U n₁ n₂.
"""
function hf_selfenergy(U::Real, ρ::AbstractMatrix)
    ComplexF64[ U*ρ[2,2]  -U*ρ[1,2];
               -U*ρ[2,1]   U*ρ[1,1] ]
end

"""
    hf_selfconsistent(model, U; Nk, max_iter, tol, mixing) → (Σ, ρ)

Self-consistent HF loop with linear mixing α:  Σ ← α Σ_new + (1−α) Σ_old.
Returns converged (Σ^HF, ρ).
"""
function hf_selfconsistent(model::QWZModel, U::Real;
                            Nk::Int       = 60,
                            max_iter::Int = 400,
                            tol::Real     = 1e-9,
                            mixing::Real  = 0.4)
    Σ = zeros(ComplexF64, 2, 2)
    for _ in 1:max_iter
        ρ     = density_matrix_hf(model, Σ; Nk)
        Σ_new = hf_selfenergy(U, ρ)
        Σ_mix = mixing * Σ_new + (1 - mixing) * Σ
        if norm(Σ_mix - Σ) < tol
            return Σ_mix, density_matrix_hf(model, Σ_mix; Nk)
        end
        Σ = Σ_mix
    end
    @warn "HF did not converge (U=$U)"
    Σ, density_matrix_hf(model, Σ; Nk)
end

"""
    hf_convergence_history(model, U; Nk, max_iter, tol, mixing)
        → (Σ, ρ, history)

Same as hf_selfconsistent but also returns the ‖ΔΣ‖ per iteration.
"""
function hf_convergence_history(model::QWZModel, U::Real;
                                  Nk::Int=60, max_iter::Int=400,
                                  tol::Real=1e-9, mixing::Real=0.4)
    Σ   = zeros(ComplexF64, 2, 2)
    hist = Float64[]
    for _ in 1:max_iter
        ρ     = density_matrix_hf(model, Σ; Nk)
        Σ_new = hf_selfenergy(U, ρ)
        Σ_mix = mixing * Σ_new + (1 - mixing) * Σ
        push!(hist, norm(Σ_mix - Σ))
        if hist[end] < tol
            return Σ_mix, density_matrix_hf(model, Σ_mix; Nk), hist
        end
        Σ = Σ_mix
    end
    @warn "HF did not converge (U=$U)"
    Σ, density_matrix_hf(model, Σ; Nk), hist
end

"""
    chern_number_hf(model, Σ_HF; Nk=60) → Int

Chern number of the effective HF Hamiltonian H₀(k) + Σ^HF (Fukui-Hatsugai-Suzuki).
Requires _link to be defined (from bulk_greens.jl).
"""
function chern_number_hf(model::QWZModel, Σ_HF::AbstractMatrix; Nk::Int=60)
    ks = range(0, 2π; length=Nk+1)[1:Nk]
    u  = [ eigvecs(bloch_hamiltonian(model, kx, ky) + Σ_HF)[:, 1]
           for kx in ks, ky in ks ]
    C  = 0.0
    for i in 1:Nk, j in 1:Nk
        ip = mod1(i+1, Nk);  jp = mod1(j+1, Nk)
        F  = _link(u[i,j], u[ip,j]) * _link(u[ip,j], u[ip,jp]) *
             conj(_link(u[i,jp], u[ip,jp])) * conj(_link(u[i,j], u[i,jp]))
        C += angle(F)
    end
    round(Int, C / (2π))
end

"""
    interacting_spectral_map(model, Σ_HF, kxs, kys_path, ωs; η)
        → Matrix{Float64}

Spectral function A_HF(k,ω) along a k-path using G_HF = [(ω+iη)I − H₀(k) − Σ^HF]⁻¹.
"""
function interacting_spectral_map(model::QWZModel, Σ_HF::AbstractMatrix,
                                   kxs, kys_path, ωs; η::Real=0.04)
    Nk = length(kxs);  Nω = length(ωs)
    A  = zeros(Float64, Nk, Nω)
    for (i, (kx, ky)) in enumerate(zip(kxs, kys_path))
        H_eff = bloch_hamiltonian(model, kx, ky) + Σ_HF
        for (j, ω) in enumerate(ωs)
            G = inv(complex(ω, η) * I - H_eff)
            A[i, j] = -imag(tr(G)) / π
        end
    end
    A
end

"""
    dos_hf(model, Σ_HF, ωs; Nk=80, η=0.04) → Vector{Float64}

Density of states with HF renormalization: ρ_HF(ω) = (1/N_k²) Σ_k A_HF(k,ω).
"""
function dos_hf(model::QWZModel, Σ_HF::AbstractMatrix, ωs;
                Nk::Int=80, η::Real=0.04)
    ks  = range(0, 2π; length=Nk+1)[1:Nk]
    ρ   = zeros(Float64, length(ωs))
    for kx in ks, ky in ks
        Es = real.(eigvals(bloch_hamiltonian(model, kx, ky) + Σ_HF))
        for (iω, ω) in enumerate(ωs)
            for E in Es
                ρ[iω] += η / ((ω - E)^2 + η^2)
            end
        end
    end
    ρ ./ (π * Nk^2)
end
