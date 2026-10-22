# Shared Monte-Carlo utilities used by the demo scripts.
# Deliberately compact — readers should be able to see every line.

using Statistics

"""
    autocorrelation(x; cmax=length(x)÷4) → ρ::Vector{Float64}

Normalised time autocorrelation `ρ(t) = (⟨x_i x_{i+t}⟩ − ⟨x⟩²) / Var(x)` for
`t = 0, 1, …, cmax`.  Plain `O(N·cmax)` estimator — fine for `cmax ≲ 10⁴`.
"""
function autocorrelation(x::AbstractVector{<:Real}; cmax::Int=length(x)÷4)
    N = length(x);   μ = mean(x);  v = var(x)
    ρ = zeros(Float64, cmax + 1)
    for t in 0:cmax
        s = 0.0
        @inbounds for i in 1:(N - t)
            s += (x[i] - μ) * (x[i+t] - μ)
        end
        ρ[t+1] = s / ((N - t) * v)
    end
    ρ
end

"""
    tau_int(x; c=6) → τ::Float64

Madras-Sokal integrated autocorrelation time, self-consistently truncated at
the first `t*` with `t* ≥ c·τ(t*)`.
"""
function tau_int(x::AbstractVector{<:Real}; c::Real=6)
    ρ = autocorrelation(x)
    τ = 0.5
    for t in 1:(length(ρ)-1)
        τ += ρ[t+1]
        if t ≥ c * τ
            return max(τ, 0.5)
        end
    end
    max(τ, 0.5)
end

"""
    block_errors(x; nblocks=20) → (mean, stderr)

Binning analysis: split `x` into `nblocks` equal blocks, return the mean and
the standard error of the block means.
"""
function block_errors(x::AbstractVector{<:Real}; nblocks::Int=20)
    N = length(x);  b = N ÷ nblocks
    μ = [mean(@view x[(k-1)*b+1:k*b]) for k in 1:nblocks]
    mean(μ), std(μ) / sqrt(nblocks)
end

"""
    jackknife(f, data; nblocks=20) → (estimate, stderr)

Delete-`b` jackknife for a derived estimator `f(::Vector) → Real`.
"""
function jackknife(f, data::AbstractVector; nblocks::Int=20)
    N = length(data);  b = N ÷ nblocks
    all_  = f(data[1:nblocks*b])
    drops = [f(vcat(data[1:(k-1)*b], data[k*b+1:nblocks*b])) for k in 1:nblocks]
    μ = mean(drops)
    σ = sqrt((nblocks - 1) / nblocks * sum((d - μ)^2 for d in drops))
    all_, σ
end

"""
    block_curve(x; max_block=length(x)÷10) → (bs, sigmas)

Variance of block means as a function of block size — used to *see* when the
naive `σ/√N` becomes correct.  Plateau height ≈ true standard error.
"""
function block_curve(x::AbstractVector{<:Real}; max_block::Int=length(x)÷10)
    N = length(x)
    bs = Int[];  σs = Float64[]
    b = 1
    while b ≤ max_block
        nb = N ÷ b
        μ  = [mean(@view x[(k-1)*b+1:k*b]) for k in 1:nb]
        push!(bs, b);  push!(σs, std(μ) / sqrt(nb))
        b = max(b + 1, round(Int, b * 1.6))
    end
    bs, σs
end
