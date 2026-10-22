using DifferentialEquations

function radial_path_angle(φ; n=4, a=1, φ₀=0)
    θₙ = 2π/n
    α = θₙ/2
    r₀ = a / (2*sin(α))
    decay = exp(-(φ-φ₀)*tan(α))
    return r₀*decay
end