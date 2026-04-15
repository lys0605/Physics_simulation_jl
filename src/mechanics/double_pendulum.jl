using DifferentialEquations

function double_pendulum!(du, u, p, t)
    m1, m2, L1, L2, g = p
    θ1, ω1, θ2, ω2 = u
    Δ = θ1 - θ2
    s, c = sin(Δ), cos(Δ)
    den1 = (m1 + m2) * L1 - m2 * L1 * c^2
    den2 = (L2 / L1) * den1

    du[1] = ω1
    du[2] = ( m2 * L1 * ω1^2 * s * c
            + m2 * g * sin(θ2) * c
            + m2 * L2 * ω2^2 * s
            - (m1 + m2) * g * sin(θ1) ) / den1
    du[3] = ω2
    du[4] = ( -m2 * L2 * ω2^2 * s * c
            + (m1 + m2) * g * sin(θ1) * c
            - (m1 + m2) * L1 * ω1^2 * s
            - (m1 + m2) * g * sin(θ2) ) / den2
end

function simulate_double_pendulum(; θ1₀=2.0, θ2₀=2.5,
                                    ω1₀=0.0, ω2₀=0.0,
                                    m1=1.0, m2=1.0,
                                    L1=1.0, L2=1.0, g=9.81,
                                    tspan=(0.0, 15.0))
    u0 = [θ1₀, ω1₀, θ2₀, ω2₀]
    prob = ODEProblem(double_pendulum!, u0, tspan, (m1, m2, L1, L2, g))
    solve(prob, Vern9(); reltol=1e-10, abstol=1e-10)
end
