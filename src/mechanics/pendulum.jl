using DifferentialEquations

function pendulum!(du, u, p, t)
    g, L, damping = p
    θ, ω = u
    du[1] = ω
    du[2] = -(g / L) * sin(θ) - damping * ω
end

function simulate_pendulum(; θ0=1.2, ω0=0.0,
                             g=9.81, L=1.0, damping=0.0,
                             tspan=(0.0, 10.0))
    prob = ODEProblem(pendulum!, [θ0, ω0], tspan, (g, L, damping))
    solve(prob, Tsit5(); reltol=1e-9, abstol=1e-9)
end

function small_angle(t; θ0=1.2, ω0=0.0, g=9.81, L=1.0)
    Ω = sqrt(g / L)
    θ0 * cos(Ω * t) + (ω0 / Ω) * sin(Ω * t)
end
