using DifferentialEquations

function spring_mass!(du, u, p, t)
    m, k, c, F, ω_drive = p
    x, v = u
    du[1] = v
    du[2] = (-k * x - c * v + F * cos(ω_drive * t)) / m
end

function simulate_spring_mass(; x0=1.0, v0=0.0,
                                m=1.0, k=4.0, c=0.15,
                                F=0.0, ω_drive=0.0,
                                tspan=(0.0, 20.0))
    prob = ODEProblem(spring_mass!, [x0, v0], tspan,
                      (m, k, c, F, ω_drive))
    solve(prob, Tsit5(); reltol=1e-9, abstol=1e-9)
end
