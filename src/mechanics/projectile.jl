using DifferentialEquations

function projectile!(du, u, p, t)
    g, drag = p
    x, y, vx, vy = u
    speed = sqrt(vx^2 + vy^2)
    du[1] = vx
    du[2] = vy
    du[3] = -drag * speed * vx
    du[4] = -g - drag * speed * vy
end

function simulate_projectile(; v0=30.0, angle_deg=55.0,
                               g=9.81, drag=0.0,
                               tspan=(0.0, 6.0))
    θ = deg2rad(angle_deg)
    u0 = [0.0, 0.0, v0 * cos(θ), v0 * sin(θ)]
    prob = ODEProblem(projectile!, u0, tspan, (g, drag))
    solve(prob, Tsit5(); reltol=1e-8, abstol=1e-8)
end
