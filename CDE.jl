using DifferentialEquations, Plots

a,b = complex(1), Complex{Float64}[0.1 0.2; 0.3 0.4]
u₀ = Complex{Float64}[0.1 0.2; 0.3 0.4]
f(u,p,t) = a*u
g(u,p,t) = b*u^2
dt = 0.01
Δt = (0.0, 1.0)
prob = SDEProblem(f,g,u₀, Δt)
sol = solve(prob, EM(), dt=dt)
plot(sol.t, real(sol[1,1,:]))
