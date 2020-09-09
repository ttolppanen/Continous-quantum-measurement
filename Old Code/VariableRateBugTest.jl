using DifferentialEquations, Plots, LinearAlgebra, Statistics
using UsefulFunctionsEtc: 𝑖, σˣᶜ, σᶻᶜ, ↑ᶜ, ᶜ

function f(du, u, p, t)
    du .= u
end
function main()
prob = ODEProblem(f, Complex{Float64}[1.0 2.0; 3.0 4.0], (0.0, 1.0))
rate(u,p,t) = norm(u[1] + u[2])
affect!(integrator) = (integrator.u[1] = complex(1.0); integrator.u[2] = complex(2.0);
                       integrator.u[3] = complex(3.0); integrator.u[4] = complex(4.0))
jump = ConstantRateJump(rate, affect!)
jump_prob = JumpProblem(prob, Direct(), jump)
sol = solve(jump_prob)
#plot(sol)
end
main()
