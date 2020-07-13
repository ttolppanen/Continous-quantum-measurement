using DifferentialEquations, Plots, LinearAlgebra
using UsefulFunctionsEtc: 𝑖, σˣᶜ, σᶻᶜ, ↑ᶜ, ᶜ, Com, AntiCom, ExpVal

Γ = 1
dρ(ρ,param,t) = ρ
ρ₀ = ᶜ[1.0 0.0; 0.0 0.0]
Δt = (0.0, 10.0)
prob = ODEProblem(dρ, ρ₀, Δt, saveat=0.1)

rate(ρ,param,t) = Γ
affect!(integrator) = (integrator.u = ᶜ[0.1 0.2; 0.3 1.0])
jump = ConstantRateJump(rate, affect!)
jump_prob = JumpProblem(prob, Direct(), jump)

ensemble_prob = EnsembleProblem(jump_prob)
sol = solve(ensemble_prob,Tsit5(),EnsembleThreads(),trajectories=2)
sum = EnsembleSummary(sol, Δt[1]:0.1:Δt[2])
plot(sum)
#plot!(sol.t, real(sol[2,2,:]), ylim=(0.0, 1.1), label = "Ground")
