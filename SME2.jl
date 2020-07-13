using DifferentialEquations, Plots, LinearAlgebra
using UsefulFunctionsEtc: 𝑖, ᶜ, Com, AntiCom, ExpVal

σ = ᶜ[0.0 0.0; 1.0 0.0]
ee = ᶜ[1.0 0.0; 0.0 0.0] #σ†σ, |e><e|
ω = 1
H =  ω * ee
p = (Γ = 1, H = H)
function f(dρ,ρ,p,t)
    dρ .= -𝑖*Com(p.H, ρ) - p.Γ/2.0*AntiCom(ee, ρ) + p.Γ*ExpVal(ρ, ee)*ρ
end
ρ₀ = ᶜ[0.3 0.0; 0.0 0.7]
Δt = (0.0, 10.0)
prob = ODEProblem(f, ρ₀, Δt, saveat=0.1)

rate(ρ,p,t) = Γ*ExpVal(ρ, ee)
affect!(integrator) = (integrator.u = ᶜ[0.0 0.0; 0.0 1.0])
jump = VariableRateJump(rate, affect!)
jump_prob = JumpProblem(prob, Direct(), jump)

sol = solve(jump_prob)
#plot(sol.t, real(sol[1,1,:]), ylim=(0.0, 1.0), label = "Excited")
#plot!(sol.t, real(sol[2,2,:]), ylim=(0.0, 1.0), label = "Ground")
