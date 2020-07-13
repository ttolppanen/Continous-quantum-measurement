using DifferentialEquations, Plots, LinearAlgebra
using UsefulFunctionsEtc: 𝑖, ᶜ, σˣᶜ, σᶻᶜ, com, antiCom, expVal, ensMean

function main()
    ω = 1
    g = 1
    Γ = 0.1
    σ = ᶜ[0.0 0.0; 1.0 0.0]
    ee = ᶜ[1.0 0.0; 0.0 0.0] #σ†σ, |e><e|
    H = ω/2.0 * σᶻᶜ + g * σˣᶜ
    dρ(ρ,param,t) = -𝑖*com(H, ρ) - Γ/2.0*antiCom(ee, ρ) + Γ*expVal(ρ, ee)*ρ
    ρ₀ = ᶜ[0.0 0.0; 0.0 1.0]
    Δt = (0.0, 30.0)
    prob = ODEProblem(dρ, ρ₀, Δt, saveat=0.1)

    rate(ρ,param,t) = Γ
    affect!(integrator) = (integrator.u = ᶜ[0.0 0.0; 0.0 1.0])
    jump = ConstantRateJump(rate, affect!)
    jump_prob = JumpProblem(prob, Direct(), jump)
    sol = ensMean(jump_prob, 5000, Δt, 0.1)
    plot(sol.t, real(sol.u[1]), ylim=(0.0, 1.1), label = "Excited")
    #plot!(sol.t, real(sol.u[4]), ylim=(0.0, 1.1), label = "Ground")
end
main()
