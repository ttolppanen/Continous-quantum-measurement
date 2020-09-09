using DifferentialEquations, Plots, LinearAlgebra, IterTools
using UsefulFunctionsEtc: 𝑖, ᶜ, com, antiCom, expVal, 𝒟, ensMean, matToComp

function main()
    ω = 1
    Γ = 1
    σ = ᶜ[0.0 0.0; 1.0 0.0]
    ee = ᶜ[1.0 0.0; 0.0 0.0] #σ†σ, |e><e|
    H = ω * ee

    dρ(ρ,param,t) = -𝑖*com(H, ρ) + Γ*𝒟(σ, ρ)
    ρ₀ = ᶜ[0.3 0.0; 0.0 0.7]
    Δt = (0.0, 10.0)
    prob = ODEProblem(dρ, ρ₀, Δt, saveat=0.1)
    solAv = solve(prob)

    dρ1(ρ,param,t) = -𝑖*com(H, ρ) - Γ/2.0*antiCom(ee, ρ) + Γ*expVal(ρ, ee)*ρ
    prob = ODEProblem(dρ1, ρ₀, Δt, saveat=0.1)
    rate(ρ,param,t) = Γ*expVal(ρ, ee)
    affect!(integrator) = (integrator.u = ᶜ[0.0 0.0; 0.0 1.0])
    jump = ConstantRateJump(rate, affect!)
    jump_prob = JumpProblem(prob, Direct(), jump)
    t, res = ensMean(jump_prob, 5000, Δt, 0.1)
    res = matToComp(res)

    plot(solAv.t, real(solAv[1,1,:]), ylim=(0.0, 1.0), label = "Excited")
    plot!(t, real(res[1,1]), ylim=(0.0, 1.0), label = "Excited for mean")
    plot!(t, x -> 0.3*ℯ^(-Γ*x), ylim=(0.0, 1.0), label = "exp(-Γt)", ls=:dash)
end
main()
