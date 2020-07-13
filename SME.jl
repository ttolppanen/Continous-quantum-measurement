using DifferentialEquations, Plots, LinearAlgebra
using UsefulFunctionsEtc: 𝑖, ᶜ, com, antiCom, expVal

function main()
    ω = 1
    Γ = 1
    σ = ᶜ[0.0 0.0; 1.0 0.0]
    ee = ᶜ[1.0 0.0; 0.0 0.0] #σ†σ, |e><e|
    H = ω * ee
    dρ(ρ,param,t) = -𝑖*com(H, ρ) - Γ/2.0*antiCom(ee, ρ) + Γ*expVal(ρ, ee)*ρ
    ρ₀ = ᶜ[0.3 0.0; 0.0 0.7]
    Δt = (0.0, 10.0)
    prob = ODEProblem(dρ, ρ₀, Δt, saveat=0.1)

    rate(ρ,param,t) = Γ*expVal(ρ, ee)
    affect!(integrator) = (integrator.u = ᶜ[0.0 0.0; 0.0 1.0])
    jump = ConstantRateJump(rate, affect!)
    jump_prob = JumpProblem(prob, Direct(), jump)

    sol = solve(jump_prob)
    println(typeof(sol.u))
    plot(sol.t, real(sol[1,1,:]), ylim=(0.0, 1.0), label = "Excited")
    plot!(sol.t, real(sol[2,2,:]), ylim=(0.0, 1.0), label = "Ground")
end
main()
