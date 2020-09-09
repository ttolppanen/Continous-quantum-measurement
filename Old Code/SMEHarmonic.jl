using DifferentialEquations, Plots, LinearAlgebra
using UsefulFunctionsEtc: 𝑖, ᶜ, expVal, lowOp, com, antiCom, ensMean, matToComp

function main()
    function emit(ρ, a, n)
        a*ρ*ad/expVal(ρ, n)
    end

    size = 3
    ω = 1
    Γ = 1
    a = lowOp(3, isComplex=true)
    ad = a'
    n = ad*a #n̂
    H = ω * (n + 0.5*I)
    dρ(ρ,param,t) = -𝑖*com(H, ρ) - Γ/2.0*antiCom(n, ρ) + Γ*expVal(ρ, n)*ρ
    ρ₀ = ᶜ[0.3 0.0 0.0; 0.0 0.3 0.0; 0.0 0.0 0.3]
    ρ₀ = ρ₀ / tr(ρ₀)
    Δt = (0.0, 10.0)
    prob = ODEProblem(dρ, ρ₀, Δt, saveat=0.1)
    rate(ρ,param,t) = Γ*expVal(ρ, n)
    affect!(integrator) = (integrator.u = emit(integrator.u, a, n))
    jump = ConstantRateJump(rate, affect!)
    jump_prob = JumpProblem(prob, Direct(), jump)

    sol = matToComp(ensMean(jump_prob, 500, Δt, 0.1))
    plot(sol.t, real(sol.u[1,1]), ylim=(0.0, 1.0), label = "0")
    plot!(sol.t, real(sol.u[2,2]), ylim=(0.0, 1.0), label = "1")
    plot!(sol.t, real(sol.u[3,3]), ylim=(0.0, 1.0), label = "2")
end
main()
