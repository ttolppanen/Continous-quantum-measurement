using DifferentialEquations, Plots, LinearAlgebra
using UsefulFunctionsEtc: 𝑖, ᶜ, com, antiCom, expVal, partialTrace, matToComp

function main()
    function handleSolutions(sol)
        A = []
        B = []
        for ρ in sol.u
            push!(A, partialTrace(ρ, 2, 2))
            push!(B, partialTrace(ρ, 2, 2, traceOverB=false))
        end
        A, B
    end
    ω = 1
    Γ = 0.5
    𝐼 = (1.0 + 0.0*im)*Matrix(I, 2, 2)
    σ = ᶜ[0.0 0.0; 1.0 0.0]
    ee = ᶜ[1.0 0.0; 0.0 0.0] #σ†σ, |e><e|
    H = ω * (kron(ee, 𝐼) + kron(𝐼, ee))
    dρ(ρ,param,t) = -𝑖*com(H, ρ) - Γ/2.0*antiCom(kron(ee, 𝐼), ρ) + Γ*expVal(ρ, kron(ee, 𝐼))*ρ - Γ/2.0*antiCom(kron(𝐼, ee), ρ) + Γ*expVal(ρ, kron(𝐼, ee))*ρ
    ρ₀ = kron(ᶜ[0.5 0.0; 0.0 0.5], ᶜ[0.4 0.0; 0.0 0.6])
    Δt = (0.0, 10.0)
    prob = ODEProblem(dρ, ρ₀, Δt, saveat=0.1)

    rate1(ρ,param,t) = Γ
    affect1!(integrator) = (integrator.u = kron(ᶜ[0.0 0.0; 0.0 1.0], partialTrace(integrator.u, 2, 2, traceOverB=false)))
    jump1 = ConstantRateJump(rate1, affect1!)
    rate2(ρ,param,t) = Γ
    affect2!(integrator) = (integrator.u = kron(partialTrace(integrator.u, 2, 2), ᶜ[0.0 0.0; 0.0 1.0]))
    jump2 = ConstantRateJump(rate2, affect2!)
    jump_prob = JumpProblem(prob, Direct(), jump1, jump2)

    sol = solve(jump_prob, Tsit5())
    A, B = handleSolutions(sol)
    A = matToComp(A)
    B = matToComp(B)
    plot(sol.t, real(A[1,1]), ylim=(0.0, 1.0), label = "A Excited")
    plot!(sol.t, real(A[2,2]), ylim=(0.0, 1.0), label = "A Ground")
    plot!(sol.t, real(B[1,1]), ylim=(0.0, 1.0), label = "B Excited")
    plot!(sol.t, real(B[2,2]), ylim=(0.0, 1.0), label = "B Ground")
end
main()
