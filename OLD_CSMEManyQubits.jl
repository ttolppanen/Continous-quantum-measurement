using DifferentialEquations, Plots, LinearAlgebra, IterTools
using UsefulFunctionsEtc: ᶜ, σᶻᶜ, σˣᶜ,smeForHD_f, smeForHD_g, photonNumber, ensMean, lowOp, makeI, kronForMany, listOfOperators

@time begin
    Γ = 0.5; ϕ = 0.0
    size = 2
    𝐼 = makeI(2)
    a = lowOp(2, isComplex=true)
    ad = a'
    n = ad*a
    n₁, n₂, n₃ = kronForMany(n, 𝐼, 1, 3), kronForMany(n, 𝐼, 2, 3), kronForMany(n, 𝐼, 3, 3)
    n =  n₁ + n₂ + n₃
    𝐻 = n
    op = listOfOperators(a, 3, 𝐼)
    f(ρ,p,t) = smeForHD_f(𝐻, ρ, op, Γ, ϕ, 3, 𝐼)
    g(ρ,p,t) = smeForHD_g(ρ, op, Γ, ϕ, 3, 𝐼)
    ρ₀ = kronForMany([ᶜ[0.0 0.0; 0.0 1.0], ᶜ[0.0 0.0; 0.0 1.0] , ᶜ[0.0 0.0; 0.0 1.0]])
    dt = 0.1
    Δt = (0.0, 10.0)
    prob = SDEProblem(f,g,ρ₀,Δt)
    sim = solve(EnsembleProblem(prob), EnsembleThreads(); trajectories= 100)
    #t, res, var = ensMean(prob, 2, Δt, dt)
    #res = photonNumber(res, n)
    #plot(t, res, ribbon=var[1,1,:])
    #plot!(t, x -> 3*exp(-x*Γ))
end
