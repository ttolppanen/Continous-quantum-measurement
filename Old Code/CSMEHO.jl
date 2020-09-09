using DifferentialEquations, Plots, LinearAlgebra
using UsefulFunctionsEtc: ᶜ, σᶻᶜ, σˣᶜ,smeForHD_f, smeForHD_g, photonNumber, ensMean, lowOp

function main()
    Γ = 1.0; ϕ = 0
    a = lowOp(3, isComplex=true)
    ad = a'
    n = ad*a
    𝐻 = (0.5*I + n)
    f(ρ,p,t) = smeForHD_f(𝐻, ρ, a, Γ, ϕ)
    g(ρ,p,t) = smeForHD_g(ρ, a, Γ, ϕ)
    ρ₀ = ᶜ[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]
    dt = 0.1
    Δt = (0.0, 10.0)
    prob = SDEProblem(f,g,ρ₀,Δt)
    t, res = ensMean(prob, 100, Δt, dt)
    res = photonNumber(res, n)
    plot(t, res)
end
@time main()
