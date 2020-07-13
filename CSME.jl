using DifferentialEquations, Plots, LinearAlgebra
using UsefulFunctionsEtc: ᶜ, σᶻᶜ, σˣᶜ,smeForHD_f, smeForHD_g, photonNumber, ensMean

function main()
    Γ = 1.0; ϕ = 0.5*π
    σ = ᶜ[0.0 0.0; 1.0 0.0]
    ee = ᶜ[1.0 0.0; 0.0 0.0] #σ†σ, |e><e|
    𝐻 = σᶻᶜ + σˣᶜ
    f(ρ,p,t) = smeForHD_f(𝐻, ρ, σ, Γ, ϕ)
    g(ρ,p,t) = smeForHD_g(ρ, σ, Γ, ϕ)
    ρ₀ = ᶜ[0.0 0.0; 0.0 1.0]
    dt = 0.1
    Δt = (0.0, 10.0)
    prob = SDEProblem(f,g,ρ₀,Δt)
    t, res = ensMean(prob, 1000, Δt, dt)
    res = photonNumber(res, σ'*σ)
    plot(t, res)
end
@time main()
