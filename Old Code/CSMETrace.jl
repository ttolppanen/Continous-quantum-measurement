using DifferentialEquations, Plots,LinearAlgebra, IterTools
using UsefulFunctionsEtc: ᶜ, σᶻᶜ, σˣᶜ,smeForHD_f, smeForHD_g, lowOp, makeI
using UsefulFunctionsEtc: kronForMany, listOfOperators, calcMeanAndVar, expVal
using UsefulFunctionsEtc: calcMean, cmrv, TimeData, calcMean, Parameters

function main(traj)
    numOfQubits = 3
    s = 2
    𝐼 = makeI(s)
    a = lowOp(s)
    ad = a'
    n = ad*a
    n =  sum(listOfOperators(n, numOfQubits, 𝐼))
    p = Parameters(Γ=0.1, ϕ=0.0, 𝐻=n, op=listOfOperators(a, numOfQubits, 𝐼), dim=s^numOfQubits)
    ρ₀ = cmrv(kronForMany([ᶜ[0.0 0.0; 0.0 1.0] for _ in 1:numOfQubits]))
    t = TimeData(0.0, 0.01, 10.0)
    prob = SDEProblem(smeForHD_f, smeForHD_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),numOfQubits))
    enProb = EnsembleProblem(prob, safetycopy=true)
    sol = solve(enProb, SRA1(), EnsembleThreads(), abstol=1e-4, reltol=1e-5,trajectories=traj, progress=true)
    trace = calcMean(sol, x -> real(tr(x)), 1:length(t.times), s^numOfQubits, traj)
    purity = calcMean(sol, x -> real(tr(x*x)), 1:length(t.times), s^numOfQubits, traj)
    hermity = calcMean(sol, x -> real(sum(x - x')), 1:length(t.times), s^numOfQubits, traj)
    plot(sol[1].t, trace)
    plot!(sol[1].t, purity)
    plot!(sol[1].t, hermity)
end
@time main(1)
