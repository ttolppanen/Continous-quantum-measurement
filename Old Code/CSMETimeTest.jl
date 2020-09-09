using DifferentialEquations, Plots,LinearAlgebra, IterTools
using UsefulFunctionsEtc: ᶜ, σᶻᶜ, σˣᶜ,smeForHD_f, smeForHD_g, lowOp, makeI
using UsefulFunctionsEtc: kronForMany, listOfOperators, calcMeanAndVar, expVal
using UsefulFunctionsEtc: calcMean, cmrv, TimeData, Parameters, makeMPA, makeVPA

function main(traj, numOfQubits)
    s = 2
    𝐼 = makeI(s)
    a = lowOp(s)
    ad = a'
    n = ad*a
    n =  sum(listOfOperators(n, numOfQubits, 𝐼))
    p = Parameters(Γ=0.5, ϕ=0.0, 𝐻=n, op=listOfOperators(a, numOfQubits, 𝐼), dim=s^numOfQubits)
    ρ₀ = cmrv(kronForMany([ᶜ[0.0 0.0; 0.0 1.0] for _ in 1:numOfQubits]))
    t = TimeData(0.0, 0.01, 10.0)
    prob = SDEProblem(smeForHD_f, smeForHD_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),numOfQubits))
    enProb = EnsembleProblem(prob, safetycopy=true)
    sol = solve(enProb, SRA1(), EnsembleThreads(), trajectories=traj)
end
println("____________________________________________________________________")
for q in 1:3
    for traj in [1, 10]
        println("q: $q, traj: $traj ")
        @time main(traj, q)
        println("")
    end
end
