using DifferentialEquations, Plots,LinearAlgebra, IterTools
using UsefulFunctionsEtc: ᶜ, σᶻᶜ, σˣᶜ,smeForHD_f, smeForHD_g, lowOp, makeI
using UsefulFunctionsEtc: kronForMany, listOfOperators, calcMeanAndVar, expVal
using UsefulFunctionsEtc: calcMean, cmrv, rvcm, TimeData, calcMean

function main(traj)
    function makeList(minTol)
        res = [0.01]
        while res[end] > minTol
            push!(res, res[end]/3)
        end
        res
    end
    function checkForPurity(solutions, traj)
        amount = 0
        for sol in solutions
            for ρ in sol.u
                ρ = rvcm(ρ, 8)
                if real(tr(ρ*ρ)) > 1.1
                    amount += 1
                    break
                end
            end
        end
        amount/traj * 100
    end
    numOfQubits = 3
    s = 2
    𝐼 = makeI(s)
    a = lowOp(s)
    ad = a'
    n = ad*a
    n =  sum(listOfOperators(n, numOfQubits, 𝐼))
    p = (Γ=0.5, ϕ=0.0, 𝐻=n, op=listOfOperators(a, numOfQubits, 𝐼), dim=s^numOfQubits)
    ρ₀ = cmrv(kronForMany([ᶜ[0.0 0.0; 0.0 1.0] for _ in 1:numOfQubits]))
    t = TimeData(0.0, 0.01, 10.0)
    #vals = makeList(0.0001)
    res = []
    vals = reverse([0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001])
    for i in vals
        prob = SDEProblem(smeForHD_f, smeForHD_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),numOfQubits))
        enProb = EnsembleProblem(prob)
        sol = solve(enProb, SRA1(), EnsembleThreads(), abstol=1e-2,reltol=1e-2,trajectories=traj, dt=i)
        push!(res, checkForPurity(sol, traj))
    end
    plot(1:7, res, ylims=(0, 100))
end
@time main(100)
