using DifferentialEquations, Plots,LinearAlgebra, IterTools
using UsefulFunctionsEtc: ᶜ, σᶻᶜ, σˣᶜ,smeForHD_f, smeForHD_g, lowOp, makeI
using UsefulFunctionsEtc: kronForMany, listOfOperators, calcMeanAndVar, expVal
using UsefulFunctionsEtc: calcMean, cmrv, TimeData, photonNumber, rvcm

function main()
    solvers = Dict(
    "SRA1" => SRA1(),
    "SRA2" => SRA2(),
    #"SRA3" => SRA3(),
    #"SOSRA" => SOSRA(),
    #"SOSRA2" => SOSRA2(),
    "EM" => EM(),
    "LambaEM" => LambaEM(),
    "EulerHeun" => EulerHeun(),
    "LambaEulerHeun" => LambaEulerHeun()
    )
    numOfQubits = 1
    s = 2
    𝐼 = makeI(s)
    a = lowOp(s, isComplex=true)
    ad = a'
    n = ad*a
    n =  sum(listOfOperators(n, numOfQubits, 𝐼))
    p = (Γ=0.5, ϕ=0.0, 𝐻=n, op=listOfOperators(a, numOfQubits, 𝐼), dim=s^numOfQubits)
    ρ₀ = cmrv(kronForMany([ᶜ[0.0 0.1; 0.1 1.0] for _ in 1:numOfQubits]))
    t = TimeData(0.0, 0.01, 1.0)
    prob = SDEProblem(smeForHD_f, smeForHD_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),numOfQubits))
    sol1 = solve(prob,SRA(),dt=t.dt, save_noise=true)
    res = [rvcm(i, s^numOfQubits) for i in sol1.u]
    res = photonNumber(res, n)
    asd = plot(sol1.t, res, label="SRA")
    for i in keys(solvers)
        prob = SDEProblem(smeForHD_f, smeForHD_g, ρ₀, t.Δt, p, noise=NoiseWrapper(sol1.W),saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),numOfQubits))
        sol = solve(prob,solvers[i],dt=t.dt)
        res = [rvcm(i, s^numOfQubits) for i in sol.u]
        res = photonNumber(res, n)
        plot!(sol.t, res, label=i)
    end
    asd
end
@time main()
