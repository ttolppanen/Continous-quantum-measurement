using DifferentialEquations, Plots,LinearAlgebra, IterTools, UsefulFunctionsEtc

function main(traj)
    numOfSys = 3
    s = 3
    𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
    𝐻 = boseHubbard(ω=1.0, U=10.0, J=1.0; n=n, a=a, 𝐼=𝐼, numOfSys=numOfSys)
    p = Parameters(Γ=1.0, ϕ=0.0, 𝐻=𝐻, op=listOfOperators(a, numOfSys, 𝐼), dim=s^numOfSys)
    e = ᶜ[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]
    g = ᶜ[1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    ρ₀ = cmrv(kronForMany([e, g, e]))
    t = TimeData(0.0, 0.01, 10.0)
    prob = SDEProblem(smeForHD_f, smeForHD_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),numOfSys))
    enProb = EnsembleProblem(prob, safetycopy=true)
    sol = solve(enProb, SRA1(), EnsembleThreads(), abstol=1e-3, reltol=1e-3, trajectories=traj,dt=t.dt)
    mean3 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n; i=3, numOfSys=numOfSys, s=s)
    mean2 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n; i=2, numOfSys=numOfSys, s=s)
    mean1 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n; i=1, numOfSys=numOfSys, s=s)
    mean = calcMean(sol, expVal, 1:length(t.times), p.dim, traj, nAll)
    plot(sol[1].t, mean, ylims=(0,4))
    plot!(sol[1].t, mean1)
    plot!(sol[1].t, mean2)
    plot!(sol[1].t, mean3)
    plot!(sol[1].t, x -> 4*exp(-x*p.Γ))
end
@time main(100)
