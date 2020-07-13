using DifferentialEquations, Plots,LinearAlgebra, IterTools, UsefulFunctionsEtc

function main(traj)
    numOfSys = 3
    s = 2
    𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
    p = Parameters(Γ=0.5, ϕ=0.0, 𝐻=nAll, op=listOfOperators(a, numOfSys, 𝐼), dim=s^numOfSys)
    ρ₀ = cmrv(kronForMany([ᶜ[0.0 0.0; 0.0 1.0] for _ in 1:numOfSys]))
    t = TimeData(0.0, 0.01, 10.0)
    prob = SDEProblem(smeForHD_f, smeForHD_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),numOfSys))
    enProb = EnsembleProblem(prob, safetycopy=true)
    sol = solve(enProb, SRA1(), EnsembleThreads(), abstol=1e-4, reltol=1e-5, trajectories=traj,dt=t.dt)
    mean1 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n;i=1, numOfSys=numOfSys, s=s)
    mean2 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n;i=2, numOfSys=numOfSys, s=s)
    mean3 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n;i=3, numOfSys=numOfSys, s=s)
    plot(sol[1].t, mean1)
    plot!(sol[1].t, mean2)
    plot!(sol[1].t, mean3)
end
@time main(1)
