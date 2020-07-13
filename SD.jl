using DifferentialEquations, Plots,LinearAlgebra, IterTools, UsefulFunctionsEtc

function main(traj)
    numOfSys = 2
    s = 3
    𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
    egOp = copy(𝐼)
    egOp[1] = -egOp[1]
    𝐻 = boseHubbard(ω=1.0, U=2.0, J=1.0; n=n, a=a, 𝐼=𝐼, numOfSys=numOfSys)
    op = [kronForMany(egOp, 𝐼, 1, numOfSys)]
    p = Parameters(Γ=1.0, ϕ=0.0, 𝐻=𝐻, op=op, dim=s^numOfSys)
    e = excitedState(s)
    g = groundState(s)
    ρ₀ = cmrv(kronForMany([[0.0 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.5], g]))
    t = TimeData(0.0, 0.01, 30.0)
    prob = SDEProblem(singleDetection_f, singleDetection_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
    enProb = EnsembleProblem(prob, safetycopy=true)
    sol = solve(enProb, SRA1(), EnsembleThreads(), abstol=1e-3, reltol=1e-3, trajectories=traj,dt=t.dt)
    mean2 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n; i=2, numOfSys=numOfSys, s=s)
    mean1 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n; i=1, numOfSys=numOfSys, s=s)
    mean = calcMean(sol, expVal, 1:length(t.times), p.dim, traj, nAll)
    conc = calcMean(sol, iConc, 1:length(t.times), p.dim, traj, s, numOfSys)
    purity = calcMean(sol, x -> real(tr(x*x)), 1:length(t.times), p.dim, traj)
    plot(sol[1].t, mean, ylims=(0,s + 0.1))
    plot!(sol[1].t, mean1)
    plot!(sol[1].t, mean2)
    plot!(sol[1].t, conc)
    plot!(sol[1].t, purity)
end
@time main(10)
