@everywhere using UsefulFunctionsEtc, DifferentialEquations, LinearAlgebra, Plots

function main(traj)#t=60 -> n≂0.546
    p = ParametersSSEDisorder(numOfSys=5,s=3,Γ=0.0,target=2,W=1.0,U=3.5,J=1.0,traj=traj,t=TimeData(0.0, 0.01, 5.0))
    ρ₀ = make1010State(p.s,p.numOfSys)
    prob = SDEProblem(sseS_f, sseS_g, ρ₀, p.t.Δt, p, saveat=p.t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
    enProb = EnsembleProblem(prob, prob_func=setNewBH_prob_func, safetycopy=true)
    sol = @time solve(enProb, SRA1(), EnsembleDistributed(), abstol=1e-3, reltol=1e-3,trajectories=traj,dt=p.t.dt)

    endMean1 = ithMeanOneSysKet(sol, p, expVal, p.op.n; i=length(p.t.times), sysᵢ=2)
    var1 = calcMeanForOneSysKet(sol, p, 2, x -> (expVal(x, p.op.n)-endMean1)^2)
    endMean2 = ithMeanOneSysKet(sol, p, expVal, p.op.n; i=length(p.t.times), sysᵢ=4)
    var2 = calcMeanForOneSysKet(sol, p, 4, x -> (expVal(x, p.op.n)-endMean2)^2)
    plot(sol[1].t, (var1 + var2)/2, yaxis=:log10)
end
@time main(16)
