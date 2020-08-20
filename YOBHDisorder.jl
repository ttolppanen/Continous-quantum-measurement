@everywhere using UsefulFunctionsEtc, DifferentialEquations, LinearAlgebra, Plots
using FileIO

function main(traj)
    p = ParametersSSEDisorder(numOfSys=5,s=4,Γ=0.0,target=2,W=1.0,U=3.5,J=1.0,traj=traj,t=TimeData(0.0, 0.01, 5.0))
    ρ₀ = make1010State(p.s,p.numOfSys)
    prob = SDEProblem(sseS_f, sseS_g, ρ₀, p.t.Δt, p, saveat=p.t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
    enProb = EnsembleProblem(prob, prob_func=setNewBH_prob_func, safetycopy=true)
    sol = @time solve(enProb, SRA1(), EnsembleDistributed(), abstol=1e-3, reltol=1e-3,trajectories=traj,dt=p.t.dt)
    save("BHDres1000.jld2", "sol", sol, "p", p, compress=true)
end
main(10)
