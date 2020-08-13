using DifferentialEquations, Plots,LinearAlgebra, IterTools, UsefulFunctionsEtc, JLD2

function main(traj)
    numOfSys = 2
    s = 3
    𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
    egOp = copy(𝐼)
    egOp[1] = -1
    𝐻 = boseHubbard(ω=1.0, U=2.0, J=1.0; n=n, a=a, 𝐼=𝐼, numOfSys=numOfSys)
    op = [(kronForMany(egOp, 𝐼, 1, numOfSys), kronForMany(1im*egOp, 𝐼, 1, numOfSys))]
    p = ParametersSSE(Γ=1.0, ϕ=0.0, 𝐻=𝐻, op=op, dim=s^numOfSys)
    ρ₀ = cvrv(kronForMany([[0.0, 0.0im, 1.0], [1.0, 0.0im, 0.0]]))
    t = TimeData(0.0, 0.01, 30.0)
    prob = SDEProblem(sseS_f, sseS_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
    enProb = EnsembleProblem(prob, safetycopy=true)
    sol = @time solve(enProb, SRA1(), EnsembleThreads(), abstol=1e-3, reltol=1e-3, trajectories=traj,dt=t.dt)
    @save "testiSave.jld2" sol numOfSys s p t traj
end
main(5)
