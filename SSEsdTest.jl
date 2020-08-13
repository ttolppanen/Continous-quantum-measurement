using DifferentialEquations, Plots,LinearAlgebra, IterTools, UsefulFunctionsEtc

function main(traj)
    numOfSys = 2
    s = 3
    Γ = 1.0
    𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
    egOp = copy(𝐼)
    egOp[1] = -1
    𝐻 = boseHubbard(ω=1.0, U=1.0, J=1.0; n=n, a=a, 𝐼=𝐼, numOfSys=numOfSys)
    op = [kronForMany(egOp, 𝐼, 1, numOfSys)]
    p = Parameters(Γ=Γ, ϕ=0.0, 𝐻=𝐻, op=op, dim=s^numOfSys)
    e = excitedState(s)
    g = groundState(s)
    ρ₀ = cmrv(kronForMany([[0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0], g]))
    t = TimeData(0.0, 0.01, 20.0)
    prob = SDEProblem(singleDetection_f, singleDetection_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
    sol = @time solve(prob, SRA1(), save_noise=true, abstol=1e-4, reltol=1e-5, trajectories=traj,dt=t.dt)
    #mean1 = calcMeanForOneSys(sol, expVal, 1:length(t.times), p.dim, traj, n; i=1, numOfSys=numOfSys, s=s)
    mean1 = [expVal(partialTrace(rvcm(ρ, 9), 3, 3), n) for ρ in sol.u]
    plot(sol.t, mean1)

    ρ₀ = cvrv(kronForMany([[0.0, 0.0im, 1.0], [1.0, 0.0im, 0.0]]))
    op = [kronForMany(sqrt(Γ)*egOp, 𝐼, 1, numOfSys), kronForMany(sqrt(Γ)*1im*egOp, 𝐼, 1, numOfSys)]
    p = ParametersSSE(Γ=Γ, ϕ=0.0, 𝐻=𝐻, op=op, dim=s^numOfSys)
    prob = SDEProblem(sseS_f, sseS_g, ρ₀, t.Δt, p, noise=NoiseWrapper(sol.W), saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
    sol = @time solve(prob, SRA1(), trajectories=traj,dt=t.dt)
    #mean1 = calcMeanForOneSysKet(sol, expVal, 1:length(t.times), p.dim, traj, n; i=1, numOfSys=numOfSys, s=s)
    mean1 = [expVal(partialTrace(rvcv(ket)*rvcv(ket)', 3, 3), n) for ket in sol.u]
    plot!(sol.t, mean1)
end
main(1)
