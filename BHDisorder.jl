@everywhere using UsefulFunctionsEtc, DifferentialEquations, LinearAlgebra, Plots

function main(traj)
    @everywhere function prob_func(prob, i, repeat)
        𝐻 = boseHubbardDisorder(Uj=prob.p.U/prob.p.J, Wj=prob.p.W/prob.p.J; n=prob.p.op.n,
                            a=prob.p.op.a, 𝐼=prob.p.op.𝐼, numOfSys=prob.p.numOfSys)
        newP = deepcopy(prob.p)
        newP.𝐻 = 𝐻
        remake(prob, p=newP)
    end
    function make101010State(s, numOfSys)
        e = complex(zeros(s))
        e[2] = 1.0
        g = complex(zeros(s))
        g[1] = 1.0
        states = [e]
        for i in 2:numOfSys
            if i % 2 == 0
                push!(states, g)
            else
                push!(states, e)
            end
        end
        cvrv(kronForMany(states))
    end
    numOfSys = 5
    s = 3
    target = 2
    op = Operators(s, numOfSys)
    egOp = copy(op.𝐼)
    egOp[1] = -1
    Γ = 0.0
    meas = [(kronForMany(sqrt(Γ)*egOp, op.𝐼, target, numOfSys), kronForMany(sqrt(Γ)*1im*egOp, op.𝐼, target, numOfSys))]
    p = ParametersSSEDisorder(Γ=Γ, W=1.0, U=3.5, J=1.0, meas=meas, op=op, numOfSys=numOfSys, s=s, traj=traj)
    ρ₀ = make101010State(s,numOfSys)
    t = TimeData(0.0, 0.01, 5.0) #t=60 -> n≂0.546
    prob = SDEProblem(sseS_f, sseS_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
    enProb = EnsembleProblem(prob, prob_func=prob_func, safetycopy=true)
    sol = @time solve(enProb, SRA1(), EnsembleDistributed(), abstol=1e-3, reltol=1e-3,trajectories=traj,dt=t.dt)
    println("Solved")
    endMean1 = ithMeanOneSysKet(sol, p, expVal, op.n; i=length(t.times), sysᵢ=2)
    var1 = calcMeanForOneSysKet(sol, x -> (expVal(x, op.n)-endMean1)^2, 1:length(t.times), p.dim, traj; i=2, numOfSys=numOfSys, s=s)
    endMean2 = ithMeanOneSysKet(sol, p, expVal, op.n; i=length(t.times), sysᵢ=4)
    var2 = calcMeanForOneSysKet(sol, x -> (expVal(x, op.n)-endMean2)^2, 1:length(t.times), p.dim, traj; i=4, numOfSys=numOfSys, s=s)
    plot(sol[1].t, (var1 + var2)/2, yaxis=:log10)
end
@time main(16)
