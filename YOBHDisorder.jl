@everywhere using UsefulFunctionsEtc, DifferentialEquations, LinearAlgebra, Plots
using JLD2

function main(traj)
    @everywhere function prob_func(prob, i, repeat)
        𝐻 = boseHubbardDisorder(Uj=prob.p.U/prob.p.J, Wj=prob.p.W/prob.p.J; n=prob.p.n,
                            a=prob.p.a, 𝐼=prob.p.𝐼, numOfSys=prob.p.numOfSys)
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
    𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
    egOp = copy(𝐼)
    egOp[1] = -1
    Γ = 0.0
    op = [(kronForMany(sqrt(Γ)*egOp, 𝐼, target, numOfSys), kronForMany(sqrt(Γ)*1im*egOp, 𝐼, target, numOfSys))]
    p = ParametersSSEDisorder(Γ=Γ, W=1.0, U=3.5, J=1.0, op=op, n=n, a=a, 𝐼=𝐼, numOfSys=numOfSys, s=s)
    ρ₀ = make101010State(s,numOfSys)
    t = TimeData(0.0, 0.01, 6.0) #t=60 -> n≂0.546
    prob = SDEProblem(sseS_f, sseS_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
    enProb = EnsembleProblem(prob, prob_func=prob_func, safetycopy=true)
    sol = @time solve(enProb, SRA1(), EnsembleDistributed(), abstol=1e-3, reltol=1e-3,trajectories=traj,dt=t.dt)
    @save "BHDres1000.jld2" sol numOfSys s p t traj
end
main(100)
