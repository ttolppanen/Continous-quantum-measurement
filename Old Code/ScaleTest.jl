using DifferentialEquations, Plots,LinearAlgebra, UsefulFunctionsEtc, TimerOutputs

function main()
    to = TimerOutput()
    numOfSys = 2
    s = 3
    timing = 0.0
    maxTime = 60.0 * 2
    while true
        while true
            𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
            egOp = copy(𝐼)
            egOp[1] = -1
            𝐻 = boseHubbard(ω=1.0, U=1.0, J=1.0; n=n, a=a, 𝐼=𝐼, numOfSys=numOfSys)
            op = [kronForMany(egOp, 𝐼, 1, numOfSys)]
            p = Parameters(Γ=1.0, ϕ=0.0, 𝐻=𝐻, op=op, dim=s^numOfSys)
            e = excitedState(s)
            g = groundState(s)
            states = [e]
            for _ in 2:numOfSys
                push!(states, g)
            end
            ρ₀ = cmrv(kronForMany(states))
            t = TimeData(0.0, 0.01, 20.0)
            prob = SDEProblem(singleDetection_f, singleDetection_g, ρ₀, t.Δt, p, saveat=t.dt, noise_rate_prototype=zeros(length(ρ₀),2))
            timing = @elapsed (@timeit to "Number of systems: $numOfSys Size: $s" solve(prob, SRA1(), abstol=1e-4, reltol=1e-5,dt=t.dt))
            if(timing ≥ maxTime)
                break
            end
            s += 1
        end
        if(s == 3)
            break
        end
        s = 3
        numOfSys += 1
    end
    show(to, sortby=:name)
end
main()
