using UsefulFunctionsEtc, DifferentialEquations, LinearAlgebra, Plots
using FileIO

function main()
    function solveMeanAndVar(r, traj, Γ, W)
        sp = StantardParameters(Γ=Γ,numOfSys=5,s=4,t=(0.0, 0.01, 5.0),traj=traj,isThisMat=false)
        p = ParametersSSEDisorder(W=W, U=3.5, J=1.0, target=3, sp=sp)
        ρ₀ = make1010State(sp.s,sp.numOfSys)
        sol = solveEnsProbSSE(p, ρ₀)
        nRes = calcMean(sol, expVal, sp.op.n𝐼[2])
        allRes = [[(i-nRes[end])^2 for i in nRes]]
        for i in 2:r
            p.𝐻 = boseHubbardDisorder(p)
            sol = solveEnsProbSSE(p, ρ₀)
            nRes = calcMean(sol, expVal, sp.op.n𝐼[2])
            push!(allRes, [(i-nRes[end])^2 for i in nRes])
        end
        calcMeanAndVar(allRes, x->x)
    end
    r = 10
    traj = 10
    t = 0.0:0.01:5.0
    finalRes = []
    res, v = solveMeanAndVar(r, traj, 0.0, 1.0)
    push!(finalRes, (res, v))
    Γ = [1.0, 5.0, 15.0]
    for γ in Γ
        res, v = solveMeanAndVar(r, traj, γ, 1.0)
        push!(finalRes, (res, v))
    end
    save("s4r10traj10.jld2", "finalRes", finalRes, compress=true)
end
@time main()
