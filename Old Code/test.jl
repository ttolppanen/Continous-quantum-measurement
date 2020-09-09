using UsefulFunctionsEtc, DifferentialEquations, LinearAlgebra, Plots

function main(traj)#t=60 -> n≂0.546
    function fluc(sol, sp)
        res = calcMean(sol, expVal, sp.op.n𝐼[3])
        [(i-res[end])^2 for i in res]
    end
    function entr(sol)
        calcMean(sol, x->vonNeumann(x, 3^3, 3^2))
    end
    sp = StantardParameters(Γ=0.0, numOfSys=5, s=4, t=(0.0, 0.01, 5.0), traj=traj,
    isThisMat=false)
    p = ParametersSSEDisorder(W=5.0, U=3.5, J=1.0, sp=sp, measOp=sp.op.n, targets=[3])
    ρ₀ = make1010State(sp.s,sp.numOfSys)
    sol = solveEnsProb(p, ρ₀)
    #res = entr(sol)
    #fot = calcMean(sol, x->expVal(x, sp.op.n𝐼[3]))
    fot = entr(sol)
    #plot(sp.t.times, res, yaxis=:log10, ylim=[10^-6, 1])
    plot(sp.t.times, fot)
end
@time main(1)
