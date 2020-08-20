using UsefulFunctionsEtc, DifferentialEquations, LinearAlgebra, Plots

function main(traj)#t=60 -> n≂0.546
    sp = StantardParameters(Γ=1.0, numOfSys=5, s=3, t=(0.0, 0.01, 5.0), traj=traj)
    p = ParametersSSEDisorder(W=1.0, U=3.5, J=1.0, sp=sp, target=2)
    ρ₀ = make1010State(sp.s,sp.numOfSys)
    sol = solveEnsembleProblem(p, ρ₀)
    nor = calcMeanKet(sol, sp, x -> norm(x))
    plot(sol[1].t, nor)
end
@time main(1)
