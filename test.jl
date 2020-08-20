using UsefulFunctionsEtc, DifferentialEquations, LinearAlgebra, Plots

function main(traj)
    sp = StantardParameters(Γ=0.02, numOfSys=2, s=3, t=(0.0, 0.01, 20.0), traj=traj, isThisMat=false)
    op = sp.op
    egOp = copy(op.𝐼)
    egOp[1] = -1
    meas = [(kronForMany(sqrt(sp.Γ)*egOp, op.𝐼, 1, sp.numOfSys), kronForMany(sqrt(sp.Γ)*1im*egOp, op.𝐼, 1, sp.numOfSys))]
    𝐻 = boseHubbard(ω=0.2, U=1.0, J=1.0; n=op.n, a=op.a, 𝐼=op.𝐼, numOfSys=sp.numOfSys)
    p = ParametersSSE(𝐻=𝐻, meas=meas, sp=sp)
    ρ₀ = cvrv(kronForMany([[0.0, 0.0im, 1.0], [1.0, 0.0im, 0.0]]))
    sol = solveEnsProbSSE(p, ρ₀)
    n = kronForMany(op.n, 1, sp)
    res = calcMean(sol, sp, expVal, n)
    plot(sp.t.times, res)
end
@time main(100)
