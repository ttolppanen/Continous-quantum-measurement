using DifferentialEquations, Plots,LinearAlgebra, IterTools, UsefulFunctionsEtc, JLD2

function main()
    @load "testiSave.jld2" sol numOfSys s p t traj
    𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
    mean2 = calcMeanForOneSysKet(sol, expVal, 1:length(t.times), p.dim, traj, n; i=2, numOfSys=numOfSys, s=s)
    mean1 = calcMeanForOneSysKet(sol, expVal, 1:length(t.times), p.dim, traj, n; i=1, numOfSys=numOfSys, s=s)
    mean = calcMeanKet(sol, expVal, 1:length(t.times), p.dim, traj, nAll)
    conc = calcMeanKet(sol, iConc, 1:length(t.times), p.dim, traj, s, numOfSys)
    purity = calcMeanKet(sol, x -> real(tr(x*x)), 1:length(t.times), p.dim, traj)
    plot(sol[1].t, mean, ylims=(0,s + 0.1))
    plot!(sol[1].t, mean1)
    plot!(sol[1].t, mean2)
    plot!(sol[1].t, conc)
    plot!(sol[1].t, purity)
end
main()
