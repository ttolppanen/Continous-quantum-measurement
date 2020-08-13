using Plots, UsefulFunctionsEtc, JLD2

function main()
    @time @load "BHDres1000.jld2" sol numOfSys s p t traj
    𝐼, a, ad, n, nAll = make_𝐼_a_ad_n_nAll(s, numOfSys)
    println("hep")
    @time begin
        endMean1 = ithMeanOneSysKet(sol, p, expVal, op.n; i=length(t.times), sysᵢ=2)
        var1 = calcMeanForOneSysKet(sol, x -> (expVal(x, n)-mean1[end])^2, 1:length(t.times), p.dim, traj; i=2, numOfSys=numOfSys, s=s)
        endMean1 = ithMeanOneSysKet(sol, p, expVal, op.n; i=length(t.times), sysᵢ=4)
        var2 = calcMeanForOneSysKet(sol, x -> (expVal(x, n)-mean2[end])^2, 1:length(t.times), p.dim, traj; i=4, numOfSys=numOfSys, s=s)
    end
    println("asd")
    plot(sol[1].t, (var1 + var2)/2, yaxis=:log10)
end
@time main()
