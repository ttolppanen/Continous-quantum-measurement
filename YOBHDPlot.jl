using Plots, UsefulFunctionsEtc, JLD2, IterTools, LinearAlgebra

function main()
    function plotVar(fileName, Γ, first=false)
        @load fileName sol p
        n2 = kronForMany(p.op.n, p.op.𝐼, 2, p.numOfSys)
        n4 = kronForMany(p.op.n, p.op.𝐼, 4, p.numOfSys)
        endMean2 = ithMeanKet(sol, p, length(p.t.times), expVal, n2)
        endMean4 = ithMeanKet(sol, p, length(p.t.times), expVal, n4)
        var, v = calcMeanAndVarKet(sol, p, x -> 0.5*((expVal(x, n2)-endMean2)^2 + (expVal(x, n4)-endMean4)^2))
        if first
            plot(sol[1].t, var, ribbon=v, yaxis=:log10, legend=:bottomright, label="Γ=$Γ")
        else
            plot!(sol[1].t, var, ribbon=v, label="Γ=$Γ")
        end
    end
    function plotSingle(fileName, Γ, i, first=false)
        @load fileName sol p
        res = singleTrajKet(sol[i], x -> norm(x))
        if first
            plot(sol[1].t, res, label="Γ=$Γ")
        else
            plot!(sol[1].t, res, label="Γ=$Γ")
        end
    end
    function plotEn(fileName, Γ, first=false)
        @load fileName sol p
        e, v = calcMeanAndVarKet(sol, p, vonNeumann, 3^3, 3^2)
        if first
            plot(sol[1].t, e, ribbon=v, label="Γ=$Γ", legend=:bottomright)
        else
            plot!(sol[1].t, e, ribbon=v, label="Γ=$Γ")
        end
    end
    function plotNorm(fileName, Γ, first=false)
        @load fileName sol p
        res, v = calcMeanAndVarKet(sol, p, x -> norm(x))
        if first
            plot(sol[1].t, res, ribbon=v, legend=:bottomright, label="Γ=$Γ")
        else
            plot!(sol[1].t, res, ribbon=v, label="Γ=$Γ")
        end
    end
    function plotF(f, data)
        pl = f(data[1][2], data[1][1], true)
        for (first, dataPair) in flagfirst(data)
            if !first
                pl = f(dataPair[2], dataPair[1])
            end
        end
        pl
    end
    function plotS(i, data)
        pl = plotSingle(data[1][2], data[1][1], i, true)
        for (first, dataPair) in flagfirst(data)
            if !first
                pl = plotSingle(dataPair[2], dataPair[1], i)
            end
        end
        pl
    end
    data = [(0.0, "CalcData/BHD200.jld2"),
    (1.0, "CalcData/BHD200_G1.jld2"),
    (5.0, "CalcData/BHD200_G5.jld2"),
    (15.0, "CalcData/BHD200_G15.jld2")]
    plotS(12, data)
    #=
    plotS(6, data)
    =#

end
@time main()
