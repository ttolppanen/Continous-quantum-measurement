using DifferentialEquations, Plots, LinearAlgebra, TimerOutputs
using UsefulFunctionsEtc: 𝑖, ᶜ, σᶻᶜ, σˣᶜ, expVal, lowOp, com, antiCom, ensMean, 𝒟, partialTrace, expVal
using UsefulFunctionsEtc: matToComp, smeForHD_f, smeForHD_g, photonNumber, lowOp, kronForMany, listOfOperators
using UsefulFunctionsEtc: calcMeanAndVar, calcMean, cmrv, rvcm, TimeData, makeI, Parameters

function main()
    q = 4
    s = 2
    u0 = cmrv(kronForMany([ᶜ[0.0 0.0; 0.0 1.0] for _ in 1:q]))
    du = zeros(length(u0),q)
    𝐼 = makeI(s)
    a = lowOp(s)
    ad = a'
    n = ad*a
    n =  sum(listOfOperators(n, q, 𝐼))
    p = Parameters(Γ=0.5, ϕ=0.0, 𝐻=n, op=listOfOperators(a, q, 𝐼), dim=s^q)
    @time smeForHD_g(du,u0,p,0.1)
end

function mainF()
    to = TimerOutput()

    q = 4
    s = 2
    u0 = cmrv(kronForMany([ᶜ[0.0 0.0; 0.0 1.0] for _ in 1:q]))
    du = similar(u0)
    𝐼 = makeI(s)
    a = lowOp(s)
    ad = a'
    n = ad*a
    n =  sum(listOfOperators(n, q, 𝐼))
    p = Parameters(Γ=0.5, ϕ=0.0, 𝐻=n, op=listOfOperators(a, q, 𝐼), dim=s^q)
    @time smeForHD_f(du,u0,p,0.1)
end
mainF()
