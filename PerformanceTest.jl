using DifferentialEquations, Plots, LinearAlgebra, TimerOutputs, UsefulFunctionsEtc

function main()
    q = 4
    s = 2
    u0 = cvrv(kronForMany([complex([0.0, 1.0]) for _ in 1:q]))
    du = zeros(length(u0),q)
    𝐼 = makeI(s)
    a = lowOp(s)
    ad = a'
    n = ad*a
    n =  sum(listOfOperators(n, q, 𝐼))
    op = [(kronForMany(a, 𝐼, 1, q), kronForMany(a, 𝐼, 1, q))]
    p = ParametersSSE(Γ=0.5, ϕ=0.0, 𝐻=n, op=op, dim=s^q)
    @time sseS_g(du,u0,p,0.1)
end

function mainF()
    to = TimerOutput()

    q = 4
    s = 2
    u0 = cvrv(kronForMany([complex([0.0, 1.0]) for _ in 1:q]))
    du = similar(u0)
    𝐼 = makeI(s)
    a = lowOp(s)
    ad = a'
    n = ad*a
    n =  sum(listOfOperators(n, q, 𝐼))
    op = [(kronForMany(a, 𝐼, 1, q), kronForMany(a, 𝐼, 1, q))]
    p = ParametersSSE(Γ=0.5, ϕ=0.0, 𝐻=n, op=op, dim=s^q)
    @time sseS_f(du,u0,p,0.1)
end
mainF()
