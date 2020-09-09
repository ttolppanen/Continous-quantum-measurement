include("./UsefulFunctions.jl")
using .UsefulFunctions
#Please comment out the two above lines ^^^ after running the code once
#To fix this you should add push!(LOAD_PATH, "path to the folder where UsefulFunctions.jl is") this line
#to the startup.jl file and then you can just use it normally using UsefulFunctions...
using DifferentialEquations, LinearAlgebra, Plots, BSON

function main()
    function collectParameters(p, r, startState, measOp)
        params = Dict(
        "Γ" => p.sp.Γ,
        "number of sites" => p.sp.numOfSys,
        "size of one site" => p.sp.s,
        "trajectories" => p.sp.traj,
        "realizations" => r,
        "time" => p.sp.t.times,
        "atol" => p.sp.atol,
        "rtol" => p.sp.rtol,
        "Is this mat" => p.sp.isThisMat,
        "H" => p.𝐻,
        "W" => p.W,
        "U" => p.U,
        "J" => p.J,
        "start state" => startState,
        "targets" => p.targets,
        "measurement operator" => measOp
        )
        params
    end
    function fileCalc(r, traj, Γ, W, U, targets, t)
        sp = StantardParameters(Γ=Γ,numOfSys=5,s=3,t=(t[1], t[2]-t[1],t[end]),traj=traj,isThisMat=false)
        p = ParametersSSEDisorder(W=W, U=U, J=1.0, measOp=sp.op.n, targets=targets, sp=sp)
        ρ₀ = make1010State(sp.s,sp.numOfSys)
        sol = solveEnsProb(p, ρ₀)
        allRes = [sol]
        for i in 2:r
            p.𝐻 = boseHubbardDisorder(p)
            sol = solveEnsProb(p, ρ₀)
            push!(allRes, sol)
        end
        params = collectParameters(p, r, rvcv(ρ₀), sp.op.n)
        allRes, params
    end
	#Calculates data and saves it in a file given realizations and trajectories.
	#The from is [[r1], [r2],...] where [r1] = [[traj1], [traj2], ...]
	#Changing the differential equations, initial state etc should be straightforward
	#with help from the examples.
	#One thing to optimize would be to forget the single trajectories and calculate
	#mean for |Ψ><Ψ| when dealing with kets.
	#Also it might be a good idea to change the file type back to jld2,
	#since it supports compression.

    r = 10
    traj = 10
    t = 0.0:0.01:5.0
    W = 1.0
    U = 3.5
    targets = [3]
	for Γ in [0, 1, 5, 15]
		data, params = fileCalc(r, traj, convert(Float64, Γ), W, U, targets, t)
		#BSON.@save "W1G$Γ.bson" data params
	end
end
@time main()
