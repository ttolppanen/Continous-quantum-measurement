include("./UsefulFunctions.jl")
using .UsefulFunctions
#Please comment out the two above lines ^^^ after running the code once
#To fix this you should add push!(LOAD_PATH, "path to the folder where UsefulFunctions.jl is") this line
#to the startup.jl file and then you can just use it normally using UsefulFunctions...
using DifferentialEquations, LinearAlgebra, Plots

function ex1()
	#s gives the number of allowed states per site. sp.op gives some useful operators such as n, n*𝐼*𝐼 etc
	sp = StantardParameters(Γ=5.0, numOfSys=3, s=3, t=(0.0, 0.01, 10.0), traj=10, isThisMat=false)
	#measOp gives the measurement operator, here it is n (not n*𝐼*𝐼, so only in the hilber space(?) of one site)
	#targets gives the targets of measurement, here we measure the middle site. [1, 5] would measure the first and last site.
    p = ParametersSSEDisorder(W=1.0, U=3.5, J=1.0, sp=sp, measOp=sp.op.n, targets=[3])
	#Ψ₀ is the inital state.
	#It is important to change the state to a real vector, so cvrv(Ψ₀) if it is ket, cmrv(ρ₀) if matrix,
	#here the change happens in the function make1010State.
 	Ψ₀ = make1010State(sp.s,sp.numOfSys)
	#The solution is a list of different trajectories and every trajectory has the solutions for each timestep.
	#[[[t1], [t2], ..., [t]], [jne], ...]
    sol = solveEnsProb(p, Ψ₀)
	#Calculates the mean. If sol = [[1,2,3], [3, 8, -1]] and x->x then the output is [2, 5, 1]
    res = calcMean(sol, x->expVal(x, sp.op.n𝐼[3]))
    plot(sp.t.times, res)
end

function ex2()
	sp = StantardParameters(Γ=5.0, numOfSys=3, s=3, t=(0.0, 0.01, 10.0), traj=10, isThisMat=true)
    p = Parameters(ϕ=1.0, 𝐻=sp.op.nAll, meas=[sp.op.n𝐼[1], sp.op.n𝐼[2]], sp=sp)#ϕ isn't used in the calculations, but it has to be there because of old notation...
    ρ₀ = kronForMany([[0.0+0im 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0],
	[1.0+0im 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
	[0.0+0im 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]])
	ρ₀ = cmrv(ρ₀)
	#solveEnsProb written open...
	prob = SDEProblem(sme_f, sme_g, ρ₀, p.sp.t.Δt, p, saveat=p.sp.t.dt, noise_rate_prototype=zeros(length(ρ₀),2*length(p.meas)))
	enProb = EnsembleProblem(prob, safetycopy=true)
	sol = solve(enProb, SRA1(), EnsembleThreads(), abstol=p.sp.atol, reltol=p.sp.rtol,trajectories=p.sp.traj,dt=p.sp.t.dt)
	sol = ensSolToListMat(sol, sp.dim)
    res = calcMean(sol, x->expVal(x, sp.op.n𝐼[3]))
    plot(sp.t.times, res, ylim=[0,3]) #Due to the Hamiltonian nothing interesting happens here
end

function ex3()#Same as above but with solveEnsProb
	sp = StantardParameters(Γ=5.0, numOfSys=3, s=3, t=(0.0, 0.01, 10.0), traj=10, isThisMat=true)
    p = Parameters(ϕ=1.0, 𝐻=sp.op.nAll, meas=[sp.op.n𝐼[1], sp.op.n𝐼[2]], sp=sp)
    ρ₀ = kronForMany([[0.0+0im 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0],
	[1.0+0im 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0],
	[0.0+0im 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0]])
	ρ₀ = cmrv(ρ₀)
	sol = solveEnsProb(p, ρ₀)
    res = calcMean(sol, x->expVal(x, sp.op.n𝐼[3]))
    plot(sp.t.times, res, ylim=[0,3]) #Due to the Hamiltonian nothing interesting happens here
end

ex3()
