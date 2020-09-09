include("./UsefulFunctions.jl")
using .UsefulFunctions
#Please comment out the two above lines ^^^ after running the code once
#To fix this you should add push!(LOAD_PATH, "path to the folder where UsefulFunctions.jl is") this line
#to the startup.jl file and then you can just use it normally using UsefulFunctions...
using Plots, BSON, IterTools, LinearAlgebra

function main()
	function vecComponentMat(ic, is, sp)::Array{Complex{Float64},2}
		m = complex(zeros(sp.s,sp.s))
		m[ic,ic] = 1.0
		kronForMany(m, is, sp)
	end
	function calcToList(val, f)
		res = []
		for i in val
			push!(res, calcMean(i, f))
		end
		res
	end
	function ploti(fileName, f1, f2, label)
		BSON.@load fileName data params
		res = calcToList(data, f1)
		m, v = calcMeanAndVar(res, f2)
		plot(params["time"], m, ribbon=v, fillalpha=0.1, label= "$label = $(params[label])")
	end
	function ploti(fileName, f1, f2, label, asd)
		BSON.@load fileName data params
		res = calcToList(data, f1)
		m, v = calcMeanAndVar(res, f2)
		plot!(params["time"], m, ribbon=v, fillalpha=0.1, label= "$label = $(params[label])")
		asd
	end
	function plotFluc(fileName, label)
		BSON.@load fileName data params
		op = Operators(params["size of one site"], params["number of sites"])
		finalRes = []
		for r in data
			n1 = calcMean(r, x->expVal(x, op.n𝐼[2]))
			n2 = calcMean(r, x->expVal(x, op.n𝐼[4]))
			n1 = [(iₙ-n1[end])^2 for iₙ in n1]
			n2 = [(iₙ-n2[end])^2 for iₙ in n2]
			n = 1/2*(n1 .+ n2)
			push!(finalRes, n)
		end
		m, v = calcMeanAndVar(finalRes, x->x)
		plot(params["time"], m, ribbon=v, fillalpha=0.1, label="$label = $(params[label])")
	end
	function plotFluc(fileName, label, asd)
		BSON.@load fileName data params
		op = Operators(params["size of one site"], params["number of sites"])
		finalRes = []
		for r in data
			n1 = calcMean(r, x->expVal(x, op.n𝐼[2]))
			n2 = calcMean(r, x->expVal(x, op.n𝐼[4]))
			n1 = [(iₙ-n1[end])^2 for iₙ in n1]
			n2 = [(iₙ-n2[end])^2 for iₙ in n2]
			n = 1/2*(n1 .+ n2)
			push!(finalRes, n)
		end
		m, v = calcMeanAndVar(finalRes, x->x)
		plot!(params["time"], m, ribbon=v, fillalpha=0.1, label="$label = $(params[label])")
		asd
	end
	function plotVonNeumann(filename, label)
		ploti(filename, x->x*x', x->vonNeumann(partialTrace(x, 3^2, 3^3)), label)
	end
	function plotVonNeumann(filename, label, asd)
		ploti(filename, x->x*x', x->vonNeumann(partialTrace(x, 3^3, 3^2)), label, asd)
	end
	function plotN(filename, label)
		op = Operators(3, 5)
		ploti(filename, x->x*x', x->expVal(x, op.n𝐼[3]), label)
	end
	function plotN(filename, label, asd)
		op = Operators(3, 5)
		ploti(filename, x->x*x', x->expVal(x, op.n𝐼[3]), label, asd)
	end
	W = 0
	#Calculates mean of x->x*x' over trajectories and then mean over these for x->tr(x).
	#Sets the label of the plot to be label= "Γ=$params["Γ"]"
	figure = ploti("testdata.bson", x->x*x', x->tr(x), "Γ")
	ploti("testdata2.bson", x->x*x', x->tr(x), "Γ", figure)#Calculates this and adds it to the same figure
	#for Γ in [1, 5, 15]
	#	plotN("W$(W)U0G$Γ.bson", "Γ", figure)
	#end
	plot!(legend = :topright, ylim=(0, Inf))
	#png("plots/Uudet/W$(W)U0_N") #Save plot, you can also do this after the plot on commandline
	figure
end

figure = main()
