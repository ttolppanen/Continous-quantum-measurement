module UsefulFunctionsEtc
    using LinearAlgebra: tr, norm
    using Statistics: mean
    using DifferentialEquations, IterTools, LinearAlgebra, Random, Distributed
    using DifferentialEquations.EnsembleAnalysis
	export 𝑖, ᶜ
    export com, antiCom, expVal, ensMean, matToComp, 𝒟, partialTrace, photonNumber
    export ℋ, kronForMany, makeI, solveOneDensity, lowOp
    export listOfOperators, smeForHD_f, smeForHD_g, calcMean, calcMeanAndVar
    export TimeData, Parameters, cmrv, rvcm, make_𝐼_a_ad_n_nAll, calcMeanForOneSys
	export calcMeanAndVarForOneSys, boseHubbard, singleDetection_f, singleDetection_g
	export iConc, excitedState, groundState, calcMeanKet, calcMeanAndVarKet
	export calcMeanForOneSysKet, calcMeanAndVarForOneSysKet, cvrv, rvcv, sse_f, sse_g
	export ParametersSSE, ParametersSSEDisorder, sseS_f, sseS_g, boseHubbardDisorder
	export Operators, ithMeanOneSysKet, make1010State, setNewBH_prob_func

    const 𝑖 = 1.0im
    const ᶜ = Complex{Float64}
    const σˣᶜ = ᶜ[0.0 1.0; 1.0 0.0]
    const σᶻᶜ = ᶜ[1.0 0.0; 0.0 -1.0]
    const ↑ᶜ = ᶜ[1.0; 0.0]
    const ↓ᶜ = ᶜ[0.0; 1.0]

    struct TimeData
    	dt::Float64
    	Δt::Tuple{Float64,Float64}
    	times::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    	function TimeData(startTime::Float64, dt::Float64, endTime::Float64)
    		new(dt, (startTime, endTime), startTime:dt:endTime)
    	end
    end
	struct Operators
		n::Array{Complex{Float64},2}
		nAll::Array{Complex{Float64},2}
		a::Array{Complex{Float64},2}
		ad::Array{Complex{Float64},2}
		𝐼::Array{Complex{Float64},2}
		function Operators(s::Int64, numOfSys::Int64)
			𝐼 = makeI(s)
		    a = lowOp(s)
		    ad = copy(a')
		    n = ad*a
		    nAll =  sum(listOfOperators(n, numOfSys, 𝐼))
			new(n, nAll, a, ad, 𝐼)
		end
	end
	mutable struct Parameters
		Γ::Float64
		ϕ::Float64
		𝐻::Array{Complex{Float64},2}
		op::Array{Array{Complex{Float64},2},1}
		dim::Int64
		mPA::Array{Array{Complex{Float64},2},1}
		vPA::Array{Float64,1}
		function Parameters(;Γ::Float64,ϕ::Float64,𝐻::Array{Complex{Float64},2},
							op::Array{Array{Complex{Float64},2},1},dim::Int64)
			new(Γ,ϕ,𝐻,op,dim,[complex(zeros(dim,dim)) for _ in 1:5],zeros((dim)^2*2))
		end
	end
	mutable struct ParametersSSE
		Γ::Float64
		ϕ::Float64
		𝐻::Array{Complex{Float64},2}
		op::Array{Tuple{Array{Complex{Float64},2},Array{Complex{Float64},2}},1}
		sumccad::Array{Complex{Float64},2}
		dim::Int64
		vPA::Array{Array{Complex{Float64},1},1}
		function ParametersSSE(;Γ::Float64,ϕ::Float64,𝐻::Array{Complex{Float64},2},
				op::Array{Tuple{Array{Complex{Float64},2},Array{Complex{Float64},2}},1},
				dim::Int64)
			new(Γ,ϕ,𝐻,op,(op[1][1] + op[1][1]'),dim,[complex(zeros(dim)) for _ in 1:4])
		end
	end
	mutable struct ParametersSSEDisorder
		Γ::Float64
		𝐻::Array{Complex{Float64},2}
		W::Float64
		U::Float64
		J::Float64
		meas::Array{Tuple{Array{Complex{Float64},2},Array{Complex{Float64},2}},1}
		target::Int64
		op::Operators
		t::TimeData
		sumccad::Array{Complex{Float64},2}
		dim::Int64
		numOfSys::Int64
		s::Int64
		traj::Int64
		vPA::Array{Array{Complex{Float64},1},1}
		function ParametersSSEDisorder(;
				Γ::Float64,W::Float64,U::Float64,J::Float64, target::Int64,
				t::TimeData, numOfSys::Int64, s::Int64, traj::Int64)

			op = Operators(s, numOfSys)
			egOp = copy(op.𝐼)
		    egOp[1] = -1
		    meas = [(kronForMany(sqrt(Γ)*egOp, op.𝐼, target, numOfSys), kronForMany(sqrt(Γ)*1im*egOp, op.𝐼, target, numOfSys))]

			new(Γ,
			boseHubbardDisorder(Wj=W/J, Uj=U/J, n=op.n, a=op.a, 𝐼=op.𝐼, numOfSys=numOfSys),
			W,
			U,
			J,
			meas,
			target,
			op,
			t,
			(meas[1][1] + meas[1][1]'),
			s^numOfSys,
			numOfSys,
			s,
			traj,
			[complex(zeros(s^numOfSys)) for _ in 1:4])
		end
	end
    function com(A::Array{Complex{Float64},2}, B::Array{Complex{Float64},2})
        A*B - B*A
    end
    function com(A::Array{Complex{Float64},2}, B::Array{Complex{Float64},2}, mPA1::Array{Complex{Float64},2}, mPA2::Array{Complex{Float64},2})
        mul!(mPA1, A, B)
        mul!(mPA2, B, A)
        mPA1 .- mPA2
    end
    function antiCom(A::Array{Complex{Float64},2}, B::Array{Complex{Float64},2})
        A*B + B*A
    end
    function antiCom(A::Array{Complex{Float64},2}, B::Array{Complex{Float64},2}, mPA1::Array{Complex{Float64},2}, mPA2::Array{Complex{Float64},2})
        mul!(mPA1, A, B)
        mul!(mPA2, B, A)
        mPA1 .+ mPA2
    end
    function expVal(s::Array{Complex{Float64},1}, op::Array{Complex{Float64},2})#Jos s on ket
        real(s' * op * s)
    end
    function expVal(ρ::Array{Complex{Float64},2}, op::Array{Complex{Float64},2})#Tiheysoperaattorille
        real(tr(op*ρ))
    end
    function expVal(ρ::Array{Complex{Float64},2}, op::Array{Complex{Float64},2}, mPA1::Array{Complex{Float64},2})#Tiheysoperaattorille
        mul!(mPA1, op, ρ)
        real(tr(mPA1))
    end
    function ensMean(prob, trajectories, Δt, dt)
        times = Δt[1]:dt:Δt[2]
        res = []
        var = []
        sol = solve(prob)
        for t in times
            push!(res, sol(t))
            push!(var, sol(t).^2)
        end
        Threads.@threads for _ in 2:trajectories
            sol = solve(prob)
            for j in 1:length(times)
                res[j] += sol(times[j])
                var[j] += sol(times[j]).^2
            end
        end
        res /= trajectories
        var /= trajectories
        matOfOnes = ones(size(res[1]))
        Threads.@threads for i in length(res)
            var[i] -= res[i]^2*matOfOnes
        end
        times, res, var
    end
    function matToComp(listOfMatrices::Array{Array{Complex{Float64},2},1})
        res = []
        l = length(listOfMatrices[1])
        for (isFirst, m) in flagfirst(listOfMatrices)
            if isFirst
                for i in 1:l
                    push!(res, [m[i]])
                end
            else
                for i in 1:l
                    push!(res[i], m[i])
                end
            end
        end
        dim = Int(sqrt(l))
        reshape(res, (dim, dim))
    end
    function 𝒟(c::Array{Complex{Float64},2}, ρ::Array{Complex{Float64},2}) #Lindblad superoperator
        c*ρ*c' - 0.5*(antiCom(c'*c, ρ))
    end
    function 𝒟(c::Array{Complex{Float64},2}, ρ::Array{Complex{Float64},2}, mPA1::Array{Complex{Float64},2}, mPA2::Array{Complex{Float64},2}, mPA3::Array{Complex{Float64},2}) #Lindblad superoperator
        mul!(mPA1, c', c)
        mPA1 .= antiCom(mPA1, ρ, mPA2, mPA3)
        mul!(mPA2, c, ρ)
        mul!(mPA3, mPA2, c')
        mPA3 .- 0.5*mPA1
    end
    function ℋ(c::Array{Complex{Float64},2}, ρ::Array{Complex{Float64},2}) #Measurement superoperator
        c*ρ + ρ*c' - expVal(ρ, c + c')*ρ
    end
	function ℋ(c::Array{Complex{Float64},2}, ρ::Array{Complex{Float64},2}, mPA1::Array{Complex{Float64},2}, mPA2::Array{Complex{Float64},2}) #Measurement superoperator
        mul!(mPA1, c, ρ)
        mul!(mPA2, ρ, c')
        mPA1 .+= mPA2
        mPA1 .- expVal(ρ, c + c', mPA2)*ρ
    end
	function boseHubbard(;ω::Float64, U::Float64, J::Float64, n::Array{Complex{Float64},2}, a::Array{Complex{Float64},2}, 𝐼::Array{Complex{Float64},2}, numOfSys::Int64)
		nᵢ = kronForMany(n, 𝐼, 1, numOfSys)
		𝐼All = kronForMany(𝐼, 𝐼, 1, numOfSys)
		H = ω*nᵢ - U*0.5*nᵢ*(nᵢ-𝐼All)
		for i in 2:numOfSys
			nᵢ = kronForMany(n, 𝐼, i, numOfSys)
			aᵢ₋₁ = kronForMany(a, 𝐼, i - 1, numOfSys)
			aᵢ = kronForMany(a, 𝐼, i, numOfSys)
			H .+= ω*nᵢ - U*0.5*nᵢ*(nᵢ-𝐼All) + J*(aᵢ₋₁*aᵢ' + aᵢ₋₁'*aᵢ)
		end
		H
	end
	function boseHubbardDisorder(;Wj::Float64, Uj::Float64, n::Array{Complex{Float64},2}, a::Array{Complex{Float64},2}, 𝐼::Array{Complex{Float64},2}, numOfSys::Int64)
		nᵢ = kronForMany(n, 𝐼, 1, numOfSys)
		𝐼All = kronForMany(𝐼, 𝐼, 1, numOfSys)
		H = rand(-Wj:0.001:Wj)*nᵢ - Uj*0.5*nᵢ*(nᵢ-𝐼All)
		for i in 2:numOfSys
			nᵢ = kronForMany(n, 𝐼, i, numOfSys)
			aᵢ₋₁ = kronForMany(a, 𝐼, i - 1, numOfSys)
			aᵢ = kronForMany(a, 𝐼, i, numOfSys)
			H .+= rand(-Wj:0.001:Wj)*nᵢ - Uj*0.5*nᵢ*(nᵢ-𝐼All) + (aᵢ₋₁*aᵢ' + aᵢ₋₁'*aᵢ)
		end
		H
	end
	function lowOp(s::Int64) #â
        a = zeros(s, s)
        for i in 1:s-1
            a[i, i + 1] = sqrt(i)
        end
        complex(a)
    end
    function partialTrace(ρ::Array{Complex{Float64},2}, aDim::Int64, bDim::Int64; traceOverB::Bool=true)::Array{Complex{Float64},2}
        if traceOverB
            A = complex(zeros(aDim, aDim))
            for i in 1:aDim
                iᵨ = 1 + (i - 1)*bDim
                for j in 1:aDim
                    jᵨ = 1 + (j - 1)*bDim
                    for d in 0:bDim-1
                        A[i, j] += ρ[iᵨ+d, jᵨ+d]
                    end
                end
            end
            A
        else
            B = complex(zeros(bDim, bDim))
            for i in 1:bDim
                for j in 1:bDim
                    for d in 0:aDim-1
                        B[i, j] += ρ[i + d*bDim, j + d*bDim]
                    end
                end
            end
            B
        end
    end
    function photonNumber(listOfρ::Array{Array{Complex{Float64},2},1}, n::Array{Complex{Float64},2})
        res = []
        for ρ in listOfρ
            push!(res, expVal(ρ, n))
        end
        res
    end
    function kronForMany(m::Union{Array{Complex{Float64},2}, Array{Complex{Float64},1}}, 𝐼, index, numOfSys)
        if index == numOfSys
            s = m
        else
            s = 𝐼
        end
        for i in reverse(1:numOfSys-1)
            if i == index
                s = kron(m, s)
            else
                s = kron(𝐼, s)
            end
        end
        s
    end
    function kronForMany(m::Union{Array{Array{Complex{Float64},2},1}, Array{Array{Complex{Float64},1},1}})
        s = m[end]
        for (isFirst, mᵢ) in flagfirst(reverse(m))
            if isFirst
            else
                s = kron(mᵢ, s)
            end
        end
        s
    end
    function makeI(size::Int64)
        𝐼 = (1.0 + 0.0*im)*Matrix(I, size, size)
    end
    function solveOneDensity(ρ::Array{Complex{Float64},2}, index::Int64, numOfSys::Int64, dim::Int64)
        sysOnLeft = index - 1
        sysOnRight = numOfSys - index
        if sysOnLeft == 0
            partialTrace(ρ, dim, dim^sysOnRight)
        elseif sysOnRight == 0
            partialTrace(ρ, dim^sysOnLeft, dim, traceOverB=false)
        else
            ρB = partialTrace(ρ, dim^sysOnLeft, dim^(sysOnRight + 1), traceOverB=false)
            partialTrace(ρB, dim, dim^sysOnRight)
        end
    end
    function listOfOperators(op::Array{Complex{Float64},2}, numOfSys::Int64, 𝐼::Array{Complex{Float64},2})::Array{Array{Complex{Float64},2},1}
        res = []
        for i in 1:numOfSys
            push!(res, kronForMany(op, 𝐼, i, numOfSys))
        end
        res
    end
    function calcMean(ensSol, f::Function, t, dim::Int64, traj::Int64, f_args...)
        mean = []
        for tᵢ in t
            sols = get_timestep(ensSol, tᵢ)
            sumₘ = 0.0
            for sol in sols
                sumₘ += f(rvcm(sol, dim), f_args...)
            end
            push!(mean, sumₘ/traj)
        end
        mean
    end
    function calcMeanAndVar(ensSol, f::Function, t, dim::Int64, traj::Int64, f_args...)
        mean = []
        var = []
        for tᵢ in t
            sols = get_timestep(ensSol, tᵢ)
            sumₘ = 0.0
            sumᵥ = 0.0
            for sol in sols
                fVal = f(rvcm(sol, dim), f_args...)
                sumₘ += fVal
                sumᵥ += fVal.^2
            end
            push!(mean, sumₘ/traj)
            push!(var, sumᵥ/traj - (sumₘ/traj).^2)
        end
        mean, var
    end
	function calcMeanForOneSys(ensSol, f::Function, t, dim::Int64, traj::Int64, f_args...;i::Int64, numOfSys::Int64, s::Int64)
        mean = []
        for tᵢ in t
            sols = get_timestep(ensSol, tᵢ)
            sumₘ = 0.0
            for sol in sols
                sumₘ += f(solveOneDensity(rvcm(sol, dim), i, numOfSys, s), f_args...)
            end
            push!(mean, sumₘ/traj)
        end
        mean
    end
    function calcMeanAndVarForOneSys(ensSol, f::Function, t, dim::Int64, traj::Int64, f_args...;i::Int64, numOfSys::Int64, s::Int64)
        mean = []
        var = []
        for tᵢ in t
            sols = get_timestep(ensSol, tᵢ)
            sumₘ = 0.0
            sumᵥ = 0.0
            for sol in sols
                fVal = f(solveOneDensity(rvcm(sol, dim), i, numOfSys, s), f_args...)
                sumₘ += fVal
                sumᵥ += fVal.^2
            end
            push!(mean, sumₘ/traj)
            push!(var, sumᵥ/traj - (sumₘ/traj).^2)
        end
        mean, var
    end
	function calcMeanKet(ensSol, f::Function, t, dim::Int64, traj::Int64, f_args...)
        mean = []
        for tᵢ in t
            sols = get_timestep(ensSol, tᵢ)
            sumₘ = 0.0
            for sol in sols
				sol = rvcv(sol)
                sumₘ += f(sol*sol', f_args...)
            end
            push!(mean, sumₘ/traj)
        end
        mean
    end
    function calcMeanAndVarKet(ensSol, f::Function, t, dim::Int64, traj::Int64, f_args...)
        mean = []
        var = []
        for tᵢ in t
            sols = get_timestep(ensSol, tᵢ)
            sumₘ = 0.0
            sumᵥ = 0.0
            for sol in sols
				sol = rvcv(sol)
                fVal = f(sol*sol', f_args...)
                sumₘ += fVal
                sumᵥ += fVal.^2
            end
            push!(mean, sumₘ/traj)
            push!(var, sumᵥ/traj - (sumₘ/traj).^2)
        end
        mean, var
    end
	function calcMeanForOneSysKet(ensSol, p, i, f::Function, f_args...)
        mean = []
        for tᵢ in 1:length(p.t.times)
            sols = get_timestep(ensSol, tᵢ)
            sumₘ = 0.0
            for sol in sols
				sol = rvcv(sol)
                sumₘ += f(solveOneDensity(sol*sol', i, p.numOfSys, p.s), f_args...)
            end
            push!(mean, sumₘ/p.traj)
        end
        mean
    end
    function calcMeanAndVarForOneSysKet(ensSol, f::Function, t, dim::Int64, traj::Int64, f_args...;i::Int64, numOfSys::Int64, s::Int64)
        mean = []
        var = []
        for tᵢ in t
            sols = get_timestep(ensSol, tᵢ)
            sumₘ = 0.0
            sumᵥ = 0.0
            for sol in sols
				sol = rvcv(sol)
                fVal = f(solveOneDensity(sol*sol', i, numOfSys, s), f_args...)
                sumₘ += fVal
                sumᵥ += fVal.^2
            end
            push!(mean, sumₘ/traj)
            push!(var, sumᵥ/traj - (sumₘ/traj).^2)
        end
        mean, var
    end
	function ithMeanOneSysKet(ensSol, p, f, f_args...; i, sysᵢ)
        mean = 0.0
        sols = get_timestep(ensSol, i)
        for sol in sols
			sol = rvcv(sol)
			mean += f(solveOneDensity(sol*sol', sysᵢ, p.numOfSys, p.s), f_args...)
        end
        mean /= p.traj
    end
	function smeForHD_f(dρ::Array{Float64,1},ρ::Array{Float64,1},p,t::Float64)
        p.mPA[4] .= rvcm(ρ, p.dim)
        p.mPA[5] .= -𝑖*com(p.𝐻, p.mPA[4], p.mPA[1], p.mPA[2])
        for σ in p.op
            p.mPA[5] .+= p.Γ*𝒟(σ*exp(𝑖*p.ϕ), p.mPA[4], p.mPA[1], p.mPA[2], p.mPA[3])
        end
        dρ .= cmrv(p.mPA[5])
    end
	function smeForHD_g(dρ::Array{Float64,2},ρ::Array{Float64,1},p,t::Float64)
	    p.mPA[4] .= rvcm(ρ, p.dim)
        for (i,σ) in enumerate(p.op)
            p.mPA[1] .= sqrt(p.Γ) * ℋ(σ*exp(𝑖*p.ϕ),p.mPA[4], p.mPA[2], p.mPA[3])
            p.vPA .= cmrv(p.mPA[1])
            for j in 1:(2*p.dim^2)
                dρ[j,i] = p.vPA[j]
            end
        end
	end
	function singleDetection_f(dρ::Array{Float64,1},ρ::Array{Float64,1},p,t::Float64) #Muista että ħ puuttuu
		p.mPA[4] .= rvcm(ρ, p.dim)
		p.mPA[5] .= -𝑖*com(p.𝐻, p.mPA[4], p.mPA[1], p.mPA[2])
		for σ in p.op
            p.mPA[5] .+= 0.5*p.Γ*𝒟(σ, p.mPA[4], p.mPA[1], p.mPA[2], p.mPA[3])
        end
        dρ .= cmrv(p.mPA[5])
	end
	function singleDetection_g(dρ::Array{Float64,2},ρ::Array{Float64,1},p,t::Float64)
	    p.mPA[4] .= rvcm(ρ, p.dim)
        for (i,σ) in enumerate(p.op)
            p.mPA[1] .= 0.5*sqrt(p.Γ) * ℋ(σ,p.mPA[4], p.mPA[2], p.mPA[3])
            p.vPA .= cmrv(p.mPA[1])
            for j in 1:(2*p.dim^2)
                dρ[j,2*i - 1] = p.vPA[j]
            end
			p.mPA[1] .= 0.5*sqrt(p.Γ) * ℋ(𝑖*σ,p.mPA[4], p.mPA[2], p.mPA[3])
            p.vPA .= cmrv(p.mPA[1])
            for j in 1:(2*p.dim^2)
				dρ[j,2*i] = p.vPA[j]
            end
        end
	end
	function sse_f(dψ::Array{Float64,1},ψ::Array{Float64,1},p,t::Float64)
		ψ .*= 1/norm(rvcv(ψ))
		p.vPA[1] .= rvcv(ψ) #p.vPA[1] = ψ
		mul!(p.vPA[2], -𝑖*p.𝐻, p.vPA[1]) #p.vPA[2] = f
		for c in p.op
			eVal = expVal(p.vPA[1], c + c')
			mul!(p.vPA[3], c, p.vPA[1]) #c*ψ
			mul!(p.vPA[4], c', p.vPA[3]) #c'*c*ψ
			p.vPA[2] .+= -1/4*p.vPA[4] + 1/8*eVal*p.vPA[3] - 1/32*eVal^2*p.vPA[1]
		end
		dψ .= cvrv(p.vPA[2])
	end
	function sse_g(dψ::Array{Float64,2},ψ::Array{Float64,1},p,t::Float64)
		p.vPA[1] .= rvcv(ψ)
		for (i,c) in enumerate(p.op)
			eVal = expVal(p.vPA[1], c + c')
			mul!(p.vPA[3], c, p.vPA[1])
			p.vPA[2] .= 1/2*p.vPA[3] - 1/4*eVal*p.vPA[1]
			g = cvrv(p.vPA[2])
			for j in 1:2*p.dim
				dψ[j,i] = g[j]
			end
		end
	end
	function sseS_f(dψ::Array{Float64,1},ψ::Array{Float64,1},p,t::Float64)
		ψ .*= 1/norm(rvcv(ψ))
		p.vPA[1] .= rvcv(ψ) #p.vPA[1] = ψ
		mul!(p.vPA[2], -𝑖*p.𝐻, p.vPA[1]) #p.vPA[2] = f
		c = p.meas[1][1]
		eVal = expVal(p.vPA[1], p.sumccad)
		mul!(p.vPA[3], c, p.vPA[1]) #c*ψ
		mul!(p.vPA[4], c', p.vPA[3]) #c'*c*ψ
		p.vPA[2] .+= -1/4*p.vPA[4] .+ 1/8*eVal*p.vPA[3] .- 1/32*eVal^2*p.vPA[1]
		dψ .= cvrv(p.vPA[2])
	end
	function sseS_g(dψ::Array{Float64,2},ψ::Array{Float64,1},p,t::Float64)
		p.vPA[1] .= rvcv(ψ)
		for (i,op) in enumerate(p.meas)
			eVal = expVal(p.vPA[1], p.sumccad)#Ei toimi monelle mittaukselle, korjaa
			mul!(p.vPA[3], op[1], p.vPA[1])
			p.vPA[2] .= 1/2*p.vPA[3] .- 1/4*eVal*p.vPA[1]
			g1 = cvrv(p.vPA[2])
			mul!(p.vPA[3], op[2], p.vPA[1])
			p.vPA[2] .= 1/2*p.vPA[3]
			g2 = cvrv(p.vPA[2])
			for j in 1:2*p.dim
				dψ[j,i*2 - 1] = g1[j]
				dψ[j,i*2] = g2[j]
			end
		end
	end
	function cmrv(m::Array{Complex{Float64},2})
        cvrv(vec(m))
    end
    function rvcm(v::Array{Float64,1}, dim::Int64)
        reshape(rvcv(v), dim, dim)
    end
    function cvrv(v::Array{Complex{Float64},1})
        vcat(real(v),imag(v))
    end
    function rvcv(v::Array{Float64,1})
        a = @view v[1:end÷2]
        b = @view v[(end÷2+1):end]
        a + 1im*b
    end
	function iConc(ρ::Array{Complex{Float64},2}, s::Int64, numOfSys::Int64)
		ρₐ = solveOneDensity(ρ, 1, numOfSys, s)
		sqrt(2*max(1-real(tr(ρₐ*ρₐ)), 0))
	end
	function excitedState(s)
	    m = complex(zeros(s,s))
	    m[end] = 1
	    m
	end
	function groundState(s)
	    m = complex(zeros(s,s))
	    m[1] = 1
	    m
	end
	function make1010State(s, numOfSys)
        e = complex(zeros(s))
        e[2] = 1.0
        g = complex(zeros(s))
        g[1] = 1.0
        states = [e]
        for i in 2:numOfSys
            if i % 2 == 0
                push!(states, g)
            else
                push!(states, e)
            end
        end
        cvrv(kronForMany(states))
    end
	function setNewBH_prob_func(prob, i, repeat)
        𝐻 = boseHubbardDisorder(Uj=prob.p.U/prob.p.J, Wj=prob.p.W/prob.p.J; n=prob.p.op.n,
                            a=prob.p.op.a, 𝐼=prob.p.op.𝐼, numOfSys=prob.p.numOfSys)
        newP = deepcopy(prob.p)
        newP.𝐻 = 𝐻
        remake(prob, p=newP)
    end

	#=function make_𝐼_a_ad_n_nAll(s::Int64, numOfSys::Int64)
		𝐼 = makeI(s)
	    a = lowOp(s)
	    ad = copy(a')
	    n = ad*a
	    nAll =  sum(listOfOperators(n, numOfSys, 𝐼))
		𝐼, a, ad, n, nAll
	end
	=#
end
