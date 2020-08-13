using Plots

function plotMat()
	function s(q,s)
		(s^q)^2
	end
	times = [7.87, 13.3, 22.9, 58.5, 119, 257, 25, 97, 240, 116, 480, 308]
	data = [s(2,3), s(2,4), s(2,5), s(2,6), s(2,7), s(2,8), s(3,3), s(3,4), s(3,5), s(4,3), s(4,4), s(5,3)]
	perm = sortperm(data)
	plot(data[perm], times[perm], label="Mat")
end

function plotVec()
	function s(q,s)
		(s^q)
	end
	times = [134, 159, 222, 268, 386, 0.5, 0.7, 1, 3.5, 9.5, 82, 98, 1, 12,
	17,	45, 183, 553, 10, 36, 187, 1019, 12, 443, 168, 4093, 1136]
	data = [s(2,10), s(2,11), s(2,12), s(2,13), s(2,14), s(2,3), s(2,4),
	s(2,5), s(2,6), s(2,7), s(2,8), s(2,9), s(3,3), s(3,4),
	s(3,5), s(3,6), s(3,7), s(3,8), s(4,3), s(4,4), s(4,5), s(4,6),
	s(5,3), s(5,4), s(6,3), s(6,4), s(7,3)]
	perm = sortperm(data)
	plot!(data[perm], times[perm], label="Vec")
end

function plotVecSame()
	function s(q,s)
		(s^q)
	end
	times = [0.5, 0.7, 1, 3.5, 9.5, 82, 1, 12,
	17, 10, 36, 12]
	data = [s(2,3), s(2,4),	s(2,5), s(2,6), s(2,7), s(2,8), s(3,3), s(3,4),
	s(3,5), s(4,3), s(4,4),
	s(5,3)]
	perm = sortperm(data)
	plot!(data[perm], times[perm], label="Vec")
end

function plotMatTest()
	function s(q,s)
		(s^q)
	end
	times = [7.87, 13.3, 22.9, 58.5, 119, 257, 25, 97, 240, 116, 480, 308]
	data = [s(2,3), s(2,4), s(2,5), s(2,6), s(2,7), s(2,8), s(3,3), s(3,4), s(3,5), s(4,3), s(4,4), s(5,3)]
	plot(data[1:6], times[1:6], yaxis=:log10, seriestype = :scatter, color=:blue1, label="Mat q = 2")
	plot!(data[7:9], times[7:9], yaxis=:log10, seriestype = :scatter, color=:orange1, label="Mat q = 3")
	plot!(data[10:11], times[10:11], yaxis=:log10, seriestype = :scatter, color=:blue4, label="Mat q = 4")
	plot!([data[12]], [times[12]], yaxis=:log10, seriestype = :scatter, color=:orange4, label="Mat q = 5")
end

function plotVecTest()
	function s(q,s)
		(s^q)
	end
	times = [134, 159, 222, 268, 386, 0.5, 0.7, 1, 3.5, 9.5, 82, 98, 1, 12,
	17,	45, 183, 553, 10, 36, 187, 1019, 12, 443, 168, 4093, 1136]
	data = [s(2,10), s(2,11), s(2,12), s(2,13), s(2,14), s(2,3), s(2,4),
	s(2,5), s(2,6), s(2,7), s(2,8), s(2,9), s(3,3), s(3,4),
	s(3,5), s(3,6), s(3,7), s(3,8), s(4,3), s(4,4), s(4,5), s(4,6),
	s(5,3), s(5,4), s(6,3), s(6,4), s(7,3)]
	plot!(data[1:12], times[1:12], yaxis=:log10, seriestype=:scatter, color=:blue1, label="Vec")
	plot!(data[13:18], times[13:18],  yaxis=:log10,seriestype=:scatter, color=:orange1, label="Vec")
	plot!(data[19:22], times[19:22],  yaxis=:log10,seriestype=:scatter, color=:red1, label="Vec")
	plot!(data[23:24], times[23:24],  yaxis=:log10,seriestype=:scatter, color=:green, label="Vec")
	#plot!(data[25:26], times[25:26], seriestype=:scatter, label="Vec")
	#plot([data[27]], [times[27]], seriestype=:scatter, label="Vec")
end

function plotVecSameTest()
	function s(q,s)
		(s^q)
	end
	times = [0.5, 0.7, 1, 3.5, 9.5, 82, 1, 12,
	17, 10, 36, 12]
	data = [s(2,3), s(2,4),	s(2,5), s(2,6), s(2,7), s(2,8), s(3,3), s(3,4),
	s(3,5), s(4,3), s(4,4), s(5,3)]
	plot!(data[1:6], times[1:6], seriestype = :scatter, color=:pink1, label="Vec q = 2")
	plot!(data[7:9], times[7:9], seriestype = :scatter, color=:green1, label="Vec q = 3")
	plot!(data[10:11], times[10:11], seriestype = :scatter, color=:pink4, label="Vec q = 4")
	plot!([data[12]], [times[12]], seriestype = :scatter, color=:green4, label="Vec q = 5")
end

#plotVecTest()
plotMatTest()
#plotVecTest()
