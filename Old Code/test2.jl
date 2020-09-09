function main()
	function containsNaN(x)
		for i in x
			if isnan(i)
				println("HEP")
			end
		end
	end
	x = [1, NaN, 3]
	containsNaN(x)
end

main()
