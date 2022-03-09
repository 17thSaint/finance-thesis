
function get_input_path(h,lambda,bets,a,which)
	file_list = readdir()
	number = length(file_list)
	for i in 1:number
		name = file_list[i]
		separated = split(name,"-")
		if length(separated) < 14 || separated[1] != "path"
			continue
		end
		h_here = parse(Float64,separated[3])
		ft_here = parse(Int,separated[5])
		timecount_here = parse(Int,separated[7])
		lambda_here = parse(Float64,separated[9])
		bets_here = parse(Float64,separated[11])
		a_here = parse(Float64,separated[13])
		which_here = parse(Float64,split(separated[15],".")[1])
		if h == h_here && lambda_here == lambda && bets == bets_here && a == a_here && which_here == which
			return true,timecount_here, ft_here
		end
	end
	return false
end


