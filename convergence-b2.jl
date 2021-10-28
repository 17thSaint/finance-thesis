import Pkg; Pkg.add("HDF5")
import Pkg; Pkg.add("Cubature")
import Pkg; Pkg.add("SpecialFunctions")
import Pkg; Pkg.add("HypergeometricFunctions")
using HDF5,SpecialFunctions,Cubature, HypergeometricFunctions

function frac_brown_wiki2(h,n,t_fin)
	times = [i*t_fin/n for i in 0:n]
	dB = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	bh = fill(0.0,n)
	for j in 1:n
		for i in 0:j-1
			function integrand(s)
				return (((times[j+1]-s)^(h-0.5))/gamma(h+0.5))*_₂F₁(h-0.5,0.5-h,h+0.5,1-times[j+1]/s)
			end
			part = hquadrature(integrand,times[i+1],times[i+2])[1]*dB[i+1]*n/t_fin
			bh[j] += part
		end
	end
	full_bh = append!([0.0],bh)
	return times, full_bh
end

function calc_b2(h,avging_counts,motion,sigma0)
	exp_squared = [0.0 for i in 1:avging_counts]
	for i in 1:avging_counts
		exp_squared[i] = (last(motion[i])^2)/avging_counts
	end
	return sum(exp_squared)
end

function write_b2_data_hdf5(h,count,data)
	println("Starting Data Write: $count")
	binary_file_pos = h5open("b2-conv-h-$h-$count.hdf5","w")
	create_group(binary_file_pos,"values")
	vals = binary_file_pos["values"]
	vals["deets"] = data
	close(binary_file_pos)
	println("Data Added, File Closed: $count")
end

times = parse(Int64,ARGS[1])
s0 = 0.5580484388471367
h = 0.75
count_end = 100
count_steps = 10
iters = Int((count_end-100)/count_steps)+1
avg_counts = [100+(i-1)*count_steps for i in 1:iters]
b2_th = (s0^2)*gamma(2-2*h)/(2*h*gamma(1.5-h)*gamma(0.5+h))
for j in 1:times
	b2_vals = [0.0 for i in 1:iters]
	println("Getting Motion")
	motion_here = [frac_brown_wiki2(h,100,1)[2] for k in 1:count_end]
	println("Finished Motion")
	for i in 1:iters
		println(100*(i-1)/iters,"%")
		avging_counts = 100+(i-1)*count_steps
		results = calc_b2(h,avging_counts,motion_here,s0)
		b2_vals[i] = results
	end
	write_b2_data_hdf5(h,j,b2_vals)
end




