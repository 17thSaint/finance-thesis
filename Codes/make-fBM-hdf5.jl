#import Pkg; Pkg.add("HDF5")
#import Pkg; Pkg.add("Cubature")
#import Pkg; Pkg.add("SpecialFunctions")
#import Pkg; Pkg.add("HypergeometricFunctions")
using HDF5, Cubature, SpecialFunctions, HypergeometricFunctions

function frac_brown_wiki2(h,n,t_fin)
	times = [i*t_fin/n for i in 0:n]
	dB = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	bh = fill(0.0,n)
	for j in 1:n
		if j%(0.05*n) == 0
			println("H=",h,", ",100*j/n,"%",", fBM")
		end
		for i in 0:j-1
			function integrand(s)
				return (((times[j+1]-s)^(h-0.5))/gamma(h+0.5))*_₂F₁(h-0.5,0.5-h,h+0.5,1-times[j+1]/s)
			end
			part = hquadrature(integrand,times[i+1],times[i+2])[1]*dB[i+1]*n/t_fin
			bh[j] += part
			#println(part,"j=",j,", ","i=",i)
		end
	end
	full_bh = append!([0.0],bh)
	return times, full_bh
end

function write_fBM_data_hdf5(h,count,data)
	println("Starting Data Write: H=$h, $count")
	binary_file_pos = h5open("fBM-h-$h-$count.hdf5","w")
	create_group(binary_file_pos,"values")
	vals = binary_file_pos["values"]
	vals["deets_v"] = data[2]
	vals["deets_t"] = data[1]
	close(binary_file_pos)
	println("Data Added, File Closed: H=$h, $count")
end

t_fin = 10
hs = [0.99,0.98,0.97,0.96,0.94,0.93,0.92,0.91]
ver = parse(Int64,ARGS[1])
h = hs[ver]
num_times = 10
for j in 1:num_times
	data_here = frac_brown_wiki2(h,2000,2000)
	write_fBM_data_hdf5(h,j,data_here)
end

# start 11:10 end 11:17  therefore 20 counts with 13 h's is 30 hours runtime
