using HDF5, Cubature, SpecialFunctions, HypergeometricFunctions

function frac_brown_wiki2(h,n,t_fin)
	times = [i*t_fin/n for i in 0:n]
	dB = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	bh = fill(0.0,n)
	for j in 1:n
		#if j%(0.05*n) == 0
		#	println("H=",h,", ",100*j/n,"%",", Noise")
		#end
		for i in 0:j-1
			function integrand(s::Float64)
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

function calc_b2(h,avging_counts,sigma0=0.55804)
	steps = 100
	final_time = 1
	motion = [frac_brown_wiki2(h,steps,final_time)[2] for k in 1:avging_counts]
	th_exp_squared = (sigma0^2)*gamma(2-2*h)/(2*h*gamma(1.5-h)*gamma(0.5+h))
	exp_squared = [0.0 for i in 1:avging_counts]
	for i in 1:avging_counts
		exp_squared[i] = (last(motion[i])^2)/avging_counts
	end
	return sum(exp_squared), th_exp_squared
end

function calc_covar(h,avging_counts,interaction_time=500,sigma0=0.55804)
	final_time = 500
	now_time = final_time
	steps = 500
	motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:avging_counts]
	covar = [0.0 for i in 1:steps-1]
	coeff = (sigma0^2)*gamma(2-2*h)/(4*h*gamma(1.5-h)*gamma(0.5+h))
	th_covar = [coeff*(now_time^(2*h)+(i*final_time/steps)^(2*h)-(abs(now_time-i*final_time/steps))^(2*h)) for i in 1:steps-1]
	for j in 1:avging_counts
		for i in 1:steps-1
			covar[i] += motion[j][interaction_time]*motion[j][i]/avging_counts
		end
	end

	return covar,th_covar
end

function write_data_hdf5(version,h,data)
	#println("Starting Data Write: H=$h, $count")
	binary_file_pos = h5open("fBM-$version-h-$h.hdf5","w")
	create_group(binary_file_pos,"values")
	vals = binary_file_pos["values"]
	vals["deets_th"] = data[2]
	vals["deets_exp"] = data[1]
	close(binary_file_pos)
	#println("Data Added, File Closed: H=$h")
end

ver = parse(Int64,ARGS[1])
vers = ["b2","covar"]
avg_counts = parse(Int64,ARGS[2])
h_start = 0.25
h_end = 0.75
dh = 0.01
h_count = Int((h_end-h_start)/dh)
hs = [h_start + dh*i for i in 0:h_count]
for i in 1:h_count+1
	if ver == 1
		data = calc_b2(hs[i],avg_counts)
	else
		data = calc_covar(hs[i],avg_counts)	
	end
	write_data_hdf5(vers[ver],hs[i],data)	
end







