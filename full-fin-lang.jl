using HDF5,SpecialFunctions,PyPlot,Statistics,MittagLeffler


function fluc_dissp_coeffs(which,beta,lambda,a,h,kBT=1)
	if which == "white"
		whi_coeff = 2*beta*lambda*kBT
		return sqrt(whi_coeff)
	elseif which == "color"
		col_coeff = kBT*a*2*gamma(1.5-h)*gamma(0.5+h)/(gamma(2*h)*gamma(2-2*h))
		return sqrt(col_coeff)
	else
		return "Error"
	end
end

function read_hdf5_data(h,count,slice=false,len=10000)
	cd("fBM-data")
	file = h5open("fBM-h-$h-$count.hdf5","r")
	data = [read(file["values"],"deets_t"), read(file["values"],"deets_v")]
	cd("..")
	if slice
		data = [read(file["values"],"deets_t")[1:len], read(file["values"],"deets_v")[1:len]]
	end
	return data
end

function get_noise(h,n,t_fin,which=rand(1:20),make_new=false)
	fBM = read_hdf5_data(h,which,true,n)[2]
	fBM = append!([fBM[1]],fBM)
	if make_new
		fBM = frac_brown_wiki2(h,n,t_fin)[2]
	end
	return [t_fin*(fBM[i+1]-fBM[i])/n for i in 1:n]
end

function get_first_guess(h,time_told,t_fin,lambda,beta,noise,x0,v0)
	time_steps = time_told-1
	second_term = [0.0 for i in 1:time_steps]
	first_term = [0.0 for i in 1:time_steps]  # convolution part
	times = [0.0 for i in 1:time_steps]
	sliced_noise = noise#[ noise[i+1] for i in 0:time_steps]
	for i in 1:time_steps
		if i%(time_steps*0.01) == 0
			println("Getting First Guess: ",100*i/time_steps," %")
		end
		time_now = i*t_fin/time_told
		times[i] = time_now
		second_term[i] = v0*(1-exp(-beta*time_now)) + x0
		noise_times = append!([0.0],[ l*t_fin/time_told for l in 1:i])
		#println(noise_times)
		for j in 1:length(noise_times)-1
			left_time = noise_times[j]
			right_time = noise_times[j+1]
			
			noises_right = sliced_noise[j+1]
			noises_left = sliced_noise[j]
			
			exp_part_left = (1-exp(-beta*(time_now - left_time)))/(beta*lambda)
			#println(exp_part_left,", ",i,", ",time_now,", ",left_time)
			exp_part_right = (1-exp(-beta*(time_now - right_time)))/(beta*lambda)
			
			first_term[i] += 0.5*(right_time-left_time)*(noises_left*exp_part_left + noises_right*exp_part_right)
			#println(first_term[i])
		end

	end
	full_rez = append!([x0],first_term + second_term)
	times = append!([0.0],times)
	return times,full_rez
	
end


function get_goft(h,config,delta_t,noise,lambda,beta,a)
	sliced_noise = [ noise[i+1] for i in 1:length(config)-2]
	velocity = [ (config[i+1] - config[i])/delta_t for i in 1:length(config)-1 ]
	accel = [ (velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1 ]
	vel_for_adding = [velocity[i] for i in 2:length(velocity)]
	g_of_t = lambda.*accel + (beta*lambda).*vel_for_adding - sliced_noise
	return g_of_t,velocity,accel
end

h = 0.75
time_steps = 10
final_time = 1
lambda = 1.0
beta = 1.0
a = 1.0
x0 = 0.0
v0 = 0.0
white_noise = get_noise(0.5,time_steps,final_time,1).*fluc_dissp_coeffs("white",beta,lambda,a,0.5)
colored_noise = get_noise(h,time_steps,final_time,2).*fluc_dissp_coeffs("color",beta,lambda,a,h)
noise = white_noise + colored_noise
times = [i*final_time/time_steps for i in 0:time_steps-1]
#guess1 = get_first_guess(h,time_steps,final_time,lambda,beta,noise,x0,v0)
#gt = get_goft(h,guess1[2],final_time/time_steps,noise,lambda,beta,a)
#plot(gt[1])

function get_fractderiv(h,delta_t,steps,og_func,f0,cap=1)
	fract_deriv = [0.0 for i in 1:steps-1]
	times = [i*delta_t for i in 1:steps-1]
	for i in 1:steps-1
		for j in 0:i
			func_here = og_func[i-j+1] - cap*f0
			fact_j = exp(sum([log(i) for i in 1:j]))
			binom_part = gamma(3-2*h)/(fact_j*gamma(3-2*h-j))
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here/(delta_t)^(2-2*h)
		end
	end
	return fract_deriv
end

function get_resids(h,config,delta_t,noise,lambda,beta,a)
	g_stuff = get_goft(h,config,delta_t,noise,lambda,beta,a)
	steps = length(config)
	final_time = steps*delta_t
	fract_deriv = [((i*delta_t)^(1-2*h))*get_fractderiv(h,delta_t,steps,config,0.0)[i] for i in 2:steps-1]
	coeff = a*gamma(2*h)
	return abs.(g_stuff[1]+coeff.*fract_deriv)
end

function move_position(num_times,chosen,step_size)
	shift_matrix = [0.0 for i = 1:num_times]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,h,num_times,chosen,step_size,delta_t,noise,top_val)
	start_resids = get_resids(h,config,delta_t,noise,lambda,beta,a)
	shift_matrix = move_position(num_times,chosen,step_size)
	new_resids = get_resids(h,config+shift_matrix,delta_t,noise,lambda,beta,a)
	exp_diff = exp.(new_resids - start_resids)
	checking = exp_diff[chosen-2] <= top_val
	if checking
		return config+shift_matrix, 1
	else
		return config, 0, new_resids
	end
end

function main_here(tol,steps,step_size,h,time_told,t_fin,lambda,gam,noise,noise_steps,x0,v0,top_val,fixed)
	println("Starting")
	running_config = append!([0.0,fixed],get_first_guess(h,time_told,t_fin,lambda,gam,noise,x0,v0,noise_steps)[2][3:time_told])
	samp_freq = Int(0.01*steps)
	time_config = fill(0.0,(time_told,Int(steps/samp_freq)))
	time_resids = fill(0.0,(time_told-2,Int(steps/samp_freq)))
	index = 1
	delta_t = t_fin/time_told
	acc_rate = 0
	num_correct = 0
	for i in 1:steps
		#if i == steps-10
		#	println("Adding More Steps")
		#	steps += 100
		#end
		
		upper = 6 + num_correct
		if upper > time_told
			upper = time_told
		end
		
		for k in 3 + num_correct:upper
			movement = acc_rej_move(running_config,h,time_told,k,step_size,delta_t,noise,noise_steps,top_val)
			running_config = movement[1]
			acc_rate += movement[2]
		end
		num_correct = 0
		
		current_res = get_resids(h,running_config,delta_t,noise,noise_steps,lambda,gam)
		
		# saving data every 10 steps
		if i%samp_freq == 0
			time_config[:,index] = [running_config[x] for x in 1:time_told]
			time_resids[:,index] = [current_res[y] for y in 1:time_told-2]
			index += 1
		end
		
		# if every time point has residuals less than tolerance then solution is found
		check_tol = [ current_res[j] < tol for j in 1:time_told-2 ]
		for i in 1:time_told-2
			if check_tol[i]
				num_correct += 1
			end
		end
		if all(check_tol)
			if i == 1
				tol *= 0.1
				println("Tolerance too low")
				continue
			end
			println("Solution Found in $i Steps")
			return running_config,time_config,time_resids,i
		end
		
		# interface data
		if i%(steps*0.01) == 0
			println("Running:"," ",100*i/steps,"%, ","Avg Res: ",mean(current_res),", Acceptance: ",acc_rate,", Number Correct: ",num_correct)
		end
	end
	
	println("No Solution")
	return time_config,time_resids
end






"fin"
