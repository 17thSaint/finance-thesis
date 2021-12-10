using HDF5,SpecialFunctions,PyPlot,Statistics,MittagLeffler

function fluc_dissp_coeffs(which,beta,lambda,a,h,kBT=1,sigma0=1)
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

# gets noise from previous stored data file or makes new if asked
function get_noise(h,n,t_fin,which=rand(1:20),make_new=false)
	fBM = read_hdf5_data(h,which,true,n+1)[2]
	if make_new
		fBM = frac_brown_wiki2(h,n,t_fin)[2]
	end
	return [t_fin*(fBM[i+1]-fBM[i])/n for i in 1:n]
end

function get_scale_inv_vals(h,lambda,gam)
	scaled_time = (lambda/gam)^(0.5/h)
	scaled_pos = fluc_dissp_coeffs("color",0,0,gam,h)*scaled_time/sqrt(gam*lambda)
	return scaled_time,scaled_pos
end

function eta(gam,h,kBT)
	#kBT = 1.0
	return sqrt(2*gam*kBT*gamma(1.5-h)*gamma(0.5+h)/(gamma(2*h)*gamma(2-2*h)))
end

function lang_soln(h,told_steps,noise_steps,noise,gam,m,t_fin,v0,which=rand(1:20))
	t_steps = told_steps-1
	times = [i*t_fin/t_steps for i in 0:t_steps]
	term_one = [0.0 for i in 1:t_steps]
	term_two = [0.0 for i in 1:t_steps]
	c_eta = eta(gam,h,1.0)
	scaled_coeffs = get_scale_inv_vals(h,m,gam)
	noise_term = noise#c_eta.*get_noise(h,Int(t_steps*noise_steps),t_fin,which)
	for i in 1:t_steps
		if i%(0.05*t_steps) == 0
			println("H=",h,", ",100*i/t_steps,"%",", Lang")
		end
		noise_times = append!([0.0],[ l*t_fin/t_steps + j*t_fin/(t_steps*noise_steps) for l in 0:i-1 for j in 1:noise_steps])
		for k in 1:length(noise_times)-1
			left_time = noise_times[k]
			right_time = noise_times[k+1]
			
			mitlef_left = mittleff(2*h,2,(-gam/m)*((left_time*scaled_coeffs[1])^(2*h)))
			mitlef_right = mittleff(2*h,2,(-gam/m)*((right_time*scaled_coeffs[1])^(2*h)))
				
			if v0 != 0.0 && k == length(noise_times)-1
				term_two[i] += v0*right_time*mitlef_right
			end
			
			if k == length(noise_times)-1
				noise_left = noise_term[1]
			else
				noise_left = noise_term[length(noise_times)-1-k]
			end
			noise_right = noise_term[length(noise_times)-k]
				
			term_one[i] += 0.5*(right_time-left_time)*c_eta*(noise_left*left_time*mitlef_left + noise_right*right_time*mitlef_right )
		end
		
	end
	coeff_term_one = (scaled_coeffs[1]^2)/(m*scaled_coeffs[2])
	full_position = append!([0.0],term_one.*coeff_term_one + term_two)
	return times,full_position
end


function get_first_guess(h,time_told,t_fin,lambda,gam,noise,x0,v0,noise_steps)
	time_steps = time_told-1
	second_term = [0.0 for i in 1:time_steps]
	first_term = [0.0 for i in 1:time_steps]  # convolution part
	times = [0.0 for i in 1:time_steps]
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	for i in 1:time_steps
		if i%(time_steps*0.01) == 0
			println("Getting First Guess: ",100*i/time_steps," %")
		end
		time_now = i*t_fin/time_told
		times[i] = time_now
		second_term[i] = v0*time_now + x0
		noise_times = append!([0.0],[ l*t_fin/time_told + k*t_fin/(time_told*noise_steps) for l in 0:i-1 for k in 1:noise_steps])
		#println(noise_times)
		for j in 1:length(noise_times)-1
			left_time = noise_times[j]
			right_time = noise_times[j+1]
			
			noises_right = noise[j+1]
			noises_left = noise[j]
			
			exp_part_left = (time_now - left_time)/first_coeff_scaled
			exp_part_right = (time_now - right_time)/first_coeff_scaled
			
			first_term[i] += 0.5*(right_time-left_time)*(noises_left*exp_part_left + noises_right*exp_part_right)
			#println(first_term[i])
		end

	end
	full_rez = append!([x0],first_term + second_term)
	times = append!([0.0],times)
	return times,full_rez
end

# g(t) is the SDE without the fractional derivative term
function get_goft(h,config,delta_t,noise,lambda,gam)
	sliced_noise = [ noise[i] for i in 2:length(config)-1]
	velocity = [ (config[i+1] - config[i])/delta_t for i in 1:length(config)-1 ]
	accel = [ (velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1 ]
	#accel[2] *= 2
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	g_of_t = first_coeff_scaled.*accel - sliced_noise
	return g_of_t,velocity,accel
end

function get_resids(h,config,delta_t,noise,lambda,gam)
	g_stuff = get_goft(h,config,delta_t,noise,lambda,gam)
	steps = length(config)
	final_time = steps*delta_t
	fract_deriv = [0.0 for i in 1:steps-2]
	for i in 1:steps-2
		top_time_now = (i+1)*delta_t
		if i%((steps-2)*0.01) == 0
			println("Getting Residuals: ",100*i/(steps-2)," %")
		end
		for j in 1:i  # using trapezoid method again for integration 
			if j == i	# doing a cutoff for t=0 which explodes
				right = g_stuff[2][j+1]*(top_time_now-delta_t*j*0.9)^(2*h-2)
			else
				right = g_stuff[2][j+1]*(top_time_now-delta_t*j)^(2*h-2)
			end
			left = g_stuff[2][j]*(top_time_now-delta_t*(j-1))^(2*h-2)
			fract_deriv[i] += 0.5*(left + right)*delta_t/gamma(2*h-1)
		end
	end
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	coeff = gam*scaled_coeffs[2]*(scaled_coeffs[1]^(2*h-2))
	return abs.(g_stuff[1]+coeff.*fract_deriv)
end

function move_position(num_times,chosen,step_size)
	shift_matrix = [0.0 for i = 1:num_times]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,h,num_times,chosen,step_size,delta_t,noise,top_val)
	start_resids = get_resids(h,config,delta_t,noise,lambda,gam)
	shift_matrix = move_position(num_times,chosen,step_size)
	new_resids = get_resids(h,config+shift_matrix,delta_t,noise,lambda,gam)
	exp_diff = exp.(new_resids - start_resids)
	checking = [ exp_diff[i] <= top_val for i in 1:length(config)-2 ] #1.000001
	if all(checking)
		return config+shift_matrix, 1, new_resids
	else
		return config, 0, start_resids
	end
	
	return "Acceptance Calculation Error"
end

function main_here(tol,steps,step_size,h,time_told,t_fin,lambda,gam,noise,x0,v0,top_val)
	println("Starting")
	# getting first config from first guess
	running_config = get_first_guess(h,time_told,t_fin,lambda,gam,noise,x0,v0,noise_steps)[2]
	# only save configuration data for every 10 attempted movements
	samp_freq = 10
	time_config = fill(0.0,(time_told,Int(steps/samp_freq)))
	time_resids = fill(0.0,(time_told-2,Int(steps/samp_freq)))
	index = 1
	delta_t = t_fin/time_told
	for i in 1:steps
		# each MC time step every time point has attempted move, except starting point
		for k in 2:time_told
			movement = acc_rej_move(running_config,h,time_told,k,step_size,delta_t,noise,top_val)
			running_config = movement[1]
		end
		
		current_res = get_resids(h,running_config,delta_t,noise,lambda,gam)
		
		# saving data every 10 steps
		if i%samp_freq == 0
			time_config[:,index] = [running_config[x] for x in 1:time_told]
			time_resids[:,index] = [current_res[y] for y in 1:time_told-2]
			index += 1
		end
		
		# if every time point has residuals less than tolerance then solution is found
		check_tol = [ current_res[j] < tol for j in 1:time_told-2 ]
		if all(check_tol)
			println("Solution Found in $i Steps")
			return running_config,time_config,time_resids,i
		end
		
		# interface data
		if i%(steps*0.01) == 0
			println("Running:"," ",100*i/steps,"%, ","Avg Res: ",mean(current_res))
		end
	end
	
	println("No Solution")
	return time_config,time_resids
end




final_time = 10
time_steps = 10
lambda = 1.0
gam = 1.0
a0 = 0.0
x0 = 0.0
v0 = 0.0
noise_steps = 1
h = 0.5
tol = 0.001
mc_steps = 100
step_size = 1
metro_val = 1.0001

noise = get_noise(h,noise_steps*time_steps,final_time,2).*fluc_dissp_coeffs("color",0,0,gam,h)
#num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,gam,noise,x0,v0,metro_val) ./ 10.0
exact = lang_soln(h,time_steps,noise_steps,noise,gam,lambda,final_time,v0)[2]
resids_exact = get_resids(h,exact,final_time/time_steps,noise,lambda,gam)
plot(resids_exact)





"fin"
