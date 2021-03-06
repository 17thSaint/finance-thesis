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
	fBM = read_hdf5_data(h,which,true,n+1)[2]
	#fBM = append!([fBM[1]],fBM)
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

function lang_soln(h,told_steps,noise,gam,m,t_fin,v0,which=rand(1:20))
	t_steps = told_steps-1
	times = [i*t_fin/t_steps for i in 0:t_steps]
	term_one = [0.0 for i in 1:t_steps]
	term_two = [0.0 for i in 1:t_steps]
	c_eta = eta(gam,h,1.0)
	scaled_coeffs = get_scale_inv_vals(h,m,gam)
	noise_term = [ noise[i+1] for i in 1:t_steps]#noise#c_eta.*get_noise(h,Int(t_steps*noise_steps),t_fin,which)
	for i in 1:t_steps
		if i%(0.05*t_steps) == 0
			println("H=",h,", ",100*i/t_steps,"%",", Lang")
		end
		noise_times = append!([0.0],[ l*t_fin/t_steps for l in 1:i])
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

function get_first_guess(h,time_told,t_fin,lambda,beta,gam,noise,x0,v0)
	time_steps = time_told-1
	second_term = [0.0 for i in 1:time_steps]
	first_term = [0.0 for i in 1:time_steps]  # convolution part
	times = [0.0 for i in 1:time_steps]
	sliced_noise = noise#[ noise[i+1] for i in 0:time_steps]
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	for i in 1:time_steps
		if i%(time_steps*0.01) == 0
			println("Getting First Guess: ",100*i/time_steps," %")
		end
		time_now = i*t_fin/time_told
		times[i] = time_now
		#second_term[i] = v0*(1-exp(-beta*time_now)) + x0
		second_term[i] = v0*time_now + x0
		noise_times = append!([0.0],[ l*t_fin/time_told for l in 1:i])
		#println(noise_times)
		for j in 1:length(noise_times)-1
			left_time = noise_times[j]
			right_time = noise_times[j+1]
			
			noises_right = sliced_noise[j+1]
			noises_left = sliced_noise[j]
			
			#exp_part_left = (1-exp(-beta*(time_now - left_time)))/(beta*lambda)
			#exp_part_right = (1-exp(-beta*(time_now - right_time)))/(beta*lambda)
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


function get_goft(h,config,delta_t,noise,lambda,beta,gam,a)
	sliced_noise = [ noise[i+1] for i in 1:length(config)-2]
	velocity = [ (config[i+1] - config[i])/delta_t for i in 1:length(config)-1 ]
	accel = [ (velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1 ]
	vel_for_adding = [velocity[i] for i in 2:length(velocity)]
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	g_of_t = first_coeff_scaled.*accel - sliced_noise#+ (beta*lambda).*vel_for_adding
	return g_of_t#,(beta*lambda).*vel_for_adding,lambda.*accel,sliced_noise
end

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

function get_resids(h,config,delta_t,noise,lambda,beta,gam,a)
	g_stuff = get_goft(h,config,delta_t,noise,lambda,beta,gam,a)
	steps = length(config)
	final_time = steps*delta_t
	fract_deriv = [((i*delta_t)^(1-2*h))*get_fractderiv(h,delta_t,steps,config,0.0)[i] for i in 2:steps-1]
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	coeff = scaled_coeffs[2]*(scaled_coeffs[1]^(2*h-2))*gam#a*gamma(2*h)
	return abs.(g_stuff+coeff.*fract_deriv)
end

function move_position(num_times,chosen,step_size)
	shift_matrix = [0.0 for i = 1:num_times]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,h,gam,num_times,chosen,step_size,delta_t,noise,top_val)
	start_resids = get_resids(h,config,delta_t,noise,lambda,beta,gam,a)
	shift_matrix = move_position(num_times,chosen,step_size)
	new_resids = get_resids(h,config+shift_matrix,delta_t,noise,lambda,beta,gam,a)
	exp_diff = exp.(new_resids - start_resids)
	checking = exp_diff[chosen-2] <= top_val
	if checking
		return config+shift_matrix, 1
	else
		return config, 0, new_resids
	end
end

function main_here(tol,steps,step_size,h,time_told,t_fin,lambda,beta,gam,a,noise,x0,v0,top_val,second)
	println("Starting")
	running_config = append!([0.0,second],get_first_guess(h,time_told,t_fin,lambda,beta,gam,noise,x0,v0)[2][3:time_told])
	samp_freq = Int(0.01*steps)
	time_config = fill(0.0,(time_told,Int(steps/samp_freq)))
	time_resids = fill(0.0,(time_told-2,Int(steps/samp_freq)))
	index = 1
	delta_t = t_fin/time_told
	acc_rate = 0
	num_wrong = [i+2 for i in 1:time_told-2]
	for i in 1:steps
		#if i == steps-10
		#	println("Adding More Steps")
		#	steps += 100
		#end
		
		upper = 4
		if length(num_wrong) < 4
			upper = length(num_wrong)
		end
		
		for k in num_wrong[1:upper]#3 + num_correct:upper
			movement = acc_rej_move(running_config,h,gam,time_told,k,step_size,delta_t,noise,top_val)
			running_config = movement[1]
			acc_rate += movement[2]
		end
		num_wrong = []
		
		current_res = get_resids(h,running_config,delta_t,noise,lambda,beta,gam,a)
		
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
				
			else
				append!(num_wrong,[i+2])
			end
		end
		if all(check_tol)
			if i == 1
				tol *= 0.1
				println("Tolerance too low")
				continue
			end
			println("Solution Found in $i Steps")
			return running_config,current_res,time_config,time_resids,i
		end
		
		# interface data
		if i%(steps*0.01) == 0
			println("Running:"," ",100*i/steps,"%, ","Avg Res: ",mean(current_res),", Acceptance: ",acc_rate,", Number Wrong: ",length(num_wrong))
		end
	end
	
	println("No Solution")
	return time_config,time_resids
end


h = 0.75
time_steps = 40
final_time = 20
lambda = 1.0
beta = 0.0
a = 1.0
gam = 0.5
x0 = 0.0
v0 = 0.0
tolerance = 0.01
mc_steps = 50000
step_size = 0.05*final_time/time_steps
metro_val = 1.000001

#white_noise = get_noise(0.5,time_steps,final_time,1).*fluc_dissp_coeffs("white",beta,lambda,a,0.5)
#colored_noise = get_noise(h,time_steps,final_time,2).*fluc_dissp_coeffs("color",beta,lambda,a,h)
#noise = white_noise + colored_noise
#noise = get_noise(h,time_steps,final_time,2).*fluc_dissp_coeffs("color",beta,lambda,a,h)
#soln = lang_soln(h,time_steps,noise,gam,lambda,final_time,v0,2)
#num_soln = main_here(tolerance,mc_steps,step_size,h,time_steps,final_time,lambda,beta,gam,a,noise,x0,v0,metro_val,soln[2][2])

end_step = 37337
samp_freq = Int(0.01*mc_steps)
which_steps = [1,5,10,25,40,72]
for i in 1:length(which_steps)
	chosens = which_steps[i]
	step_number = samp_freq*chosens
	plot(soln[1][2:end-1],num_soln[4][:,chosens],label="$step_number")
end
#plot(soln[1],num_soln[1],"-p",label="Num")
#plot(soln[1],soln[2],"-k",label="Exact")
xlabel("Time")
ylabel("Position")
title("Residuals for MC Sim of Fractional Langevin Equation, H=$h")
legend()





"fin"
