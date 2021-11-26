include("simulate-frac-brownian.jl")

function fluc_dissp_coeffs(which,beta,lambda,a,h,kBT=1,sigma0=1) #sigma0=0.55804
	if which == "white"
		whi_coeff = 2*beta*lambda*kBT/(sigma0^2)
		return sqrt(whi_coeff)
	elseif which == "color"
		col_coeff = kBT*a*2*gamma(1.5-h)*gamma(0.5+h)/((sigma0^2)*gamma(2*h)*gamma(2-2*h))
		return sqrt(col_coeff)
	else
		return "Error"
	end
end

function get_first_guess(h,time_told,t_fin,beta,lambda,a,white_noise,colored_noise,x0,v0)
	time_steps = time_told-1
	second_term = [0.0 for i in 1:time_steps]
	first_term = [0.0 for i in 1:time_steps]
	times = [0.0 for i in 1:time_steps]
	for i in 1:time_steps
		#println("Getting First Guess: ",100*i/time_steps," %")
		time_now = i*t_fin/time_steps
		times[i] = time_now
		if beta == 0
			second_term[i] = v0*time_now + x0
		else
			second_term[i] = v0*(1-exp(-beta*time_now))/beta + x0
		end
		for j in 1:100
			local_time_change = t_fin/time_steps/100
			noise_time = time_now - (101-j)*local_time_change
			#noises = white_noise[100*(i-1)+j] - colored_noise[100*(i-1)+j]
			noises = white_noise[100*i-j+1] - colored_noise[100*i-j+1]
			if beta == 0
				#exp_part = (time_now-noise_time)/lambda
				exp_part = noise_time/lambda
			else
				#exp_part = (1 - exp(-beta*(time_now-noise_time)))/(beta*lambda)
				exp_part = (1 - exp(-beta*noise_time))/(beta*lambda)
			end
			first_term[i] += noises*exp_part
		end
	end
	full_rez = append!([x0],first_term + second_term)
	times = append!([0.0],times)
	return times,full_rez
end

function get_goft(config,delta_t,white_noise,colored_noise,beta,lambda,x0,v0,a0)
	sliced_white = append!([white_noise[1]],[ white_noise[i*100] for i in 1:length(config)-1])
	sliced_color = append!([colored_noise[1]],[ colored_noise[i*100] for i in 1:length(config)-1])
	velocity = append!([v0],[ (config[i+1] - config[i])/delta_t for i in 1:length(config)-1 ])
	accel = append!([a0],[ (velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1 ])
	g_of_t = lambda.*accel + (beta*lambda).*velocity - sliced_white + sliced_color
	return g_of_t,velocity
end

function get_residuals(h,config,delta_t,white_noise,colored_noise,beta,lambda,a,x0,v0,a0)
	#println("Getting g(t)")
	g_stuff = get_goft(config,delta_t,white_noise,colored_noise,beta,lambda,x0,v0,a0)
	steps = length(config)
	final_time = steps*delta_t
	fract_deriv = [0.0 for i in 1:steps-1]
	for i in 1:steps-1
		top_time_now = i*delta_t
		percent = 100*i/(steps-1)
		#println("Calculating Fractional Derivative: $percent %")
		for j in 1:i
			if j == i
				right = g_stuff[2][j+1]*(top_time_now-delta_t*j/2)^(2*h-2)
			else
				right = g_stuff[2][j+1]*(top_time_now-delta_t*j)^(2*h-2)
			end
			left = g_stuff[2][j]*(top_time_now-delta_t*(j-1))^(2*h-2)
			fract_deriv[i] += 0.5*(left + right)*delta_t/gamma(2*h-1)
		end
	end
	fract_deriv = append!([0.0],fract_deriv)
	return abs.(g_stuff[1]+a.*fract_deriv)
end

function move_position(num_times,chosen,step_size)
	shift_matrix = [0.0 for i = 1:num_times]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,h,num_times,chosen,step_size,delta_t,white_noise,colored_noise)
	start_resids = get_residuals(h,config,delta_t,white_noise,colored_noise,beta,lambda,a,x0,v0,a0)
	shift_matrix = move_position(num_times,chosen,step_size)
	new_resids = get_residuals(h,config+shift_matrix,delta_t,white_noise,colored_noise,beta,lambda,a,x0,v0,a0)
	exp_diff = exp.(new_resids - start_resids)
	checking = [ exp_diff[i] <= 1.001 for i in 1:length(config) ] 
	if all(checking)
		return config+shift_matrix, 1, new_resids
	else
		return config, 0, start_resids
	end
	
	return "Acceptance Calculation Error"
end

function main_here(tol,steps,step_size,h,time_told,t_fin,beta,lambda,a,white_noise,colored_noise,x0,v0,a0)
	println("Starting")
	running_config = get_first_guess(h,time_told,t_fin,beta,lambda,a,white_noise,colored_noise,x0,v0)[2]
	samp_freq = 10
	time_config = fill(0.0,(time_told,Int(steps/samp_freq)))
	time_resids = fill(0.0,(time_told,Int(steps/samp_freq)))
	index = 1
	delta_t = t_fin/time_told
	for i in 1:steps
		for k in 2:time_told
			movement = acc_rej_move(running_config,h,time_told,k,step_size,delta_t,white_noise,colored_noise)
			running_config = movement[1]
		end
		
		current_res = get_residuals(h,running_config,delta_t,white_noise,colored_noise,beta,lambda,a,x0,v0,a0)
		
		if i%samp_freq == 0
			time_config[:,index] = [running_config[x] for x in 1:time_told]
			time_resids[:,index] = [current_res[y] for y in 1:time_told]
			index += 1
		end
		#=
		check_triv = [ running_config[j] == 0.0 for j in 2:time_told ]
		if all(check_triv)
			println("Trivial Solution Found in $i Steps")
			running_config = get_first_guess(h,time_told,t_fin,beta,lambda,a,white_noise,colored_noise,x0,v0)[2]
		end
		=#
		check_tol = [ current_res[j] < tol for j in 2:time_told ]
		if all(check_tol)
			println("Solution Found in $i Steps")
			return running_config,time_config,time_resids,i
		end
	
		if i%(steps*0.01) == 0
			println("Running:"," ",100*i/steps,"%, ","Avg Res: ",mean(current_res))
		end
		#println("Data Added",DateTime(now()))
	end
	
	println("No Solution")
	return time_config,time_resids
end

# for step_size = 0.0001, tough to get avg tol below 0.135

tol = 0.1
mc_steps = 3000
step_size = 0.00001
final_time = 10
time_steps = 10
beta = 0
lambda = 1
a = 1
x0 = 0.0
v0 = 0.0
a0 = 0.0

h = 0.5


#white_noise = noise(0.5,100*time_steps,final_time,1).*fluc_dissp_coeffs("white",beta,lambda,a,h)
#colored_noise = -noise(h,100*time_steps,final_time,2).*fluc_dissp_coeffs("color",beta,lambda,a,h)
#letsgo = main_here(tol,mc_steps,step_size,h,time_steps,final_time,beta,lambda,a,white_noise,colored_noise,x0,v0,a0)
#plot(letsgo[1])

#test_first = get_first_guess(h,time_steps,final_time,beta,lambda,a,white_noise,colored_noise,x0,v0)
#plot(test_first[1],test_first[2])

analytic_soln = lang_soln(h,time_steps,100,a,lambda,final_time,v0,2)


#other_resids = get_residuals(h,analytic_soln[4],final_time/time_steps,white_noise,colored_noise,beta,lambda,a,x0,v0,a0)
#plot(other_resids)

#=
for i in 1:Int(round(letsgo[4]/100))
	#plot(letsgo[3][:,i*10])
	plot(letsgo[2][:,i*10])
end
=#








"fin"
