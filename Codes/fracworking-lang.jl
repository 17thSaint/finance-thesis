using HDF5,SpecialFunctions,PyPlot,Statistics,MittagLeffler,LaTeXStrings

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
	cd("..")
	cd("fBM-data")
	file = h5open("fBM-h-$h-$count.hdf5","r")
	data = [read(file["values"],"deets_t"), read(file["values"],"deets_v")]
	cd("..")
	cd("Codes")
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
	scaled_time = 1#(lambda/gam)^(0.5/h)
	scaled_pos = 1#fluc_dissp_coeffs("color",0,0,gam,h)*scaled_time/sqrt(gam*lambda)
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
	noise_term = [ noise[i*noise_steps+1] for i in 1:t_steps]#noise#c_eta.*get_noise(h,Int(t_steps*noise_steps),t_fin,which)
	for i in 1:t_steps
		if i%(0.05*t_steps) == 0
			println("H=",h,", ",100*i/t_steps,"%",", Lang")
		end
		noise_times = append!([0.0],[ l*t_fin/t_steps + j*t_fin/(t_steps*1) for l in 0:i-1 for j in 1:1])
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
	sliced_noise = [ noise[i*noise_steps+1] for i in 0:time_steps]
	noise_steps_local = 1
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	for i in 1:time_steps
		if i%(time_steps*0.01) == 0
			println("Getting First Guess: ",100*i/time_steps," %")
		end
		time_now = i*t_fin/time_told
		times[i] = time_now
		second_term[i] = v0*time_now + x0
		noise_times = append!([0.0],[ l*t_fin/time_told + k*t_fin/(time_told*noise_steps_local) for l in 0:i-1 for k in 1:noise_steps_local])
		#println(noise_times)
		for j in 1:length(noise_times)-1
			left_time = noise_times[j]
			right_time = noise_times[j+1]
			
			noises_right = sliced_noise[j+1]
			noises_left = sliced_noise[j]
			
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
function get_goft(h,config,delta_t,noise,noise_steps,lambda,gam)
	sliced_noise = [ noise[i*noise_steps+1] for i in 1:length(config)-2]
	velocity = [ (config[i+1] - config[i])/delta_t for i in 1:length(config)-1 ]
	accel = [ (velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1 ]
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	g_of_t = first_coeff_scaled.*accel - sliced_noise
	return g_of_t,velocity,accel
end

function get_fractderiv(h,delta_t,steps,og_func,f0,fracderiv_steps,binom_part1,cap=1)
	fract_deriv = [0.0 for i in 1:steps-1]
	times = [i*delta_t for i in 1:steps-1]
	for i in 1:steps-1
		top = i
		if i >= 169
			top = 169
		end
		for j in 0:top#fracderiv_steps*i
			func_here = og_func[fracderiv_steps*i-j+1] - cap*f0
			#fact_j = exp(sum([log(i) for i in 1:j]))
			#binom_part = gamma(3-2*h)/(fact_j*gamma(3-2*h-j))
			binom_part = binom_part1[j+1]
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here/(delta_t/fracderiv_steps)^(2-2*h)
		end
	end
	return fract_deriv
end

function get_resids(h,config,delta_t,noise,noise_steps,lambda,gam,binom_part1)
	g_stuff = get_goft(h,config,delta_t,noise,noise_steps,lambda,gam)
	#println("Got G(t)")
	steps = length(config)
	final_time = steps*delta_t
	fract_deriv = [get_fractderiv(h,delta_t,steps,config,0.0,noise_steps,binom_part1)[i] for i in 2:steps-1]
	#println("Got Fract Deriv")
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	coeff = gam*scaled_coeffs[2]*(scaled_coeffs[1]^(2*h-2))
	return abs.(g_stuff[1]+coeff.*fract_deriv)#,g_stuff[1],coeff.*fract_deriv
end

function move_position(num_times,chosen,step_size)
	shift_matrix = [0.0 for i = 1:num_times]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,h,num_times,chosen,step_size,delta_t,noise,noise_steps,top_val,binom_part1)
	start_resids = get_resids(h,config,delta_t,noise,noise_steps,lambda,gam,binom_part1)
	shift_matrix = move_position(num_times,chosen,step_size)
	#println(shift_matrix,config)
	new_resids = get_resids(h,config+shift_matrix,delta_t,noise,noise_steps,lambda,gam,binom_part1)
	exp_diff = exp.(new_resids - start_resids)
	#println(length(config)-2,", ",length(exp_diff))
	#checking = [ exp_diff[i] <= top_val for i in 1:length(config)-2 ] #1.000001
	checking = exp_diff[chosen-2] <= top_val
	if checking#all(checking)
		return config+shift_matrix, 1#, new_resids, start_resids, shift_matrix
	else
		return config, 0, new_resids#, start_resids, shift_matrix
	end
	
	return "Acceptance Calculation Error"
end

function main_here(tol,steps,step_size,h,time_told,t_fin,lambda,gam,noise,noise_steps,x0,v0,top_val,fixed)
	#println("Starting")
	# getting first config from first guess
	running_config = append!([0.0,fixed],get_first_guess(h,time_told,t_fin,lambda,gam,noise,x0,v0,noise_steps)[2][3:time_told])
	# only save configuration data for every 10 attempted movements
	samp_freq = 10
	time_config = fill(0.0,(time_told,Int(steps/samp_freq)))
	time_resids = fill(0.0,(time_told-2,Int(steps/samp_freq)))
	index = 1
	index2 = 0
	delta_t = t_fin/time_told
	acc_rate = 0
	num_wrong = [i+2 for i in 1:time_told-2]
	factorial_top = time_told
	if time_told > 169
		factorial_top = 169
	end
	fact_j = [exp(sum([log(i) for i in 1:j])) for j in 0:factorial_top]
	binom_part1 = gamma(3-2*h).*[1/(fact_j[i]*gamma(3-2*h-(i-1))) for i in 1:length(fact_j)]
	for i in 1:steps
		# each MC time step every time point has attempted move, except starting point
		#=
		for k in 3:time_told
			movement = acc_rej_move(running_config,h,time_told,k,step_size,delta_t,noise,noise_steps,top_val)
			running_config = movement[1]
			acc_rate += movement[2]
			#println(movement[2],", ",movement[3],movement[4],movement[5])
		end
		=#
		
		upper = 4
		if length(num_wrong) < 4
			upper = length(num_wrong)
		end
		
		for k in num_wrong[1:upper]#3 + num_correct:upper
			movement = acc_rej_move(running_config,h,time_told,k,step_size,delta_t,noise,1,top_val,binom_part1)
			running_config = movement[1]
			acc_rate += movement[2]
			index2 += 1
		end
		num_wrong = []
		
		current_res = get_resids(h,running_config,delta_t,noise,noise_steps,lambda,gam,binom_part1)
		
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
			println("Solution Found in $i Steps")
			return running_config,time_config,time_resids,i
		end
		
		# interface data
		if i%(steps*0.01) == 0
			println("Running:"," ",100*i/steps,"%, ","Acceptance: ",acc_rate,"/",index2,", Number Wrong: ",length(num_wrong))
		end
	end
	
	println("No Solution")
	return time_config,time_resids
end

function get_sd(soln) # squared displacement
	time_steps = length(soln)
	start = soln[1]
	sd = [0.0 for i in 1:time_steps]
	for i in 1:time_steps
		sd[i] = abs2(soln[i] - start)
	end
	return sd
end



# boundary points b/c 2nd order SDE

final_time = 10000
time_steps = 10000
times = [i*final_time/time_steps for i in 0:time_steps-1]
lambda = 5.0
bets = 1.0
a = 1.0
gam = bets*lambda - a
a0 = 0.0
x0 = 0.0
v0 = 0.0
noise_steps = 1
tol = 0.0001
mc_steps = 1000000
metro_val = 1.000001
step_size = 0.001#0.008*final_time/time_steps
h = 0.5

msds = [[0.0 for i in 1:time_steps] for j in 0:0]
for i in 1:9
	white_noise = get_noise(0.5,noise_steps*time_steps,final_time,Int(2*i)).*fluc_dissp_coeffs("white",bets,lambda,a,h) 
	colored_noise = get_noise(h,noise_steps*time_steps,final_time,Int(2*i+1)).*fluc_dissp_coeffs("color",bets,lambda,a,h) 
	noise = white_noise - colored_noise#get_noise(h,noise_steps*time_steps,final_time,1).*fluc_dissp_coeffs("color",0,0,gam,h)
	exact = lang_soln(h,time_steps,noise_steps,noise,gam,lambda,final_time,v0,1)
	msds[1] += get_sd(exact[2])/9
end
#num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,gam,noise,noise_steps,x0,v0,metro_val,exact[2][2])

time_start = 10^2#2*10^(-1)
time_end = 4*10^2#2*10^0
time_start2 = 2.5*10^1
time_end2 = 10^2
times_squared = [time_start + i*(time_end-time_start)/10 for i in 0:10]
times_squared2 = [time_start2 + i*(time_end2-time_start2)/10 for i in 0:10]
yaxis_times_linear = [(times_squared2[i])*0.04 for i in 1:length(times_squared)] #[(times_squared[i]^2)*0.001 for i in 1:length(times_squared)]
yaxis_times_squared2 = [(times_squared2[i]^2)*(10^-1) for i in 1:length(times_squared2)]
#plot(times_squared,yaxis_times_squared,"-k")
plot(times_squared2,yaxis_times_squared2,"-k",label=latexstring("\$ t^2 \$"))
plot(times_squared,yaxis_times_linear,"-r",label=latexstring("\$ t \$"))
legend()

plot(times,msds[1])
xlabel("Time")
ylabel("MSD")
title("Mean Squared Deviation for No Memory")
xscale("log")
yscale("log")
#plot(times,num_soln[1],"-p")

#=
hs = [0.55,0.6,0.65,0.7,0.75]
avg_time = [0.0 for i in 1:length(hs)]
errors = [0.0 for i in 1:length(hs)]
for i in 1:length(hs)
	h = hs[i]
	local_avg = [0.0 for i in 1:25]
	for j in 2:6
		println(i,", ",j)
		noise = get_noise(h,noise_steps*time_steps,final_time,j).*fluc_dissp_coeffs("color",0,0,gam,h)
		exact = lang_soln(h,time_steps,noise_steps,noise,gam,lambda,final_time,v0,j)
		for k in 1:5
			num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,gam,noise,noise_steps,x0,v0,metro_val,exact[2][2])
			local_avg[i] = num_soln[4]
		end
	end
	avg_time[i] = mean(local_avg)
	errors[i] = std(local_avg)
end

errorbar(hs,avg_time,yerr=[errors,errors],fmt="-o")
xlabel("Hurst Parameter, H")
ylabel("MC Steps to Solution")
=#

	
		




#= Optimization of step size
step_start = 0.0001
step_end = 0.003
step_sizes = [(step_start + i*(step_end-step_start)/10) for i in 0:10]
avg_times = [0.0 for i in 1:length(step_sizes)]
errors = [0.0 for i in 1:length(step_sizes)]
for i in 1:length(step_sizes)
	step_size = step_sizes[i]*final_time/time_steps
	local_avg = [0.0 for i in 1:5]
	for j in 1:5
		println(i,", ",j)
		num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,gam,noise,noise_steps,x0,v0,metro_val,exact[2][2])
		local_avg[j] = num_soln[4]
	end
	println(mean(local_avg),", ",std(local_avg))
	avg_times[i] = mean(local_avg)
	errors[i] = std(local_avg)
end

errorbar(step_sizes,avg_times,yerr=[errors,errors],fmt="-o")
xlabel("Step Size (Percent of Time Step)")
ylabel("MC Steps to Solution")
title("Optimization of Max Step Size")
=#
#=
samp_freq = 10
which_steps = [1,100,500,740]
for i in 1:length(which_steps)
	chosens = which_steps[i]
	step_number = samp_freq*chosens
	plot(exact[1][2:end-1],num_soln[3][:,chosens],label="$step_number")
end
#plot(exact[1],num_soln[1],"-p",label="Num")
#plot(exact[1],exact[2],label="Exact")
xlabel("Time")
ylabel("Position")
title("Residuals of MC Simulation of Fractional Langevin Equation, H=$h")
legend()
=#


"fin"
