#import Pkg; Pkg.add("HDF5")
#import Pkg; Pkg.add("SpecialFunctions")
#import Pkg; Pkg.add("Statistics")
#import Pkg; Pkg.add("MittagLeffler")
using HDF5,SpecialFunctions,PyPlot,Statistics,MittagLeffler,LaTeXStrings


function fluc_dissp_coeffs(which,bets,lambda,a,h,kBT=1,sigma0=1)
	if which == "white"
		whi_coeff = 2*bets*lambda*kBT
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

function get_fractderiv(h,delta_t,steps,og_func,f0,noise_steps,cap=1)
	fract_deriv = [0.0 for i in 1:steps-1]
	times = [i*delta_t for i in 1:steps-1]
	for i in 1:steps-1
		for j in 0:noise_steps*i
			func_here = og_func[noise_steps*i-j+1] - cap*f0
			fact_j = exp(sum([log(i) for i in 1:j]))
			binom_part = gamma(3-2*h)/(fact_j*gamma(3-2*h-j))
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here/(delta_t/noise_steps)^(2-2*h)
		end
	end

	return fract_deriv
end


function get_first_guess(h,time_told,t_fin,lambda,bets,noise,x0,v0,noise_steps)
	time_steps = time_told-1
	second_term = [0.0 for i in 1:time_steps]
	first_term = [0.0 for i in 1:time_steps]  # convolution part
	times = [0.0 for i in 1:time_steps]
	sliced_noise = [ noise[i*noise_steps+1] for i in 0:time_steps]
	noise_steps_local = 1
	#scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	#first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	for i in 1:time_steps
		if i%(time_steps*0.01) == 0
			println("Getting First Guess: ",100*i/time_steps," %")
		end
		time_now = i*t_fin/time_told
		times[i] = time_now
		second_term[i] = v0*(1-exp(-bets*time_now)) + x0
		noise_times = append!([0.0],[ l*t_fin/time_told + k*t_fin/(time_told*noise_steps_local) for l in 0:i-1 for k in 1:noise_steps_local])
	
		for j in 1:length(noise_times)-1
			left_time = noise_times[j]
			right_time = noise_times[j+1]
			
			noises_right = sliced_noise[j+1]
			noises_left = sliced_noise[j]
			
			#needs scaling term
			exp_part_left = (1-exp(-bets*(time_now - left_time)))/(bets*lambda)  
			exp_part_right = (1-exp(-bets*(time_now - right_time)))/(bets*lambda)
			
			first_term[i] += 0.5*(right_time-left_time)*(noises_left*exp_part_left + noises_right*exp_part_right)
		end

	end
	full_rez = append!([x0],first_term + second_term)
	times = append!([0.0],times)
	return times,full_rez
end

# g(t) is the SDE without the fractional derivative term
function get_goft(h,config,delta_t,noise,noise_steps,lambda,bets)
	sliced_noise = [ noise[i*noise_steps+1] for i in 1:length(config)-2]
	velocity = [ (config[i+1] - config[i])/delta_t for i in 1:length(config)-1 ]
	accel = [ (velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1 ]
	#scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	#first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	g_of_t = lambda.*accel - sliced_noise + (bets*lambda).*(velocity[1:end-1])
	return g_of_t,velocity,accel,sliced_noise
end

function get_resids(h,config,delta_t,noise,noise_steps,lambda,bets,a)
	g_stuff = get_goft(h,config,delta_t,noise,noise_steps,lambda,bets)
	#println("Got G(t)")
	steps = length(config)
	final_time = steps*delta_t
	fract_deriv = [((i*delta_t)^(1-2*h))*get_fractderiv(h,delta_t,steps,config,0.0,noise_steps)[i] for i in 2:steps-1]
	#println("Got Fract Deriv")
	#scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	coeff = a*gamma(2*h)#*scaled_coeffs[2]*(scaled_coeffs[1]^(2*h-2))
	return abs.(g_stuff[1]-coeff.*fract_deriv)#,g_stuff[1],coeff.*fract_deriv
end

function move_position(num_times,chosen,step_size)
	shift_matrix = [0.0 for i = 1:num_times]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,h,lambda,bets,a,num_times,chosen,step_size,delta_t,noise,noise_steps,top_val)
	start_resids = get_resids(h,config,delta_t,noise,noise_steps,lambda,bets,a)
	shift_matrix = move_position(num_times,chosen,step_size)
	#println(shift_matrix,config)
	new_resids = get_resids(h,config+shift_matrix,delta_t,noise,noise_steps,lambda,bets,a)
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

function main_here(tol,steps,step_size,h,time_told,t_fin,lambda,bets,a,noise,noise_steps,x0,v0,top_val)
	#println("Starting")
	# getting first config from first guess
	running_config = append!([0.0],get_first_guess(h,time_told,t_fin,lambda,bets,noise,x0,v0,noise_steps)[2][2:time_told])
	# only save configuration data for every 10 attempted movements
	samp_freq = Int(0.01*steps)
	time_config = fill(0.0,(time_told,Int(steps/samp_freq)))
	time_resids = fill(0.0,(time_told-2,Int(steps/samp_freq)))
	index = 1
	delta_t = t_fin/time_told
	acc_rate = 0
	num_wrong = [i+2 for i in 1:time_told-2]
	for i in 1:steps
		
		upper = 4
		if length(num_wrong) < 4
			upper = length(num_wrong)
		end
		
		for k in num_wrong[1:upper]
			movement = acc_rej_move(running_config,h,lambda,bets,a,time_told,k,step_size,delta_t,noise,1,top_val)
			running_config = movement[1]
			acc_rate += movement[2]
		end
		num_wrong = []
		
		current_res = get_resids(h,running_config,delta_t,noise,noise_steps,lambda,bets,a)
		
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
		if i%(steps*0.001) == 0
			println("Running:"," ",100*i/steps,"%, ","Acceptance: ",acc_rate,", Number Wrong: ",length(num_wrong),", H = $h")
		end
	end
	
	println("No Solution")
	return time_config,time_resids
end

function get_avg_price_jump(series)
	len = length(series)
	avg_jump = [0.0 for i in 1:len-1]
	for i in 1:len-1
		avg_jump[i] = abs(series[i+1]-series[i])
	end
	return mean(avg_jump),std(avg_jump)
end

function auto_correlation(energies, delta_t)
    average_energy = mean(energies)
    
    points = Int(floor(length(energies)-delta_t))
    
    energy_fluctuations = energies.-average_energy
    
    autocorrelation_top = [0.0 for i in 1:points]
    autocorrelation_bottom = mean(energy_fluctuations.^2)
    error_bot = std(energy_fluctuations.^2)
    
    for i in 1:points
        autocorrelation_top[i] = energy_fluctuations[i]*energy_fluctuations[i+delta_t]
    end
    
    full_top = mean(autocorrelation_top)
    rezz = full_top/autocorrelation_bottom
    
    error_top = std(autocorrelation_top)
    full_error = rezz*sqrt((error_top/full_top)^2 + (error_bot/autocorrelation_bottom)^2)
    #println(full_error)
    return rezz,full_error
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

function write_msd_data_hdf5(data,h)
	println("Starting Data Write")
	binary_file_pos = h5open("msd-vs-h-$h.hdf5","w")
	create_group(binary_file_pos,"all-data")
	alldata = binary_file_pos["all-data"]
	alldata["msd"] = data
	close(binary_file_pos)
	println("Data Added, File Closed: $h")
end

function read_msd_hdf5_data(h)
	file = h5open("msd-vs-h-$h.hdf5","r")
	data = read(file["all-data"],"msd")
	return data
end



final_time = 1
time_steps = 10
times = [i*final_time/time_steps for i in 0:time_steps-1]
a = 1.0
a0 = 0.0
x0 = 0.0
v0 = 0.0
noise_steps = 1

tol = 0.001
mc_steps = 5000000
metro_val = 1.000001
step_size = 0.0001#0.005
lambda = 5.0
#h = 0.55
bets = 1.0


#= expected value at given time
selected_time = time_steps
white_count = 1
lambdas = [1.0 + i*10/10 for i in 0:10]
hs = [0.55,0.75,0.95]
counts_hs = [20,20,10]
#=
expected_price = [[0.0 for i in 1:length(lambdas)] for j in 1:length(hs)]
errors = [[0.0 for i in 1:length(lambdas)] for j in 1:length(hs)]
for l in 1:length(hs)
	h = hs[l]
	count = counts_hs[l]
	for i in 1:length(lambdas)
		println(l,", ",i)
		lambda = lambdas[i]
		local_expected_prices = []
		for j in 1:count
			colored_noise = get_noise(h,time_steps,final_time,j).*fluc_dissp_coeffs("color",bets,lambda,a,h)
			for k in 1:white_count
				white_noise = get_noise(0.5,time_steps,final_time,k).*fluc_dissp_coeffs("white",bets,lambda,a,0.5)
				noise = white_noise + colored_noise
				num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,bets,a,noise,1,x0,v0,metro_val)
				append!(local_expected_prices,abs(num_soln[1][end]))
			end
		end
		expected_price[l][i] = mean(local_expected_prices)
		errors[l][i] = std(local_expected_prices)
	end
	errorbar(lambdas,expected_price[l],yerr=[errors[l],errors[l]],label="$h")
end
legend()
=#
for i in 1:length(hs)
	h = hs[i]
	plot(lambdas,errors[i],label="$h")
end
#ylabel(latexstring("\$ \\mathbb{E} \$ [P(t=1)]"))
#title("Expected Value of Price at t = 1 vs Liquidity")
title("Variance of Expected Value vs Liquidity, range H")
xlabel(latexstring("Market Liquidity, \$ \\lambda \$"))


# errors = [[0.011123156136320799, 0.005113516536587368, 0.003229451912248103, 0.0023643448165093068, 0.0018743155764085202, 0.0015401585147341167, 0.0013137741339704504, 0.0011389413702820618, 0.001010775957045089, 0.0009048780704512181, 0.0008234191901862014],[0.010690379373428099, 0.005397942836544762, 0.0037371227273192815, 0.0028384420253396688, 0.002243368887067092, 0.0018625804539577845, 0.001584797556060686, 0.0013813023255170999, 0.0012193851039699265, 0.0010958038260582461, 0.000994683473608056],[0.014114381638318695, 0.007186096756963081, 0.004661006566537109, 0.0034527772937272797, 0.0027357279004645934, 0.0022620764057069624, 0.0019287586043179414, 0.0016841306345919099, 0.0014978804670156202, 0.0013424575102295678, 0.0012196607255481897]]
# expected_price = [[0.020083749435273208, 0.011822895704126878, 0.009161463126122284, 0.0077296992383412534, 0.0068070261296105895, 0.00614618526352113, 0.005650622970025069, 0.005248321812828266, 0.004931353414387341, 0.004659398794969061, 0.0044310333860214015], [0.016375740156689002, 0.009920234258719204, 0.007817103162634124, 0.0067191781031230345, 0.006019442233483335, 0.00549891647905151, 0.005101905307731918, 0.004775247323544507, 0.004507313034934853, 0.004282254473072314, 0.00408722643500412], [0.02037188071247686, 0.011972135885754893, 0.009333159912950674, 0.007911678668445245, 0.006960053662540694, 0.006282737274525685, 0.005774241292109119, 0.005365755113685405, 0.005036202348615111, 0.004758759318343516, 0.004516900220451618]]

=#

# Mean squared displacement
hs = [0.55,0.75,0.95]
#lambdas = [1.0,2.0,5.0]
counts_hs = [3,3,3]
selected = parse(Int64,ARGS[1])
msds = [[0.0 for j in 1:time_steps] for i in 1:length(hs)]
white_noise = get_noise(0.5,time_steps,final_time,1).*fluc_dissp_coeffs("white",bets,lambda,a,0.5)
h = hs[selected]
count = counts_hs[selected]
for j in 1:count
	colored_noise = get_noise(h,time_steps,final_time,j).*fluc_dissp_coeffs("color",bets,lambda,a,h)
	noise = white_noise + colored_noise
	num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,bets,a,noise,1,x0,v0,metro_val)
	msds[selected] += get_sd(num_soln[1])/count
end
write_msd_data_hdf5(msds[selected],h)

#= Price correlation to history
hs = [0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.975]
white_noise = get_noise(0.5,time_steps,final_time,1).*fluc_dissp_coeffs("white",bets,lambda,a,0.5)

dts = [i for i in 1:4]
selected_corrs = [[0.0 for i in 1:length(hs)] for j in 1:4]
errors = [[0.0 for i in 1:length(hs)] for j in 1:4]
for i in 1:length(hs)
	h = hs[i]
	colored_noise = get_noise(h,time_steps,final_time,3).*fluc_dissp_coeffs("color",bets,lambda,a,h)
	noise = white_noise + colored_noise
	num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,bets,a,noise,1,x0,v0,metro_val)
	for j in 1:4
		rezzs = auto_correlation(num_soln[1],dts[j])
		selected_corrs[j][i] = rezzs[1]
		errors[j][i] = rezzs[2]
	end
end
for i in 1:4
	sep = dts[i]
	errorbar(hs,selected_corrs[i],yerr=[errors[i],errors[i]],fmt="-o",label="$sep")
end
legend()
xlabel("Hurst parameter")
ylabel("History Correlation")
=#
#=  Volatility vs Hurst, no relation
hs = [0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.975]
white_noise = get_noise(0.5,time_steps,final_time,1).*fluc_dissp_coeffs("white",bets,lambda,a,0.5)
avg_jumps = [[0.0 for i in 1:length(hs)] for j in 1:3]
errors = [[0.0 for i in 1:length(hs)] for j in 1:3]
for j in 1:3
	for i in 1:length(hs)
		h = hs[i]
		colored_noise = get_noise(h,time_steps,final_time,j+1).*fluc_dissp_coeffs("color",bets,lambda,a,h)
		noise = white_noise + colored_noise
		num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,bets,a,noise,1,x0,v0,metro_val)
		dats = get_avg_price_jump(num_soln[1])
		avg_jumps[j][i] = dats[1]
		errors[j][i] = dats[2]
	end
	errorbar(hs,avg_jumps[j],yerr=[errors[j],errors[j]],fmt="-o")
end
=#
#= Volatility vs frequency of transactions
lambdas = [1.0,2.0,5.0]
hs = [0.55,0.75,0.95]
betss = [0.01 + i*10/10 for i in 0:9]
avg_jumps = [[0.0 for i in 1:length(betss)] for j in 1:length(lambdas)]
errors = [[0.0 for i in 1:length(betss)] for j in 1:length(lambdas)]
for j in 1:length(lambdas)
	#lambda = lambdas[j] 
	h = hs[j]
	for i in 1:length(betss)
		println(j,", ",i)
		bets = betss[i]
		white_noise = get_noise(0.5,time_steps,final_time,1).*fluc_dissp_coeffs("white",bets,lambda,a,0.5)
		colored_noise = get_noise(h,time_steps,final_time,2).*fluc_dissp_coeffs("color",bets,lambda,a,h)
		noise = white_noise + colored_noise
		num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,bets,a,noise,1,x0,v0,metro_val)
		rezzs = get_avg_price_jump(num_soln[1])
		avg_jumps[j][i] = rezzs[1]
		errors[j][i] = rezzs[2]
	end
	errorbar(betss,avg_jumps[j],yerr=[errors[j],errors[j]],fmt="-o",label="$h")
end
legend()
xlabel("Frequency of Transactions")
ylabel("AVG Price Jump")
title("Volatility vs Transaction Frequency for range Liquidity")
=#

#= Volat and Liquidity
hs = [0.7,0.725,0.775,0.8,0.55,0.675,0.75,0.875,0.975]
lambdas = [10 + i*80/10 for i in 0:10]
avg_jumps = [[0.0 for i in 1:length(lambdas)] for j in 1:length(hs)]
errors = [[0.0 for i in 1:length(lambdas)] for j in 1:length(hs)]
for j in 1:length(hs)
	h = hs[j]
	for i in 1:length(lambdas)
		println(j,", ",i)
		lambda = lambdas[i]
		white_noise = get_noise(0.5,time_steps,final_time,1).*fluc_dissp_coeffs("white",bets,lambda,a,0.5)
		data_here = [0.0 for p in 1:20]
		for k in 1:10
			colored_noise = get_noise(h,time_steps,final_time,k).*fluc_dissp_coeffs("color",bets,lambda,a,h)
			noise = white_noise + colored_noise
			num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,bets,a,noise,1,x0,v0,metro_val)
			data_here[k] = get_avg_price_jump(num_soln[1])[1]
		end
		avg_jumps[j][i] = mean(data_here)
		errors[j][i] = std(data_here)
		
	end
	errorbar(lambdas,avg_jumps[j],yerr=[errors[j],errors[j]],fmt="-o",label="$h")
end
legend()
xlabel("Market Liquidity")
ylabel("AVG Price Movement")
title("Volatility vs Market Liquidity")
=#
#= Paths with ranging memory strength
h = 0.75
as = [0.001,2.5,5.0]
	for i in 1:length(as)
		println(i)
		a = as[i]
		white_noise = get_noise(0.5,time_steps,final_time,1).*fluc_dissp_coeffs("white",bets,lambda,a,0.5)
		colored_noise = get_noise(h,time_steps,final_time,2).*fluc_dissp_coeffs("color",bets,lambda,a,h)
		noise = white_noise + colored_noise
		num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,bets,a,noise,1,x0,v0,metro_val)
		plot(times,num_soln[1],"-p",label="$a")
	end
legend()
xlabel("Time")
ylabel("Price")
title("Paths for range of Memory Strength")
=#


"fin"
