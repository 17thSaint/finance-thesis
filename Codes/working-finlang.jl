#import Pkg; Pkg.add("HDF5")
#import Pkg; Pkg.add("SpecialFunctions")
#import Pkg; Pkg.add("Statistics")
#import Pkg; Pkg.add("MittagLeffler")
using HDF5,SpecialFunctions,PyPlot,Statistics,MittagLeffler,LaTeXStrings
include("find-input.jl")

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
	scaled_time = (lambda/gam)^(0.5/h)
	scaled_pos = fluc_dissp_coeffs("color",0,0,gam,h)*scaled_time/sqrt(gam*lambda)
	return scaled_time,scaled_pos
end

function get_fractderiv(h,delta_t,steps,og_func,f0,noise_steps,binom_part1,cap=1)
	fract_deriv = [0.0 for i in 1:steps-1]
	times = [i*delta_t for i in 1:steps-1]
	for i in 1:steps-1
		top = i
		if i >= 169
			top = 169
		end
		for j in 0:top
			func_here = og_func[noise_steps*i-j+1] - cap*f0
			#fact_j = exp(sum([log(i) for i in 1:j]))
			#binom_part = gamma(3-2*h)/(fact_j*gamma(3-2*h-j))
			binom_part = binom_part1[j+1]
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
		#if i%(time_steps*0.01) == 0
		#	println("Getting First Guess: ",100*i/time_steps," %")
		#end
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

function get_resids(h,config,delta_t,noise,noise_steps,lambda,bets,a,binom_part1,upper_bound)
	g_stuff = get_goft(h,config,delta_t,noise,noise_steps,lambda,bets)
	#println("Got G(t)")
	steps = upper_bound#length(config)
	final_time = steps*delta_t
	frac_part = get_fractderiv(h,delta_t,steps,config,0.0,noise_steps,binom_part1)
	time_coeff_fractderiv = [1.0 for i in 2:steps-1]#[((i*delta_t)^(1-2*h)) for i in 2:steps-1]
	fract_deriv = [time_coeff_fractderiv[i-1]*frac_part[i] for i in 2:steps-1]
	#println("Got Fract Deriv")
	#scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	coeff = a#*gamma(2*h)#*scaled_coeffs[2]*(scaled_coeffs[1]^(2*h-2))
	return abs.(g_stuff[1][1:upper_bound-2]-coeff.*fract_deriv)#,g_stuff[1],coeff.*fract_deriv
end

function move_position(num_times,chosen,step_size)
	shift_matrix = [0.0 for i = 1:num_times]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,h,lambda,bets,a,num_times,chosen,step_size,delta_t,noise,noise_steps,top_val,binom_part1)
	upper_bound = chosen + 10
	if upper_bound > length(config)
		upper_bound = length(config)
	end
	start_resids = get_resids(h,config,delta_t,noise,noise_steps,lambda,bets,a,binom_part1,upper_bound)
	shift_matrix = move_position(num_times,chosen,step_size)
	#println(shift_matrix,config)
	new_resids = get_resids(h,config+shift_matrix,delta_t,noise,noise_steps,lambda,bets,a,binom_part1,upper_bound)
	exp_diff = exp.(new_resids - start_resids)
	#println(length(config)-2,", ",length(exp_diff))
	#checking = [ exp_diff[i] <= top_val for i in 1:length(config)-2 ] #1.000001
	checking = exp_diff[chosen-2] <= top_val
	if checking#all(checking)
		return config+shift_matrix, 1#, new_resids, start_resids, shift_matrix
	else
		return config, 0, exp_diff[chosen-2]#, start_resids, shift_matrix
	end
	
	return "Acceptance Calculation Error"
end

function main_here(tol,steps,step_size,h,time_told,t_fin,lambda,bets,a,noise,noise_steps,x0,v0,top_val,starting_input=[false,[0.0]])
	#println("Starting")
	# getting first config from first guess
	if starting_input[1]
		running_config = append!(starting_input[2],get_first_guess(h,time_told,t_fin,lambda,bets,noise,x0,v0,noise_steps)[2][length(starting_input[2])+1:time_told])
	else
		running_config = append!([0.0],get_first_guess(h,time_told,t_fin,lambda,bets,noise,x0,v0,noise_steps)[2][2:time_told])
	end
	# only save configuration data for every 10 attempted movements
	samp_freq = Int(0.01*steps)
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
		
		upper = 4
		if length(num_wrong) < 4
			upper = length(num_wrong)
		end
		
		for k in num_wrong[1:upper]
			movement = acc_rej_move(running_config,h,lambda,bets,a,time_told,k,step_size,delta_t,noise,1,top_val,binom_part1)
			running_config = movement[1]
			acc_rate += movement[2]
			index2 += 1
		end
		num_wrong = []
		
		current_res = get_resids(h,running_config,delta_t,noise,noise_steps,lambda,bets,a,binom_part1,time_told)
		
		#= saving data every 10 steps
		if i%samp_freq == 0
			time_config[:,index] = [running_config[x] for x in 1:time_told]
			time_resids[:,index] = [current_res[y] for y in 1:time_told-2]
			index += 1
		end
		=#
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
			println("Running:"," ",100*i/steps,"%, ","Acceptance: ",acc_rate,"/",index2,", Number Wrong: ",length(num_wrong),", H = $h")
			#plot(running_config)
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

function get_alpha_coeff(alpha)
	num = gamma((1+alpha)/2)*gamma((3-alpha))
	denom = gamma(2-alpha)*gamma(alpha)
	result = sqrt(num/denom)
	return result
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

function write_data_hdf5(ver,data,h,final_time,time_count,lambda,bets,a,avg_count=1,which=1)
	cd("..")
	cd("fin-data")
	println("Starting Data Write")
	if ver == "msd"
		binary_file_pos = h5open("$ver-h-$h-ft-$final_time-points-$time_count-lam-$lambda-bet-$bets-a-$a-errors-v0-avg-$avg_count.hdf5","w")
		create_group(binary_file_pos,"all-data")
		alldata = binary_file_pos["all-data"]
		alldata["msd"] = data[1]
		alldata["errors"] = data[2]
	else
		binary_file_pos = h5open("$ver-h-$h-ft-$final_time-points-$time_count-lam-$lambda-bet-$bets-a-$a-which-$which.hdf5","w")
		create_group(binary_file_pos,"all-data")
		alldata = binary_file_pos["all-data"]
		alldata["path"] = data
	end
	close(binary_file_pos)
	cd("..")
	cd("Codes")
	println("Data Added, File Closed: $h")
end

function read_msd_hdf5_data(ver,h,final_time,time_count,lambda,bets,a,avg_count=1,which=1)
	cd("..")
	cd("fin-data")
	if ver == "msd"
		file = h5open("$ver-h-$h-ft-$final_time-points-$time_count-lam-$lambda-bet-$bets-a-$a-errors-v0-avg-$avg_count.hdf5","r")
		data = [read(file["all-data"],ver),read(file["all-data"],"errors")]
	else
		file = h5open("$ver-h-$h-ft-$final_time-points-$time_count-lam-$lambda-bet-$bets-a-$a-which-$which.hdf5","r")
		data = read(file["all-data"],ver)
	end
	cd("..")
	cd("Codes")
	return data
end



final_time = 500
time_steps = 500
times = [i*final_time/time_steps for i in 0:time_steps-1]
a = 1.0
a0 = 0.0
x0 = 0.0
lambda = 50.0
bets = 1.0
kbt = 1.0
v0 = sqrt(kbt/lambda)
noise_steps = 1

tol = 0.01
mc_steps = 5000000
metro_val = 1.000001
step_size = 0.1#0.00001
#h = 0.75
#alph = 2-2*h


#= expected value at given time
selected_time = time_steps
white_count = 7
#lambdas = [1.0 + i*10/10 for i in 0:10]
as = [1.5*bets*lambda + i*5.5*bets*lambda/10 for i in 0:10]
as_rel = [as[i]/(bets*lambda) for i in 1:length(as)]
hs = [0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975]
counts_hs = [20,20,10,20,10,20,10,20,10,20,10,20,10,10,10,10,10,10,9]

expected_price = [[0.0 for i in 1:length(as)] for j in 1:length(hs)]
errors = [[0.0 for i in 1:length(as)] for j in 1:length(hs)]
for l in 1:length(hs)
	h = hs[l]
	count = counts_hs[l]
	for i in 1:1#length(as)
		println(l,", ",i)
		#lambda = lambdas[i]
		a = as[3]#as[i]
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
	#total_counts = count*white_count
	#errorbar(lambdas,expected_price[l],yerr=[errors[l],errors[l]],label="$total_counts")
end


for i in 1:1#length(hs)
	h = hs[i]
	plot(hs,[errors[i][1] for i in 1:length(hs)])
end
#legend()

#ylabel(latexstring("\$ \\mathbb{E} \$ [P(t=1)]"))
ylabel("Volatility")
#title("Expected Value of Price at t = 1 vs Liquidity")
title("Volatility vs Hurst Parameter")
#xlabel(latexstring("Market Liquidity, \$ \\lambda \$"))
#xlabel(latexstring("\$ a / \\beta\\lambda\$"))
xlabel("Hurst parameter, H")

=#
# volat_vs_h_lambda1.0 = [0.03532629342614326,0.028564641189068728,0.03525452953354956,0.025626495202568073,0.026894175157000696,0.025600327767250124,0.02628260457135696,0.027594002666838115,0.028244514978930213,0.024490799420920556,0.03496750211959252,0.027226628326523727,0.02666934331336109,0.022956596077158116,0.018288178038532128,0.02969355060967904,0.022257833983462684,0.022508391541704496,0.019910118360160132]
# hs = hs = [0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975]

# errors = [[0.011123156136320799, 0.005113516536587368, 0.003229451912248103, 0.0023643448165093068, 0.0018743155764085202, 0.0015401585147341167, 0.0013137741339704504, 0.0011389413702820618, 0.001010775957045089, 0.0009048780704512181, 0.0008234191901862014],[0.010690379373428099, 0.005397942836544762, 0.0037371227273192815, 0.0028384420253396688, 0.002243368887067092, 0.0018625804539577845, 0.001584797556060686, 0.0013813023255170999, 0.0012193851039699265, 0.0010958038260582461, 0.000994683473608056],[0.014114381638318695, 0.007186096756963081, 0.004661006566537109, 0.0034527772937272797, 0.0027357279004645934, 0.0022620764057069624, 0.0019287586043179414, 0.0016841306345919099, 0.0014978804670156202, 0.0013424575102295678, 0.0012196607255481897]]
# expected_price = [[0.020083749435273208, 0.011822895704126878, 0.009161463126122284, 0.0077296992383412534, 0.0068070261296105895, 0.00614618526352113, 0.005650622970025069, 0.005248321812828266, 0.004931353414387341, 0.004659398794969061, 0.0044310333860214015], [0.016375740156689002, 0.009920234258719204, 0.007817103162634124, 0.0067191781031230345, 0.006019442233483335, 0.00549891647905151, 0.005101905307731918, 0.004775247323544507, 0.004507313034934853, 0.004282254473072314, 0.00408722643500412], [0.02037188071247686, 0.011972135885754893, 0.009333159912950674, 0.007911678668445245, 0.006960053662540694, 0.006282737274525685, 0.005774241292109119, 0.005365755113685405, 0.005036202348615111, 0.004758759318343516, 0.004516900220451618]]



# Mean squared displacement
hs = [0.95,0.75,0.55]
#lambdas = [1.0,2.0,5.0]
counts_hs = [10 for i in 1:length(hs)]
#=
lambda = 5.0#[0.1 + i*1.9/20 for i in 0:19]
as = [0.1 + i*2.8/20 for i in 0:19]#[0.1 + i*4.8/20 for i in 0:19]#1.0 .* lambdas
betss = [0.1 + i*0.5/20 for i in 0:19]
alphas = sort([2-2*hs[i] for i in 1:length(hs)])

alphas_lwbig_tcbig = []
as_lwbig_tcbig = []
betss_lwbig_tcbig = []

alphas_lwbig_twbig = []
as_lwbig_twbig = []
betss_lwbig_twbig = []

alphas_lwbig_teq = []
as_lwbig_teq = []
betss_lwbig_teq = []


alphas_lcbig_tcbig = []
as_lcbig_tcbig = []
betss_lcbig_tcbig = []

alphas_lcbig_twbig = []
as_lcbig_twbig = []
betss_lcbig_twbig = []

alphas_lcbig_teq = []
as_lcbig_teq = []
betss_lcbig_teq = []

for k in 1:length(betss)
	bets = betss[k]
	tw = 1/bets
	for j in 1:length(alphas)
		alpha = alphas[j]
		for i in 1:length(as)
			a = as[i]
			tc = lambda/a
			lw = sqrt(2)/(bets*sqrt(lambda))
			lc = (sqrt(2*lambda)/a)*get_alpha_coeff(alpha)
			if lw > lc
				if tw > tc
					append!(alphas_lwbig_twbig,[alpha])
					append!(as_lwbig_twbig,[a])
					append!(betss_lwbig_twbig,[bets])
				elseif tw < tc
					append!(alphas_lwbig_tcbig,[alpha])
					append!(as_lwbig_tcbig,[a])
					append!(betss_lwbig_tcbig,[bets])
				elseif tw == tc
					append!(alphas_lwbig_teq,[alpha])
					append!(as_lwbig_teq,[a])
					append!(betss_lwbig_teq,[bets])
				end
			elseif lw < lc
				if tw > tc
					println(alpha,", ",a,", ",bets)
					append!(alphas_lcbig_twbig,[alpha])
					append!(as_lcbig_twbig,[a])
					append!(betss_lcbig_twbig,[bets])
				elseif tw < tc
					append!(alphas_lcbig_tcbig,[alpha])
					append!(as_lcbig_tcbig,[a])
					append!(betss_lcbig_tcbig,[bets])
				elseif tw == tc
					append!(alphas_lcbig_teq,[alpha])
					append!(as_lcbig_teq,[a])
					append!(betss_lcbig_teq,[bets])
				end
			end
		end
	end
end
#=
scatter(alphas_lwbig_twbig,betss_lwbig_twbig,c="b",label="lw>lc,tw>tc")
scatter(alphas_lcbig_twbig,betss_lcbig_twbig,c="r",label="lw<lc,tw>tc")
scatter(alphas_lwbig_tcbig,betss_lwbig_tcbig,c="g",label="lw>lc,tw<tc")
scatter(alphas_lcbig_tcbig,betss_lcbig_tcbig,c="c",label="lw<lc,tw<tc")
scatter(alphas_lwbig_teq,betss_lwbig_teq,c="k",label="lw>lc,tw=tc")
scatter(alphas_lcbig_teq,betss_lcbig_teq,c="m",label="lw<lc,tw=tc")
=#
#
scatter3D(as_lwbig_twbig,alphas_lwbig_twbig,betss_lwbig_twbig,c="b",label="lw>lc,tw>tc")
scatter3D(as_lcbig_twbig,alphas_lcbig_twbig,betss_lcbig_twbig,c="r",label="lw<lc,tw>tc")
scatter3D(as_lwbig_tcbig,alphas_lwbig_tcbig,betss_lwbig_tcbig,c="g",label="lw>lc,tw<tc")
scatter3D(as_lcbig_tcbig,alphas_lcbig_tcbig,betss_lcbig_tcbig,c="c",label="lw<lc,tw<tc")
scatter3D(as_lwbig_teq,alphas_lwbig_teq,betss_lwbig_teq,c="k",label="lw>lc,tw=tc")
scatter3D(as_lcbig_teq,alphas_lcbig_teq,betss_lcbig_teq,c="m",label="lw<lc,tw=tc")
#
xlabel("A")
ylabel("Alpha")
zlabel("Beta")
legend()
=#
#selected = parse(Int64,ARGS[1])
#h = hs[selected]
#as = [0.1,1.0,2.5,5.0]

#
msds = [[0.0 for j in 1:time_steps] for i in 1:length(hs)]
errors = [[0.0 for j in 1:time_steps] for i in 1:length(hs)]
for i in 1:length(hs)
	h = hs[i]
	count = counts_hs[i]
	local_msds = fill(0.0,(time_steps,count))
	for j in 1:count
		white_noise = get_noise(0.5,time_steps,final_time,j).*fluc_dissp_coeffs("white",bets,lambda,a,0.5)
		colored_noise = get_noise(h,time_steps,final_time,j).*fluc_dissp_coeffs("color",bets,lambda,a,h)
		noise = white_noise - colored_noise
		#
		cd("..")
		cd("fin-data")
		check_ifinput = get_input_path(h,lambda,bets,a,j)
		cd("..")
		cd("Codes")
		if check_ifinput[1]
			og_dats = read_msd_hdf5_data("path",h,check_ifinput[3],check_ifinput[2],lambda,bets,a,1,j)
		else
			og_dats = [0.0]
		end
		starting_input = [check_ifinput[1],og_dats]
		#
		num_soln = main_here(tol,mc_steps,step_size,h,time_steps,final_time,lambda,bets,a,noise,1,x0,v0,metro_val)
		local_msds[:,j] = get_sd(num_soln[1])
		#plot(times,num_soln[1],label="$h")
		write_data_hdf5("path",num_soln[1],h,final_time,time_steps,lambda,bets,a,1,j)
	end
	#
	msds[i] = [mean(local_msds[k,:]) for k in 1:time_steps]
	errors[i] = [std(local_msds[k,:]) for k in 1:time_steps]
	write_data_hdf5("msd",[msds[i],errors[i]],h,final_time,time_steps,lambda,bets,a,count)
	alph = round(2-2*h,digits=2)
	#errorbar(times,msds[i],yerr=[errors[i],errors[i]],label="$alph")
	plot(times,msds[i],label="$alph")
	xscale("log")
	yscale("log")
	#
end
legend()
#

#plot(times,msds[1])

#= including comparison lines for MSD
for i in 4:6
	which_bets = 2
	tcomp = "tw < tc"
	if which_bets == 1
		tcomp = "tw > tc"
	end
	bets = betss[which_bets]
	lcomp = "lw < lc"
	if i > 3
		lcomp = "lw > lc"
	end
	h = hs[i]
	alph = round(2-2*h,digits=2)
	dats = read_msd_hdf5_data("msd",h,final_time,time_steps,lambda,bets,a,5)[1]
	plot(times,dats,label="$alph")
	title("MSD for $tcomp and $lcomp")
	#title(latexstring("Mean Squared Deviation \$ \\alpha = $alph \$ with \$ v_0 \$"))
end
=#
#=
xscale("log")
yscale("log")
time_start_cubed_short = 4*10^(-1)
time_end_cubed_short = 6*10^(-1)
time_start_squared_short = 2*10^(-1)
time_end_squared_short = 5*10^(-1)
times_squared_shorttime = [time_start_squared_short + i*(time_end_squared_short-time_start_squared_short)/10 for i in 0:10]
times_cubed_shorttime = [time_start_cubed_short + i*(time_end_cubed_short-time_start_cubed_short)/10 for i in 0:10]
yaxis_times_squared_shorttime = [(times_squared_shorttime[i]^2)*(10^(-0.5)) for i in 1:length(times_squared_shorttime)]
yaxis_times_cubed_shorttime = [(times_cubed_shorttime[i]^3)*(10^(-1.0)) for i in 1:length(times_cubed_shorttime)]
plot(times_squared_shorttime,yaxis_times_squared_shorttime,"-k",label=latexstring("\$ t^2 \$"))
plot(times_cubed_shorttime,yaxis_times_cubed_shorttime,"-r",label=latexstring("\$ t^3 \$"))
legend()
=#
#=
hs = [0.95,0.75,0.55]
h = hs[3]
alph = round(2-2*h,digits=2)
dats = read_msd_hdf5_data("msd",h,final_time,time_steps,lambda,bets,a,5)[1]
plot(times,dats)#,label="$alph")
title(latexstring("Mean Squared Deviation \$ \\alpha = $alph \$ with \$ v_0 \$"))
=#
#=
time_start_longtime = 1*10^(1)#2.5*10^(1)
time_end_longtime = 4*10^(1)#10^(2)
times_squared_longtime = [time_start_longtime + i*(time_end_longtime-time_start_longtime)/10 for i in 0:10]
#alph = 0
#yaxis_times_linear = [(times_squared2[i])*0.04 for i in 1:length(times_squared)] #[(times_squared[i]^2)*0.001 for i in 1:length(times_squared)]
#yaxis_times_2malph_longtime = [(times_squared_longtime[i]^(2-alph))*(10^-0.5) for i in 1:length(times_squared_longtime)]
#yaxis_times_2palph_longtime = [(times_squared_longtime[i]^(2+alph))*(10^1.7) for i in 1:length(times_squared_longtime)]
yaxis_times_squared_longtime = [(times_squared_longtime[i]^(2))*(10^0.3) for i in 1:length(times_squared_longtime)]
#plot(times_squared_longtime,yaxis_times_2malph_longtime,"-r",label=latexstring("\$ t^{2-\\alpha} \$"))
plot(times_squared_longtime,yaxis_times_squared_longtime,"-k",label=latexstring("\$ t^{2} \$"))
#plot(times_squared_longtime,yaxis_times_2palph_longtime,"-r",label=latexstring("\$ t^{2+\\alpha} \$"))
#plot(times_squared_longtime,yaxis_times_squared_longtime,"-k",label=latexstring("\$ t^{2} \$"))

legend()
xscale("log")
yscale("log")
xlabel("Time")
ylabel("MSD")
#title(latexstring("Mean Squared Deviation for range \$ \\alpha \$"))
=#

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
