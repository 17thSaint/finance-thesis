using HDF5,SpecialFunctions,PyPlot

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

function get_first_guess(h,time_told,t_fin,lambda,gam,noise,x0,v0,a0,noise_steps)
	time_steps = time_told-1
	second_term = [0.0 for i in 1:time_steps]
	first_term = [0.0 for i in 1:time_steps]  # convolution part
	times = [0.0 for i in 1:time_steps]
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	for i in 1:time_steps
		#println("Getting First Guess: ",100*i/time_steps," %")
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
function get_goft(h,config,delta_t,noise,lambda,gam,x0,v0,a0)
	println("Geting g(t)")
	sliced_noise = [ noise[i] for i in 2:length(config)-1]
	println("Got Noise")
	velocity = [ (config[i+1] - config[i])/delta_t for i in 1:length(config)-1 ]
	println("Got Velocity")
	accel = [ (velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1 ]
	#accel[2] *= 2
	println("Got Acceleration")
	scaled_coeffs = get_scale_inv_vals(h,lambda,gam)
	first_coeff_scaled = lambda*scaled_coeffs[2]/(scaled_coeffs[1]^2)
	g_of_t = first_coeff_scaled.*accel - sliced_noise
	return g_of_t,velocity,accel
end


final_time = 100
time_steps = 1000
lambda = 1.0
gam = 1.0
a0 = 0.0
x0 = 0.0
v0 = 0.0
noise_steps = 1
h = 0.3

noise = get_noise(h,noise_steps*time_steps,final_time,1).*fluc_dissp_coeffs("color",0,0,gam,h)
check_first = get_first_guess(h,time_steps,final_time,lambda,gam,noise,x0,v0,a0,noise_steps)
resids = get_goft(h,check_first[2],final_time/time_steps,noise,lambda,gam,x0,v0,a0)
plot(resids[1])



"fin"
