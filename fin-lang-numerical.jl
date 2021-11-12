include("simulate-frac-brownian.jl")

function fluc_dissp_coeffs(which,beta,lambda,a,h,kBT=1,sigma0=0.55804)
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

function get_first_guess(h,time_steps,t_fin,beta,lambda,a,white_noise,colored_noise,x0,v0)
	second_term = [0.0 for i in 1:time_steps]
	first_term = [0.0 for i in 1:time_steps]
	times = [0.0 for i in 1:time_steps]
	for i in 1:time_steps
		println(100*i/time_steps," %")
		time_now = i*t_fin/time_steps
		times[i] = time_now
		second_term[i] = v0*(1-exp(-beta*time_now))/beta + x0
		for j in 1:100
			local_time_change = t_fin/time_steps/100
			noise_time = time_now - (101-j)*local_time_change
			noises = white_noise[100*(i-1)+j] - colored_noise[100*(i-1)+j]
			exp_part = 1 - exp(-beta*(time_now-noise_time))
			first_term[i] += noises*exp_part/(beta*lambda)
		end
	end
	full_rez = append!([x0],first_term + second_term)
	times = append!([0.0],times)
	return times,full_rez
end

h = 0.75
final_time = 10
time_steps = 100
beta = 1
lambda = 1
a = 0.5
x0 = 10.0
v0 = 0.01
white_noise = noise(0.5,100*time_steps,final_time,1).*fluc_dissp_coeffs("white",beta,lambda,a,h)
colored_noise = noise(h,100*time_steps,final_time,2).*fluc_dissp_coeffs("color",beta,lambda,a,h)
first_try = get_first_guess(h,time_steps,final_time,beta,lambda,a,white_noise,colored_noise,x0,v0)
plot(first_try[1],first_try[2])




"fin"
