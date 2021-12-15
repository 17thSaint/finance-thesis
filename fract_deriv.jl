#using SpecialFunctions,PyPlot,CurveFit

function get_fractderiv(h,delta_t,steps,og_func,noise_steps,lim_val=0.1)
	fract_deriv = [0.0 for i in 1:steps-1]
	for i in 1:steps-1
		for j in 0:noise_steps*i
			func_here = og_func[noise_steps*i-j+1]
			fact_j = exp(sum([log(i) for i in 1:j]))
			binom_part = gamma(3-2*h)/(fact_j*gamma(3-2*h-j))
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here/(delta_t/noise_steps)^(2-2*h)
		end
	end
	return fract_deriv
end


time_steps = 10
dt = 1/time_steps
h = 0.3
alpha = 2-2*h
times = [i*dt for i in 1:time_steps-1]

noises = [10]
#power_vals = [0.0 for j in 1:length(hs)]
pow = 2
for j in 1:length(noises)
	noise_steps = noises[j]
	#noise_times = [i*dt/noise_steps for i in 1:time_steps*noise_steps]
	og_func = [(t*dt/noise_steps)^pow for t in 1:time_steps*noise_steps]
	rezz = get_fractderiv(h,dt,time_steps,og_func,noise_steps)
	plot(times,rezz,label="Steps=$noise_steps")
	#plot(noise_times,og_func,label="Steps=$noise_steps")
end
theory = [(factorial(pow)*times[i]^(pow-alpha))/gamma(pow+1-alpha) for i in 1:time_steps-1]
plot(times,theory,label="Theory")
#plot(hs,hs,label="Theory")
#plot(hs,power_vals)
legend()

