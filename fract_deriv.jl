#using SpecialFunctions,PyPlot,CurveFit,MittagLeffler

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
alpha = 0.1
h = 0.5*(2 - alpha)
times = [i*dt for i in 1:time_steps-1]

noise_steps = 5
#power_vals = [0.0 for j in 1:length(hs)]
#pow = 1
#=
for j in 1:length(noises)
	noise_steps = noises[j]
	#noise_times = [i*dt/noise_steps for i in 1:time_steps*noise_steps]
	og_func = [(t*dt/noise_steps)^pow for t in 1:time_steps*noise_steps]
	rezz = get_fractderiv(h,dt,time_steps,og_func,noise_steps)
	plot(times,rezz)
	#plot(noise_times,og_func,label="Steps=$noise_steps")
end
=#
#og_func = [(t*dt/noise_steps)^pow for t in 1:time_steps*noise_steps]
og_func = [sum([((t*dt/noise_steps)^i)/factorial(i) for i in 1:10]) for t in 1:time_steps*noise_steps]
#og_func = [exp(t*dt/noise_steps) for t in 1:time_steps*noise_steps]
rezz = get_fractderiv(h,dt,time_steps,og_func,noise_steps)
scatter(times,rezz,label="Calculated")
#theory = [sum([(factorial(j)*times[i]^(j-alpha))/gamma(j+1-alpha)/factorial(j) for j in 1:10]) for i in 1:time_steps-1]
theory = [mittleff(1,2-alpha,times[i])*times[i]^(1-alpha) for i in 1:time_steps-1]
plot(times,theory,"r",label="Theory")
#plot(hs,hs,label="Theory")
#plot(hs,power_vals)
legend()


