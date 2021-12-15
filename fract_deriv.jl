using SpecialFunctions,PyPlot,CurveFit

function get_fractderiv(h,delta_t,steps,og_func,noise_steps,lim_val=0.1)
	fract_deriv = [0.0 for i in 1:steps-1]
	#times = [0.0 for i in 1:steps-1]
	for i in 1:steps-1
		#times[i] = i*delta_t
		noise_times = append!([0.0],[ l*delta_t + j*delta_t/noise_steps for l in 0:i-1 for j in 1:noise_steps])
		for j in 1:length(noise_times)-1
			func_here = og_func[j]
			binom_part = gamma(3-2*h)/(factorial(j)*gamma(3-2*h-j))
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here
		end
	end
	return fract_deriv#,times
end


dt = 1.0
time_steps = 5
h = 0.3
alpha = 2-2*h
times = [i*dt for i in 1:time_steps-1]

noises = [1,2,3]
#power_vals = [0.0 for j in 1:length(hs)]
pow = 3
for j in 1:length(noises)
	noise_steps = noises[j]
	og_func = [(t*dt)^pow for t in 1:time_steps*noise_steps]
	rezz = get_fractderiv(h,dt,time_steps,og_func,noise_steps)
	plot(times,rezz,label="Steps=$noise_steps")
	#fit_vals = power_fit(rezz[2],rezz[1])
	#power_vals[j] = 0.5*(fit_vals[2]+2-pow)
end
theory = [(factorial(pow)*times[i]^(pow-alpha))/gamma(pow+1-alpha) for i in 1:time_steps-1]
plot(times,theory,label="Theory")
#plot(hs,hs,label="Theory")
#plot(hs,power_vals)
legend()

