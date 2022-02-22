using SpecialFunctions,PyPlot,CurveFit,MittagLeffler,LaTeXStrings

function get_fractderiv(h,delta_t,steps,og_func,f0,noise_steps,cap=1)
	fract_deriv = [0.0 for i in 1:steps-1]
	times = [i*delta_t for i in 1:steps-1]
	for i in 1:steps-1
		top = noise_steps*i
		if top > 150
			top = 150
		end
		for j in 0:top
			func_here = og_func[noise_steps*i-j+1] - cap*f0
			fact_j = exp(sum([log(i) for i in 1:j]))
			binom_part = gamma(3-2*h)/(fact_j*gamma(3-2*h-j))
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here/(delta_t/noise_steps)^(2-2*h)
		end
	end
	#rl_shift = [cap*f0*(times[i]^(2*h-2))/gamma(2*h-1) for i in 1:steps-1]
	#fract_deriv -= rl_shift
	return fract_deriv
end



time_steps = 100
final_time = 100
dt = final_time/time_steps
times = [i*dt for i in 1:time_steps-1]
alpha = 0.7
h = 0.5*(2 - alpha)
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
pow = 3
f0 = 1.0
#og_func = [((1*t*dt/noise_steps)^pow)*1.0 for t in 1:time_steps*noise_steps]
#og_func = [0.5 for t in 1:time_steps*noise_steps]
og_func = [exp(1.0*t*dt/noise_steps) for t in 1:time_steps*noise_steps]
rezz_Cap = get_fractderiv(h,dt,time_steps,og_func,f0,noise_steps)
#rezz_RL = get_fractderiv(h,dt,time_steps,og_func,f0,noise_steps,0)
theory_exp_Cap = [mittleff(1,2-alpha,times[i]*1.0)*times[i]^(1-alpha) for i in 1:time_steps-1]
#theory_exp_RL = [mittleff(1,1-alpha,times[i]*1.0)*times[i]^(-alpha) for i in 1:time_steps-1]
scatter(times,rezz_Cap,label="Num Cap")
plot(times,theory_exp_Cap,"-g",label="Theory Cap")
#scatter(times,rezz_RL,label="Num RL")
#theory_tpow_Cap = [1.0*factorial(pow)*((1.0*times[i])^(pow-alpha))/gamma(pow+1-alpha) for i in 1:time_steps-1]
#plot(times,theory_exp_RL,"-r",label="Theory RL")

#plot(times,theory_tpow_Cap,label="$alpha")#"Theory Cap")
#plot(hs,hs,label="Theory")
#plot(hs,power_vals)
legend()
#xlabel("Time")
#ylabel("Fractional Derivative of Position")
#title(latexstring("\$ D^{$alpha} [ e^t ] \$"))




"fin"
