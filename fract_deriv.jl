using SpecialFunctions,PyPlot,CurveFit

function get_fractderiv(h,delta_t,steps,og_func,lim_val=0.1)
	fract_deriv = [0.0 for i in 1:steps-1]
	times = [0.0 for i in 1:steps-1]
	for i in 1:steps-1
		times[i] = i*delta_t
		for j in 0:i
			func_here = og_func[i-j+1]
			binom_part = gamma(3-2*h)/(factorial(j)*gamma(3-2*h-j))
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here
		end
	end
	return fract_deriv,times
end


dt = 1.0
time_steps = 20

hs = [0.2,0.3,0.4,0.45,0.55,0.6,0.7,0.8]
power_vals = [0.0 for j in 1:length(hs)]
pow = 8
for j in 1:length(hs)
	h = hs[j]
	og_func = [(t*dt)^pow for t in 1:time_steps]
	rezz = get_fractderiv(h,dt,time_steps,og_func)
	fit_vals = power_fit(rezz[2],rezz[1])
	power_vals[j] = 0.5*(fit_vals[2]+2-pow)
end
plot(hs,hs,label="Theory")
plot(hs,power_vals)
legend()

