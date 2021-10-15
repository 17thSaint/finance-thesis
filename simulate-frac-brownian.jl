using Cubature, SpecialFunctions, HypergeometricFunctions, PyPlot, MittagLeffler


function frac_brown_wiki2(h,n,t_fin)
	times = [i*t_fin/n for i in 0:n]
	dB = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	bh = fill(0.0,n)
	for j in 1:n
		if j%(0.05*n) == 0
			println("H=",h,", ",100*j/n,"%",", Noise")
		end
		for i in 0:j-1
			function integrand(s)
				return (((times[j+1]-s)^(h-0.5))/gamma(h+0.5))*_₂F₁(h-0.5,0.5-h,h+0.5,1-times[j+1]/s)
			end
			part = hquadrature(integrand,times[i+1],times[i+2])[1]*dB[i+1]*n/t_fin
			bh[j] += part
			#println(part,"j=",j,", ","i=",i)
		end
	end
	full_bh = append!([0.0],bh)
	return times, full_bh
end

function noise(h,n,t_fin)
	fBM = frac_brown_wiki2(h,n,t_fin)[2]
	return [t_fin*(fBM[i+1]-fBM[i])/n for i in 1:n]
end

function eta(gam,h,kBT)
	#kBT = 1.0
	return sqrt(2*gam*kBT*gamma(1.5-h)*gamma(0.5+h)/(gamma(2*h)*gamma(2-2*h)))
end

function lang_soln(h,t_steps,noise_steps,gam,m,t_fin)
	times = [i*t_fin/t_steps for i in 0:t_steps]
	term_one = [0.0 for i in 1:t_steps]
	#term_two = [0.0 for i in 1:t_steps]
	noise_term = noise(h,Int(t_steps*noise_steps),t_fin)
	c_eta = eta(gam,h,1)
	for i in 1:t_steps
		if i%(0.05*t_steps) == 0
			println("H=",h,", ",100*i/t_steps,"%",", Lang")
		end
		noise_times = [ times[i] + j*t_fin/(t_steps*noise_steps) for j in 0:noise_steps ]
		if i != 1
			term_one[i] += term_one[i-1]
		end
		for k in 0:noise_steps-1
			time_now = noise_times[k+1]
			mitlef = mittleff(2*h,2,(-gam/m)*(time_now^(2*h)))
			#if k == 0
			#	term_two[i] += v0*time_now*mitlef
			#end
			term_one[i] += c_eta * noise_term[noise_steps*i-k] * time_now * mitlef
		end
	end
	full_position = append!([0.0],term_one./m) #+ term_two)
	return times,full_position
end

# 10000:20 bad around 8, 5000:20 bad around 16, 2500:20 good

function return_lang_soln(h,t_steps,noise_steps,gam,m,t_fin)
	if t_fin/(t_steps*noise_steps) < 1.0/125.0
		proper_noise = 125.0*t_fin/t_steps
		println("STOP Bad Parameters: Make Noise Steps=$proper_noise")
		return 0,0
	else
		return lang_soln(h,t_steps,noise_steps,gam,m,t_fin)
	end
end

la0p25 = return_lang_soln(0.25,100,50,1,1,40)
plot(la0p25[1],la0p25[2],label="0.25")
#la0p5 = lang_soln(0.5,100,50,1,1,40)
#la0p75 = lang_soln(0.75,100,50,1,1,40)






#=
n_here = 500
t_fin_here = 10
for i in 1:3
	h_here = i/4
	plot(frac_brown_wiki2(h_here,n_here,t_fin_here)[1],frac_brown(h_here,n_here,t_fin_here)[2],label="$h_here")
end
legend()
=#



"fin"
