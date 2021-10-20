using Cubature, SpecialFunctions, HypergeometricFunctions, PyPlot, MittagLeffler, Statistics


function frac_brown_wiki2(h,n,t_fin)
	times = [i*t_fin/n for i in 0:n]
	dB = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	bh = fill(0.0,n)
	for j in 1:n
		#if j%(0.05*n) == 0
		#	println("H=",h,", ",100*j/n,"%",", Noise")
		#end
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

function calc_b2(h,avging_counts,motion,sigma0)
	#steps = 100
	#final_time = 1
	#motion = [frac_brown_wiki2(h,steps,final_time)[2] for k in 1:avging_counts]
	th_exp_squared = (sigma0^2)*gamma(2-2*h)/(2*h*gamma(1.5-h)*gamma(0.5+h))
	exp_squared = [0.0 for i in 1:avging_counts]
	for i in 1:avging_counts
		exp_squared[i] = (last(motion[i])^2)/avging_counts
	end
	return sum(exp_squared), th_exp_squared
end

function calc_covar(h,avging_counts,motion,sigma0)
	final_time = 1
	steps = 100
	#motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:avging_counts]
	covar = [0.0 for i in 1:steps-1]
	coeff = (sigma0^2)*gamma(2-2*h)/(4*h*gamma(1.5-h)*gamma(0.5+h))
	th_covar = [coeff*(final_time^(2*h)+(i*final_time/steps)^(2*h)-(abs(final_time-i*final_time/steps))^(2*h)) for i in 1:steps-1]
	for j in 1:avging_counts
		for i in 1:steps-1
			covar[i] += motion[j][steps]*motion[j][i]/avging_counts
		end
	end

	return covar,th_covar
end

function calc_sigma(avging_counts,motion)
	h = 0.5
	steps = 100
	final_time = 1
	motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:avging_counts]
	covar = [0.0 for i in 1:steps-1]
	th_covar = [i*final_time/steps for i in 1:steps-1]
	for j in 1:avging_counts
		for i in 1:steps-1
			covar[i] += motion[j][steps]*motion[j][i]/avging_counts
		end
	end
	sigma0 = sqrt.(abs.(covar./th_covar))
	return mean(sigma0), std(sigma0)
end

#=
num = 50
sigs = [0.0 for i in 1:num]
stds = [0.0 for i in 1:num]
for i in 1:num
	println(100*(i-1)/num,"%")
	here = calc_sigma(50)
	sigs[i] = here[1]
	stds[i] = here[2]
end
=#
s0 = 0.5580484388471367#mean(sigs)
#println("Sigma=",s0," +/- ",std(sigs))




#=  Figure of covariance
h = 0.75
steps = 100
final_time = 1
avg_counts = 200
motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:avg_counts]
here = calc_covar(h,avg_counts,motion,s0)
plot(1:steps-1,here[1],1:steps-1,here[2])
=#

#=  Convergence of B^2 for number of counts
h = 0.25
steps = 50
final_time = 1
num_counts = 60
b2_vals = [0.0 for i in 1:num_counts]
th_b2 = [0.0 for i in 1:num_counts]
counts = [0.0 for i in 1:num_counts]
motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:Int(10+(num_counts-1)*10)]
for i in 1:num_counts
	println(i)
	here = Int(10+(i-1)*10)
	counts[i] = here 
	this = calc_b2(h,here,motion,s0)
	b2_vals[i] = this[1]
	th_b2[i] = this[2]
end
plot(counts,b2_vals,counts,th_b2)
=#



counts = 10
avging_counts = 200
h_vals = [0.25 + (i-1)*0.5/counts for i in 1:counts]
b2 = [0.0 for i in 1:counts]
th_b2 = [0.0 for i in 1:counts]
for i in 1:counts
	h = 0.25 + (i-1)*0.5/counts
	println("Motion for i=",i)
	motion_here = [frac_brown_wiki2(h,100,1)[2] for k in 1:avging_counts]
	println("B^2 Calc for i=",i)
	here = calc_b2(h,avging_counts,motion_here,s0)
	b2[i] = here[1]
	th_b2[i] = here[2]
end
plot(h_vals,b2,h_vals,th_b2)
#plot(h_vals,th_b2_vals_b2,label="B^2 TH")
#plot(h_vals,b2_vals_covar,label="Covariance EXP")
#plot(h_vals,th_b2_vals_covar,label="Covariance TH")
	

#= Figure sigma0calc-b2-covar.png shows that the two calculation methods are quite close to each other, but their variances are quite large

h = 0.5
counts = [Int(10+(i-1)*10) for i in 1:20]
s0_b2 = [0.0 for i in 1:20]
s0_covar = [0.0 for i in 1:20]
for i in 1:20
	println("Motion for i=",i)
	motion_here = [frac_brown_wiki2(h,100,1)[2] for k in 1:counts[i]]
	println("Sigma Calculation for i=",i)
	s0_b2[i] = calc_sigma0_b2(h,counts[i],motion_here)
	s0_covar[i] = calc_sigma0_covar(h,counts[i],motion_here)[1]
end
plot(counts,s0_b2,counts,s0_covar)

=#

#=  Figure named expB2-vs-fBMsteps-counts2040 shows that fBM steps doesn't seem to impact calculation much either, but increasing counts dampens variation

h = 0.5
final_time = 1
num_counts = 10
th_exp_squared = [gamma(2-2*h)/(2*h*gamma(1.5-h)*gamma(0.5+h)) for i in 1:num_counts]
exp_squared = [0.0 for i in 1:num_counts]
step_vals = [100+(i-1)*50 for i in 1:num_counts]
for i in 1:num_counts
	steps = 100+(i-1)*50
	motion = [frac_brown_wiki2(h,steps,final_time)[2] for k in 1:40]
	for j in 1:40
		exp_squared[i] += (last(motion[j])^2)/40
	end
end

plot(step_vals,exp_squared,step_vals,th_exp_squared)

=#


#= Figure named expB2-vs-avgingcounts-400steps.png shows that number of averaging steps has no impact on accuracy of calculated value therefore steps=400 is good, take avging counts=20

h = 0.5
final_time = 1
steps = 400
end_step = Int(steps)
num_counts = 12
final_counts = Int(10+(num_counts-1)*5)
motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:final_counts]
th_exp_squared = [gamma(2-2*h)/(2*h*gamma(1.5-h)*gamma(0.5+h)) for i in 1:num_counts]
exp_squared = [0.0 for i in 1:num_counts]
counts = [10+(i-1)*5 for i in 1:num_counts]
for i in 1:num_counts
	avg_counts = 10+(i-1)*5
	for j in 1:Int(avg_counts)
		exp_squared[i] += (motion[j][end_step]^2)/avg_counts
	end
end

plot(counts,exp_squared,counts,th_exp_squared)

=#

#=   Calculate Covariance
h = 0.5
final_time = 1
steps = 200
averaging_counts = 50
motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:averaging_counts]
end_step = Int(steps)
covar = [0.0 for i in 1:end_step-1]
th_covar = [0.5*(1^(2*h)+(i/end_step)^(2*h)-(abs(1-i/end_step))^(2*h)) for i in 1:end_step-1]
for j in 1:averaging_counts
	for i in 1:end_step-1
		println(i/(end_step-1))
		covar[i] += motion[j][end_step]*motion[j][i]/averaging_counts
	end
end
plot(1:end_step-1,th_covar,1:end_step-1,covar)
=#

# 10000:20 bad around 8, 5000:20 bad around 16, 2500:20 good
# t_fin/(t_steps*noise_steps) = 1/200

#=
time_steps = 200
noise_steps = 100
final_time = time_steps/2
for i in 1:3
	h = i/4
	dats = lang_soln(h,time_steps,noise_steps,1,1,final_time)
	plot(dats[1],dats[2],label="$h")
end
legend()
=#



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
