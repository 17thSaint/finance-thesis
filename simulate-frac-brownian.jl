using HDF5, Cubature, SpecialFunctions, HypergeometricFunctions, PyPlot, MittagLeffler, Statistics, FinancialToolbox


function frac_brown_wiki2(h,n,t_fin)
	times = [i*t_fin/n for i in 0:n]
	dB = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	bh = fill(0.0,n)
	for j in 1:n
		if j%(0.05*n) == 0
			println("H=",h,", ",100*j/n,"%",", Noise")
		end
		for i in 0:j-1
			function integrand(s::Float64)
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

function auto_correlation(energies, delta_t)
    average_energy = mean(energies)
    
    points = Int(floor(length(energies)-delta_t))
    
    energy_fluctuations = energies.-average_energy
    
    autocorrelation_top = [0.0 for i in 1:points]
    autocorrelation_bottom = mean(energy_fluctuations.^2)
    
    for i in 1:points
        autocorrelation_top[i] = energy_fluctuations[i]*energy_fluctuations[i+delta_t]
    end
    
    return (mean(autocorrelation_top)/autocorrelation_bottom)
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
function noise(h,n,t_fin,which=rand(1:20),make_new=false)
	fBM = read_hdf5_data(h,which,true,n+1)[2]
	if make_new
		fBM = frac_brown_wiki2(h,n,t_fin)[2]
	end
	return [t_fin*(fBM[i+1]-fBM[i])/n for i in 1:n]
end

# fluctuation dissipation coefficient
function eta(gam,h,kBT)
	#kBT = 1.0
	return sqrt(2*gam*kBT*gamma(1.5-h)*gamma(0.5+h)/(gamma(2*h)*gamma(2-2*h)))
end

# returns calculated solution to equation 
# m*d^2 x/dt^2 + gamma*D^(alpha) [x(t)] = f_H
# x(t) = 1/m * Convolution(f_H(t), t*E_(2H,2) (-gamma/m * t^(2H))) + v0*t*E_(2H,2) (-gamma/m * t^(2H))
# this assumes that x0 = 0
# uses trapezoid rule for convolution integral
function lang_soln(h,t_steps,noise_steps,gam,m,t_fin,v0,which=rand(1:20))
	times = [i*t_fin/t_steps for i in 0:t_steps]
	term_one = [0.0 for i in 1:t_steps]
	term_two = [0.0 for i in 1:t_steps]
	c_eta = eta(gam,h,1.0)
	noise_term = c_eta.*noise(h,Int(t_steps*noise_steps),t_fin,which)
	for i in 1:t_steps
		if i%(0.05*t_steps) == 0
			println("H=",h,", ",100*i/t_steps,"%",", Lang")
		end
		noise_times = append!([0.0],[ l*t_fin/t_steps + j*t_fin/(t_steps*noise_steps) for l in 0:i-1 for j in 1:noise_steps])
		for k in 1:length(noise_times)-1
			left_time = noise_times[k]
			right_time = noise_times[k+1]
			
			mitlef_left = mittleff(2*h,2,(-gam/m)*(left_time^(2*h)))
			mitlef_right = mittleff(2*h,2,(-gam/m)*(right_time^(2*h)))
				
			if v0 != 0.0 && k == length(noise_times)-1
				term_two[i] += v0*right_time*mitlef_right
			end
			
			if k == length(noise_times)-1
				noise_left = noise_term[1]
			else
				noise_left = noise_term[length(noise_times)-1-k]
			end
			noise_right = noise_term[length(noise_times)-k]
				
			term_one[i] += 0.5*(right_time-left_time)*c_eta*(noise_left*left_time*mitlef_left + noise_right*right_time*mitlef_right )
		end
		
	end
	full_position = append!([0.0],term_one./m + term_two)
	return times,full_position
end





time_steps = 1000
gam = 1.0
mass = 1.0
final_time = 10
v0 = 0.0
h = 0.5
#solns = [lang_soln(i,time_steps,10,gam,mass,final_time,v0,2) for i in [0.3,0.5,0.75]]
#gen_noise = [noise(i,Int(time_steps*100),final_time,1)*eta(gam,i,1.0) for i in [0.3,0.5,0.75]]

soln = lang_soln(h,time_steps,100,gam,mass,final_time,v0,2)


#=
function pair_covar(motion_1,motion_2,count_inter_times)
	len = length(motion_1)
	covar = [[0.0 for i in 1:len] for j in 1:count_inter_times]
	for j in 1:count_inter_times
		interaction_time = Int(j*len/count_inter_times)
		for i in 1:len
			covar[j][i] = motion_1[interaction_time]*motion_2[i]
		end
	end
	return covar
end
len = length(solns[1][2])-1
corrs = [0.0 for i in 1:Int(0.75*len)-1]
dts = [1+(i-1)*1 for i in 1:Int(0.75*len)-1]
for i in 1:Int(0.75*len)-1
	println(i/(Int(0.75*len)-1))
	#corrs[i] = cor(solns[1][2][1+(i-1)*10:i*10],solns[1][2][92:101])
	#corrs[i] = cor(gen_noise[1][1+(i-1)*1000:i*1000],gen_noise[1][9001:10000])
	corrs[i] = auto_correlation(solns[1][2],dts[i])
end
plot(dts,corrs)
=#
#=	This is the if statement loop for associativity of convolution
ex_mit = lang_soln(h,time_steps,100,a,lambda,final_time,v0,"mitag",2)
ex_noise = lang_soln(h,time_steps,100,a,lambda,final_time,v0,"blah",2)
			if shifted == "mitag"
				mitlef_left = mittleff(2*h,2,(-gam/m)*(left_time^(2*h)))
				mitlef_right = mittleff(2*h,2,(-gam/m)*(right_time^(2*h)))
				
				if v0 != 0.0 && k == length(noise_times)-1
					term_two[i] += v0*right_time*mitlef_right
				end
				
				if k == length(noise_times)-1
					noise_left = noise_term[1]
				else
					noise_left = noise_term[length(noise_times)-1-k]
				end
				noise_right = noise_term[length(noise_times)-k]
				term_one[i] += 0.5*(right_time-left_time)*c_eta*(noise_left*left_time*mitlef_left + noise_right*right_time*mitlef_right )
				#println("Mitag: ",left_time)
			elseif shifted == "blah"
				mitlef_right = mittleff(2*h,2,(-gam/m)*((times[i+1]-left_time)^(2*h)))
				mitlef_left = mittleff(2*h,2,(-gam/m)*((times[i+1]-right_time)^(2*h)))
				#println("Noise: ",mitlef_left,", ",mitlef_right)
				if v0 != 0.0 && k == length(noise_times)-1
					term_two[i] += v0*right_time*mitlef_right
				end
		
				if k == 1
					noise_left = noise_term[1]
				else
					noise_left = noise_term[k-1]
				end
				noise_right = noise_term[k]
				term_one[i] += 0.5*(right_time-left_time)*c_eta*(noise_left*(times[i+1]-right_time)*mitlef_left + noise_right*(times[i+1]-left_time)*mitlef_right )
				#println("Noise: ",left_time)
			else
				return "Error"
			end
			=#



#=  Look at volatility of price, approximation used is weird and likely wrong
lam = 100
beta = 100
m = lam
volat = [0.0 for i in 1:10]
as = [0.0 for i in 1:10]
h = 0.5
for i in 1:10
	a = (i-1)*(10*1)/10
	println("A=$a")
	as[i] = a
	gam = lam*beta - a
	local_volat = 0
	for k in 1:20
		res = lang_soln(h,100,100,gam,m,1,0,k)
		local_volat += (std(res[4])^2)/20
	end
	volat[i] = local_volat
end
d = [volat[i]*0.5*(beta-as[i]/lam)^2 for i in 1:10]
d_th = [(lam*beta-as[i])/(lam^2) for i in 1:10]
plot(as,d,as,d_th)

# larger mass leads to less volatility but not linear relationship
=#

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

function calc_covar(h,avging_counts,motion,interaction_time=10000,sigma0=0.55804)
	final_time = 10000
	now_time = final_time
	steps = 10000
	#motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:avging_counts]
	covar = [0.0 for i in 1:steps-1]
	#coeff = (sigma0^2)*gamma(2-2*h)/(4*h*gamma(1.5-h)*gamma(0.5+h))
	#th_covar = [coeff*(now_time^(2*h)+(i*final_time/steps)^(2*h)-(abs(now_time-i*final_time/steps))^(2*h)) for i in 1:steps-1]
	for j in 1:avging_counts
		for i in 1:steps-1
			covar[i] += motion[j][interaction_time]*motion[j][i]/avging_counts
		end
	end

	return covar#,th_covar
end
#=
# looking at covariance of noise which for h=0.5 should be dirac delta, trying to find max
# of set to find on average what lambda^2 D should be
h = 0.5
motion_here = [noise(h,9999,1,i) for i in 1:20]
#motion_here = [read_hdf5_data(h,i,true,100)[2] for i in 1:20]
maxes = [0.0 for i in 1:5000]
cross_vals = [0.0 for i in 1:5000]
for i in 1:5000
	println(i)
	crossing = Int(1 + (i-1)*2)
	cross_vals[i] = crossing
	covars = calc_covar(h,20,motion_here,0.55,crossing)
	maxes[i] = maximum(covars[1])
end
plot(cross_vals,maxes)
=#
function calc_sigma(avging_counts,motion)
	h = 0.5
	steps = 200
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
	here = calc_sigma(50,1)
	sigs[i] = here[1]
	stds[i] = here[2]
end
plot(sigs)
=#
#s0 = 0.5580484388471367#mean(sigs)
#println("Sigma=",s0," +/- ",std(sigs))

	


#=  Figure of covariance
h = 0.75
steps = 100
final_time = 1
avg_counts = 20
#motion = [frac_brown_wiki2(h,steps,final_time)[2] for i in 1:avg_counts]
motion = [read_hdf5_data(h,i,true,steps)[2] for i in 1:20]
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


#=
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
=#	

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

#=
function asset_price(n,t_fin,volat,corr,r,start_price)
	dB = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	dW = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	asset = append!([start_price],[0.0 for i in 1:n])
	times = [i*t_fin/n for i in 0:n]
	for i in 2:n+1
		if i%(0.05*n) == 0
			println("H=",h,", ",100*(i-1)/n,"%",", Asset Price")
		end
		change_asset = r*asset[i-1]*(t_fin/n)+volat[i-1]*asset[i-1]*(corr*dW[i-1]+sqrt(1-corr^2)*dB[i-1])
		asset[i] = asset[i-1] + change_asset
	end
	return times,asset
end

function get_imp_vol(asset_now,r,t_mat,volat,points)
	range_strike = [0.7*asset_now + i*0.6*asset_now/points for i in 0:points-1]
	implied_vol = [0.0 for i in 1:points]
	for i in 1:points
		option_price = blkprice(asset_now,range_strike[i],r,t_mat,volat)
		implied_vol[i] = blsimpv(asset_now,range_strike[i],r,t_mat,option_price)
	end
	return range_strike,implied_vol
end

function get_volatility(h,t_steps,gam,m,t_fin,volat0)
	v0 = log(volat0)
	volat = append!([volat0],[0.0 for i in 1:t_steps-1])
	base = lang_soln(h,t_steps,100,gam,m,t_fin,v0)
	for i in 1:t_steps-1
		if i%(0.05*t_steps) == 0
			println("H=",h,", ",100*i/t_steps,"%",", Volatility")
		end
		volat[i+1] = exp((base[2][i+1]-base[2][i])/(t_fin/t_steps))
	end
	return base,volat
end
=#



"fin"
