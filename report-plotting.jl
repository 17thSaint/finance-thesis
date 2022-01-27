using HDF5,SpecialFunctions,PyPlot,LaTeXStrings,MittagLeffler

function frac_deriv_poly(a,time,sigma,alpha)
       val = gamma(sigma+1)*((time-a)^(sigma-alpha))/gamma(sigma-alpha+1)
       return val
end

function read_hdf5_data(h,version)
	#cd("fBM-data")
	#file = h5open("fBM-$version-h-$h.hdf5","r")
	#data = [read(file["values"],"deets_exp"), read(file["values"],"deets_th")]
	file = h5open("fraclang-soln-h-$h-$version.hdf5","r")
	data = [read(file["values"],"soln"),read(file["values"],"config_time"),read(file["values"],"resids_time"),read(file["values"],"soln_step")]
	#cd("..")
	return data
end

function get_ml(h2,xs)
	count = length(xs)
	ml = [mittleff(h2,2,xs[j]) for j in 1:count]
	return ml,h2
end

function get_bincoef(alph,n)
	top = gamma(alph+1)
	bot = gamma(n+1)*gamma(alph-n+1)
	return top/bot
end

ns = [i for i in 2:7]
alphs = [0.1 + i*.89/50 for i in 0:50]
bin = [[get_bincoef(alphs[i],ns[j]) for i in 1:51] for j in 1:6]
for i in 1:6
	lab = ns[i]
	plot(alphs,bin[i],label="$lab")
end
legend()
xlabel(latexstring("\$ \\alpha \$"))
title(latexstring("Binomial Coefficient(\$ \\alpha, n \$)"))




#=  Mittag-leffler function
final_x = -8.1
start_x = -1.1
count = 50
xs = [start_x + i*(final_x-start_x)/count for i in 0:count]
hs = [1.01 + i*0.98/5.0 for i in 0:5]
for i in 1:6
	mls = get_ml(hs[i],xs)
	lab = round(mls[2],digits=2)
	plot(xs,mls[1],label=latexstring("\$ E_{$lab,2} \$"))
end
legend()
xlabel("x")
title("Mittag-Leffler Function")
=#
#= plotting covar
h = 0.73
xs = [i/500 for i in 1:499]
dats = read_hdf5_data(h,"covar")
plot(xs,dats[1],label="Sim")
plot(xs,dats[2],label="Th")
xlabel("Time (t)")
ylabel(latexstring("\$ \\langle B(t) B(t=1.0) \\rangle  \$"))
title("Covariance of fBM, H=$h")
legend()
=#

#=  b2 plotting fBM
h_start = 0.25
h_end = 0.75
dh = 0.01
h_count = Int((h_end-h_start)/dh)
hs = [h_start + dh*i for i in 0:h_count]
b2s = fill(0.0,(h_count+1,2))
for i in 1:h_count+1
	b2s[i,:] = read_hdf5_data(hs[i],"covar")
end
labels = ["Sim","Th"]
for i in 1:2
	lab = labels[i]
	if i == 1
		plot(hs,b2s[:,i],"-p",label="$lab")
	else
		plot(hs,b2s[:,i],label="$lab")
	end
end
legend()
xlabel("Hurst parameter, H")
ylabel(latexstring("\$ \\mathbb{E}[B_{H}^{2}(t=1)]  \$"))
=#

#=
a = 0
pow = 2
time_start = 0
dt = 0.001
times_count = 4000
times = [time_start+i*dt for i in 0:times_count-1]

for alp in [round(-0.2*i,digits=2) for i in 1:4]
	vals = [frac_deriv_poly(a,i,pow,alp) for i in times]
	plot(times,vals,label="$alp")
end
plot(times,[frac_deriv_poly(a,i,pow,-1) for i in times],"-k",label="-1.0")
legend()
xlabel("Time")
ylabel(latexstring("\$ D^{\\alpha} t^2 \$"))
=#




"fin"
