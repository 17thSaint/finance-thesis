using HDF5,SpecialFunctions,PyPlot,LaTeXStrings

function frac_deriv_poly(a,time,sigma,alpha)
       val = gamma(sigma+1)*((time-a)^(sigma-alpha))/gamma(sigma-alpha+1)
       return val
end

function read_hdf5_data(h,version)
	cd("fBM-data")
	file = h5open("fBM-$version-h-$h.hdf5","r")
	data = [read(file["values"],"deets_exp"), read(file["values"],"deets_th")]
	cd("..")
	return data
end

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
#=
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
