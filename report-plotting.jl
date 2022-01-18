using SpecialFunctions,PyPlot,LaTeXStrings

function frac_deriv_poly(a,time,sigma,alpha)
       val = gamma(sigma+1)*((time-a)^(sigma-alpha))/gamma(sigma-alpha+1)
       return val
end

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





"fin"
