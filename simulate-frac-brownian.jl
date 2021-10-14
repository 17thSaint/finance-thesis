using Cubature, SpecialFunctions, HypergeometricFunctions, PyPlot, MittagLeffler


function frac_brown_wiki2(h,n,t_fin)
	times = [i*t_fin/n for i in 0:n]
	dB = append!([0.0],[rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n])
	bh = fill(0.0,n+1)

	for j in 0:n
		for i in 0:j-1
			function integrand(s)
				return (((times[j+1]-s)^(h-0.5))/gamma(h+0.5))*_₂F₁(h-0.5,0.5-h,h+0.5,1-times[j+1]/s)
			end
			part = hquadrature(integrand,times[i+1],times[i+2])[1]*dB[i+1]*n/t_fin
			bh[j+1] += part
			#println(part,"j=",j,", ","i=",i)
		end
	end
	return times, bh
end

function noise(h,n,t_fin)
	fBM = frac_brown_wiki2(h,n,t_fin)[2]
	return [t_fin*(fBM[i+1]-fBM[i])/n for i in 1:n]
end

function eta(gam,h,kBT)
	#kBT = 1.0
	return sqrt(2*gam*kBT*gamma(1.5-h)*gamma(0.5+h)/(gamma(2*h)*gamma(2-2*h)))
end

function position(h,n,t_fin,gam,m,v0)
	times = [i*t_fin/n for i in 0:n]
	here = [0.0 for i in 0:n]
	for i in 1:n
		println(i)
		for j in 0:i-1
			mitlef = mittleff(2*h,2,(-gam/m)*(times[i+1]^(2*h)))
			here[i+1] += eta(gam,h,1.0)*noise(h,n,t_fin)[n-j]*times[i+1]*mitlef
		end
	end
	
	sec_term = v0.*times.*[mittleff(2*h,2,(-gam/m)*(times[i]^(2*h))) for i in 1:n+1]
	
	return times, here./m + sec_term
end

for i in 1:3
	test1 = 0
	test1 = position(i/4,50,1,1,1,0)
	plot(test1[1],test1[2],label="$i")
end
legend()





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
