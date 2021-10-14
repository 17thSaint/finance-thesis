using PyPlot
function w_colored_noise(h,n,t_fin,cutoff)
	times = [i*t_fin/n for i in 0:n]
	noise = [0.0 for i in 0:n]
	for i in 0:n
		for j in 0:cutoff
			noise[i+1] += rand(Float64)*rand(-1:2:1)*(0.9^(j*h)) * sin(2*pi*( times[i+1]*(0.9^(-j))+rand(Float64) )) + rand(Float64)*rand(-1:2:1)*(0.9^(-j*h))*sin(2*pi*(times[i+1]*(0.9^j)+rand(Float64)))
		end
	end
	return times,noise
end

for i in 1:10
	omega = 10*i*i
	println("Running Cutoff=",omega)
	here = w_colored_noise(0.5,100,1,Int(omega))
	plot(here[1],here[2],label="$omega")
end
legend()

