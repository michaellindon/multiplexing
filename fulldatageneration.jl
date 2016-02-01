srand(2)
ł=0.1
p=2
ν=p+0.5
ρ²=3.0
d=p+1
Tobs=1
g=realization(ł,ρ²,collect(0:Tobs/1000:Tobs))
μg=0.0

#Create 0 branch
f₀=realization(ł,ρ²,collect(0:Tobs/1000:Tobs))
Λ=2000
Tₚ=rand(PPProcess(Λ),0,Tobs)
gₚ=FFBS2(g,Tₚ,ł,ρ²)
f₀ₚ=FFBS2(f₀,Tₚ,ł,ρ²)
μf₀=0.0
ξ₀ₐₐ=Dict()
ξ₀ₐᵣ=Dict()
ξ₀ᵣᵣ=Dict()
for t in Tₚ
	yf=rand(Normal(μf₀+f₀ₚ[t][1],1))
	if(yf>0)
		yg=rand(Normal(μg+gₚ[t][1],1))
		if(yg>0)
			ξ₀ₐₐ[t]=Dict("yf"=>yf,"yg"=>yg)
		else
			ξ₀ₐᵣ[t]=Dict("yf"=>yf,"yg"=>yg)
		end
	else
		ξ₀ᵣᵣ[t]=Dict("yf"=>yf)
	end
end

figure()
subplot(411)
plot(convert(Array{Float64,1},sort(collect(keys(f₀)))),convert(Array{Float64},[Λ*Φ(μf₀+f₀[key][1]) for key in sort(collect(keys(f₀)))]),c="blue")

#=plot(sort(collect(keys(f₀ₚ))),[f₀ₚ[key][1] for key in sort(collect(keys(f₀ₚ)))])=#
#=figure()=#
#=subplot(311)=#
#=plot(sort(collect(keys(f₀ₚ))),[Φ(gₚ[key][1])*Λ*Φ(f₀ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₀ₐₐ)),-0*ones(length(ξ₀ₐₐ)),c="green",marker="|",linestyle="None",markersize=20)=#
#=subplot(312)=#
#=plot(sort(collect(keys(f₀ₚ))),[(1-Φ(gₚ[key][1]))*Λ*Φ(f₀ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₀ₐᵣ)),-0*ones(length(ξ₀ₐᵣ)),c="green",marker="|",linestyle="None",markersize=20)=#
#=subplot(313)=#
#=plot(sort(collect(keys(f₀ₚ))),[Λ-Λ*Φ(f₀ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₀ᵣᵣ)),-0*ones(length(ξ₀ᵣᵣ)),c="green",marker="|",linestyle="None",markersize=20)=#
Tₚ=rand(PPProcess(Λ),0,Tobs)
f₀ₚ=FFBS2(f₀,Tₚ,ł,ρ²)
ξ₀ₐ=Dict()
ξ₀ᵣ=Dict()
for t in Tₚ
	yf=rand(Normal(μf₀+f₀ₚ[t][1],1))
	if(yf>0)
		ξ₀ₐ[t]=Dict("yf"=>yf)
	else
		ξ₀ᵣ[t]=Dict("yf"=>yf)
	end
end
plot(convert(Array{Float64,1},sort(collect(keys(ξ₀ₐ)))),-0*ones(length(ξ₀ₐ)),c="blue",marker="|",linestyle="None",markersize=20)



#Create 1 branch
f₁=realization(ł,ρ²,collect(0:Tobs/1000:Tobs))
Tₚ=rand(PPProcess(Λ),0,Tobs)
gₚ=FFBS2(g,Tₚ,ł,ρ²)
f₁ₚ=FFBS2(f₁,Tₚ,ł,ρ²)
μf₁=0.0
ξ₁ₐₐ=Dict()
ξ₁ₐᵣ=Dict()
ξ₁ᵣᵣ=Dict()
for t in Tₚ
	yf=rand(Normal(μf₁+f₁ₚ[t][1],1))
	if(yf>0)
		yg=rand(Normal(μg+gₚ[t][1],1))
		if(yg>0)
			ξ₁ₐᵣ[t]=Dict("yf"=>yf,"yg"=>yg)
		else
			ξ₁ₐₐ[t]=Dict("yf"=>yf,"yg"=>yg)
		end
	else
		ξ₁ᵣᵣ[t]=Dict("yf"=>yf)
	end
end

subplot(412)
plot(convert(Array{Float64,1},sort(collect(keys(f₁)))),convert(Array{Float64},[Λ*Φ(μf₁+f₁[key][1]) for key in sort(collect(keys(f₁)))]),c="red")
#=plot(sort(collect(keys(f₁ₚ))),[f₁ₚ[key][1] for key in sort(collect(keys(f₁ₚ)))])=#
#=plot(sort(collect(keys(gₚ))),[gₚ[key][1] for key in sort(collect(keys(gₚ)))])=#
Tₚ=rand(PPProcess(Λ),0,Tobs)
f₁ₚ=FFBS2(f₁,Tₚ,ł,ρ²)
ξ₁ₐ=Dict()
ξ₁ᵣ=Dict()
for t in Tₚ
	yf=rand(Normal(μf₁+f₁ₚ[t][1],1))
	if(yf>0)
		ξ₁ₐ[t]=Dict("yf"=>yf)
	else
		ξ₁ᵣ[t]=Dict("yf"=>yf)
	end
end
plot(convert(Array{Float64,1},sort(collect(keys(ξ₁ₐ)))),-0*ones(length(ξ₁ₐ)),c="red",marker="|",linestyle="None",markersize=20)

subplot(413)
plot(convert(Array{Float64,1},sort(collect(keys(g)))),convert(Array{Float64},[Φ(μg+g[key][1]) for key in sort(collect(keys(g)))]),c="green")
subplot(414)
plot(convert(Array{Float64,1},sort(collect(keys(f₀)))),convert(Array{Float64},[Λ*Φ(μf₀+f₀[key][1]) for key in sort(collect(keys(f₀)))]),c="blue",linestyle="--")
plot(convert(Array{Float64,1},sort(collect(keys(f₁)))),convert(Array{Float64},[Λ*Φ(μf₁+f₁[key][1]) for key in sort(collect(keys(f₁)))]),c="red",linestyle="--")
fp₁=FFBS2(f₁,collect(0:Tobs/1000:Tobs),ł,ρ²)
fp₀=FFBS2(f₀,collect(0:Tobs/1000:Tobs),ł,ρ²)
gp=FFBS2(g,collect(0:Tobs/1000:Tobs),ł,ρ²)
plot(convert(Array{Float64,1},sort(collect(keys(gp)))),convert(Array{Float64},[Φ(μg+g[key][1])*Λ*Φ(μf₀+fp₀[key][1])+(1-Φ(μg+g[key][1]))*Λ*Φ(μf₁+fp₁[key][1]) for key in sort(collect(keys(g)))]),c="purple")
plot(sort(collect(keys(ξ₀ₐₐ))),0*ones(length(ξ₀ₐₐ)),c="blue",marker="|",linestyle="None",markersize=20)
plot(sort(collect(keys(ξ₁ₐₐ))),0*ones(length(ξ₁ₐₐ)),c="red",marker="|",linestyle="None",markersize=20)
#=figure()=#
#=subplot(311)=#
#=plot(sort(collect(keys(f₁ₚ))),[Φ(gₚ[key][1])*Λ*Φ(f₁ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₁ₐᵣ)),-0*ones(length(ξ₁ₐᵣ)),c="green",marker="|",linestyle="None",markersize=20)=#
#=subplot(312)=#
#=plot(sort(collect(keys(f₁ₚ))),[(1-Φ(gₚ[key][1]))*Λ*Φ(f₁ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₁ₐₐ)),-0*ones(length(ξ₁ₐₐ)),c="green",marker="|",linestyle="None",markersize=20)=#
#=subplot(313)=#
#=plot(sort(collect(keys(f₁ₚ))),[Λ-Λ*Φ(f₁ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₁ᵣᵣ)),-0*ones(length(ξ₁ᵣᵣ)),c="green",marker="|",linestyle="None",markersize=20)=#


ξ=Dict()
for t in keys(ξ₀ₐₐ)
	ξ[t]=Dict("yf"=>ξ₀ₐₐ[t]["yf"],"yg"=>ξ₀ₐₐ[t]["yg"],"γ"=>0)
end
for t in keys(ξ₁ₐₐ)
	ξ[t]=Dict("yf"=>ξ₁ₐₐ[t]["yf"],"yg"=>ξ₁ₐₐ[t]["yg"],"γ"=>1)
end

Ts=union(collect(keys(ξ₁ₐₐ)),collect(keys(ξ₀ₐₐ)),collect(keys(ξ₀ₐ)),collect(keys(ξ₁ₐ)))
Td=union(collect(keys(ξ₁ₐₐ)),collect(keys(ξ₀ₐₐ)))
f₁=FFBS2(f₁,Ts,ł,ρ²)
f₀=FFBS2(f₀,Ts,ł,ρ²)
g=FFBS2(g,Td,ł,ρ²)
