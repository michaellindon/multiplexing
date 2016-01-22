srand(2)
ł=0.5
p=1
ν=p+0.5
λ=sqrt(2*ν)/ł
ρ²=10.0
q=2*ρ²*√π*λ^(2*p+1)*gamma(p+1)/gamma(p+0.5)
d=p+1
Tobs=1
g=realization(λ,q,collect(0:Tobs/1000:Tobs))
μg=0.0

#Create 0 branch
f₀=realization(λ,q,collect(0:Tobs/1000:Tobs))
Λ=4000
Tₚ=rand(PPProcess(Λ),0,Tobs)
gₚ=FFBS2(g,Tₚ,λ,q)
f₀ₚ=FFBS2(f₀,Tₚ,λ,q)
μf₀=0.0
ξ₀ₐₐ=Dict()
ξ₀ₐᵣ=Dict()
ξ₀ᵣ=Dict()
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
		ξ₀ᵣ[t]=Dict("yf"=>yf)
	end
end

figure()
subplot(311)
plot(sort(collect(keys(f₀ₚ))),[f₀ₚ[key][1] for key in sort(collect(keys(f₀ₚ)))])
#=figure()=#
#=subplot(311)=#
#=plot(sort(collect(keys(f₀ₚ))),[Φ(gₚ[key][1])*Λ*Φ(f₀ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₀ₐₐ)),-0*ones(length(ξ₀ₐₐ)),c="green",marker="|",linestyle="None",markersize=20)=#
#=subplot(312)=#
#=plot(sort(collect(keys(f₀ₚ))),[(1-Φ(gₚ[key][1]))*Λ*Φ(f₀ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₀ₐᵣ)),-0*ones(length(ξ₀ₐᵣ)),c="green",marker="|",linestyle="None",markersize=20)=#
#=subplot(313)=#
#=plot(sort(collect(keys(f₀ₚ))),[Λ-Λ*Φ(f₀ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₀ᵣ)),-0*ones(length(ξ₀ᵣ)),c="green",marker="|",linestyle="None",markersize=20)=#



#Create 1 branch
f₁=realization(λ,q,collect(0:Tobs/1000:Tobs))
Λ=1000
Tₚ=rand(PPProcess(Λ),0,Tobs)
gₚ=FFBS2(g,Tₚ,λ,q)
f₁ₚ=FFBS2(f₁,Tₚ,λ,q)
μf₁=0.0
ξ₁ₐₐ=Dict()
ξ₁ₐᵣ=Dict()
ξ₁ᵣ=Dict()
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
		ξ₁ᵣ[t]=Dict("yf"=>yf)
	end
end

subplot(312)
plot(sort(collect(keys(f₁ₚ))),[f₁ₚ[key][1] for key in sort(collect(keys(f₁ₚ)))])
subplot(313)
plot(sort(collect(keys(gₚ))),[gₚ[key][1] for key in sort(collect(keys(gₚ)))])

#=figure()=#
#=subplot(311)=#
#=plot(sort(collect(keys(f₁ₚ))),[Φ(gₚ[key][1])*Λ*Φ(f₁ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₁ₐᵣ)),-0*ones(length(ξ₁ₐᵣ)),c="green",marker="|",linestyle="None",markersize=20)=#
#=subplot(312)=#
#=plot(sort(collect(keys(f₁ₚ))),[(1-Φ(gₚ[key][1]))*Λ*Φ(f₁ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₁ₐₐ)),-0*ones(length(ξ₁ₐₐ)),c="green",marker="|",linestyle="None",markersize=20)=#
#=subplot(313)=#
#=plot(sort(collect(keys(f₁ₚ))),[Λ-Λ*Φ(f₁ₚ[key][1]) for key in sort(collect(keys(gₚ)))])=#
#=plot(collect(keys(ξ₁ᵣ)),-0*ones(length(ξ₁ᵣ)),c="green",marker="|",linestyle="None",markersize=20)=#


ξ=Dict()
for t in keys(ξ₀ₐₐ)
	ξ[t]=Dict("yf"=>ξ₀ₐₐ[t]["yf"],"yg"=>ξ₀ₐₐ[t]["yg"],"γ"=>0)
end
for t in keys(ξ₁ₐₐ)
	ξ[t]=Dict("yf"=>ξ₁ₐₐ[t]["yf"],"yg"=>ξ₁ₐₐ[t]["yg"],"γ"=>1)
end

T=union(collect(keys(ξ₁ₐₐ)),collect(keys(ξ₀ₐₐ)))
