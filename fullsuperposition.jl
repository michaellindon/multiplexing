importall Base.Random
importall Distributions
using PyPlot

include("poissonpointprocess.jl")
include("analyticalEigen.jl")
include("statespace.jl")
include("fulldatageneration.jl")


niter=40
iter=1
gᵧ=1
σ²ₘ=3
accepted=0
r=2
gᵧTrace=zeros(niter)
ρ²Trace=zeros(niter)
łTrace=zeros(niter)
μTrace=zeros(niter)
logodds=1
logoddsTrace=zeros(niter)
σ²=1.0
μ=0.0
for key in keys(f₀)
	f₀[key]=zeros(d,1)
end
for key in keys(f₁)
	f₁[key]=zeros(d,1)
end
for key in keys(g)
	g[key]=zeros(d,1)
end

for iter=1:niter
	println(iter)
	Ξₚ=PPProcess(Λ)

	#Sample ξ₀ₐ
	for t in keys(ξ₀ₐ)
			ξ₀ₐ[t]["yf"]=rand(Truncated(Normal(f₀[t][1],1),0,Inf))
	end

	#Sample ξ₀ᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₀,Tₚ,ł,ρ²)
	ξ₀ᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(fₚ[t][1])))==1)
			ξ₀ᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(fₚ[t][1],1),-Inf,0)))
		end
	end

	#Sample ξ₁ₐ
	for t in keys(ξ₁ₐ)
			ξ₁ₐ[t]["yf"]=rand(Truncated(Normal(f₁[t][1],1),0,Inf))
	end

	#Sample ξ₁ᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₁,Tₚ,ł,ρ²)
	ξ₁ᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(fₚ[t][1])))==1)
			ξ₁ᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(fₚ[t][1],1),-Inf,0)))
		end
	end

	#Sample ξ₀ₐᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,ł,ρ²)
	fₚ=FFBS2(f₀,Tₚ,ł,ρ²)
	ξ₀ₐᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(gₚ[t][1]))*Φ(fₚ[t][1])))==1)
			ξ₀ₐᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(fₚ[t][1],1),0,Inf)),"yg"=>rand(Truncated(Normal(gₚ[t][1],1),-Inf,0)))
		end
	end

	#sample ξ₀ᵣᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₀,Tₚ,ł,ρ²)
	ξ₀ᵣᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(fₚ[t][1])))==1)
			ξ₀ᵣᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(fₚ[t][1],1),-Inf,0)))
		end
	end

	#Sample ξ₁ₐᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,ł,ρ²)
	fₚ=FFBS2(f₁,Tₚ,ł,ρ²)
	ξ₁ₐᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli(Φ(gₚ[t][1])*Φ(fₚ[t][1])))==1)
			ξ₀ₐᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(fₚ[t][1],1),0,Inf)),"yg"=>rand(Truncated(Normal(gₚ[t][1],1),0,Inf)))
		end
	end

	#sample ξ₁ᵣᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₁,Tₚ,ł,ρ²)
	ξ₁ᵣᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(fₚ[t][1])))==1)
			ξ₁ᵣᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(fₚ[t][1],1),-Inf,0)))
		end
	end

	ξ₀ₐₐ=Dict{Float64,Dict{UTF8String,Float64}}()
	ξ₁ₐₐ=Dict{Float64,Dict{UTF8String,Float64}}()
	gₚ=FFBS2(g,Td,ł,ρ²)
	f₀ₚ=FFBS2(f₀,Td,ł,ρ²)
	f₁ₚ=FFBS2(f₁,Td,ł,ρ²)
	for t in Td
		denominator=Φ(gₚ[t][1])*Φ(f₀ₚ[t][1])+(1-Φ(gₚ[t][1]))*Φ(f₁ₚ[t][1])
		if(rand(Bernoulli((1-Φ(gₚ[t][1]))*Φ(f₁ₚ[t][1])/denominator))==1)
			ξ₁ₐₐ[t]=Dict("yf"=>rand(Truncated(Normal(f₁ₚ[t][1],1),0,Inf)), "yg"=>rand(Truncated(Normal(gₚ[t][1],1),-Inf,0)))
		else
			ξ₀ₐₐ[t]=Dict("yf"=>rand(Truncated(Normal(f₀ₚ[t][1],1),0,Inf)), "yg"=>rand(Truncated(Normal(gₚ[t][1],1),0,Inf)))
		end
	end

	yf0=Dict{Float64,Float64}()
	yf1=Dict{Float64,Float64}()
	yg=Dict{Float64,Float64}()
	for t in keys(ξ₀ₐ)
		yf0[t]=ξ₀ₐ[t]["yf"]
	end
	for t in keys(ξ₀ᵣ)
		yf0[t]=ξ₀ᵣ[t]["yf"]
	end
	for t in keys(ξ₁ₐ)
		yf1[t]=ξ₁ₐ[t]["yf"]
	end
	for t in keys(ξ₁ᵣ)
		yf1[t]=ξ₁ᵣ[t]["yf"]
	end
	for t in keys(ξ₀ₐₐ)
		yf0[t]=ξ₀ₐₐ[t]["yf"]
		yg[t]=ξ₀ₐₐ[t]["yg"]
	end
	for t in keys(ξ₀ₐᵣ)
		yf0[t]=ξ₀ₐᵣ[t]["yf"]
		yg[t]=ξ₀ₐᵣ[t]["yg"]
	end
	for t in keys(ξ₁ₐₐ)
		yf1[t]=ξ₁ₐₐ[t]["yf"]
		yg[t]=ξ₁ₐₐ[t]["yg"]
	end
	for t in keys(ξ₁ₐᵣ)
		yf1[t]=ξ₁ₐᵣ[t]["yf"]
		yg[t]=ξ₁ₐᵣ[t]["yg"]
	end
	for t in keys(ξ₀ᵣᵣ)
		yf0[t]=ξ₀ᵣᵣ[t]["yf"]
	end
	for t in keys(ξ₁ᵣᵣ)
		yf1[t]=ξ₁ᵣᵣ[t]["yf"]
	end
	g=FFBS(yg,collect(keys(yg)),μ,σ²,ł,ρ²)
	f₀=FFBS(yf0,collect(keys(yf0)),μ,σ²,ł,ρ²)
	f₁=FFBS(yf1,collect(keys(yf1)),μ,σ²,ł,ρ²)
	#=łₚ=ł+rand(Normal(0,0.03));	λₚ=sqrt(2*ν)/łₚ;	qₚ=2*ρ²*√π*λₚ^(2*p+1)*gamma(p+1)/gamma(p+0.1);=#
	#=if(log(rand(Uniform(0,1)))<glogdensity(g,gᵧ,λₚ,qₚ)-glogdensity(g,gᵧ,ł,ρ²))=#
		#=ł=łₚ=#
		#=λ=λₚ=#
		#=q=qₚ=#
		#=accepted=accepted+1=#
	#=end=#
	#=ρ²,q=rho(ρ²,g,gᵧ,ł,ρ²)=#
	#=n1=length(g)=#
	#=μ=rand(Normal((n1/σ²)*mean([y[key]-g[key][1] for key in keys(g)])*(1/((n1/σ²)+(1/σ²ₘ))),sqrt(1/((n1/σ²)+(1/σ²ₘ)))))=#
	#=ρ²Trace[iter]=ρ²=#
	#=łTrace[iter]=ł=#
	#=μTrace[iter]=μ=#
	#=logoddsTrace[iter]=logodds=#
	#=gᵧTrace[iter]=gᵧ=#
	#=plot(sort(collect(keys(g))),Φ([g[key][1] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);=#
	subplot(411)
	plot(convert(Array{Float64,1},sort(collect(keys(f₀)))),convert(Array{Float64},[Λ*Φ(f₀[key][1]) for key in sort(collect(keys(f₀)))]),c="grey",alpha=0.1)
	#=plot(sort(collect(keys(f₀))),([f₀[key][1] for key in sort(collect(keys(f₀)))]),c="grey",alpha=0.1);=#
	subplot(412)
	plot(convert(Array{Float64,1},sort(collect(keys(f₁)))),convert(Array{Float64},[Λ*Φ(f₁[key][1]) for key in sort(collect(keys(f₁)))]),c="grey",alpha=0.1)
	#=plot(sort(collect(keys(f₁))),([f₁[key][1] for key in sort(collect(keys(f₁)))]),c="grey",alpha=0.1);=#
	subplot(413)
	plot(convert(Array{Float64,1},sort(collect(keys(g)))),convert(Array{Float64},[Φ(g[key][1]) for key in sort(collect(keys(g)))]),c="grey",alpha=0.1)
	#=plot(sort(collect(keys(g))),([g[key][1] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);=#
	#=subplot(414)=#
	#=plot(convert(Array{Float64,1},sort(collect(keys(g)))),convert(Array{Float64},[Φ(g[key][1])*Λ*Φ(f₀[key][1])+(1-Φ(g[key][1]))*Λ*Φ(f₁[key][1]) for key in sort(collect(keys(g)))]),c="grey",alpha=0.1)=#
end
figure()
subplot(511)
plot(ρ²Trace)
subplot(512)
plot(μTrace)
subplot(513)
plot(łTrace)
subplot(514)
plot(gᵧTrace)
ylim(-0.5,1.5)
subplot(515)
plot(logoddsTrace)
mean(ρ²Trace)
mean(μTrace)
mean(łTrace)

#=tdataStream=open("tjuliaHigh.txt","w")=#
#=gdataStream=open("gjuliaHigh.txt","w")=#
#=gammadataStream=open("gammajuliaHigh.txt","w")=#
#=zgdataStream=open("zgjuliaHigh.txt","w")=#
#=writedlm(tdataStream,t)=#
#=writedlm(gdataStream,gTrace)=#
#=writedlm(zgdataStream,zgTrace)=#
#=writedlm(gammadataStream,γTrace)=#
#=close(tdataStream)=#
#=close(gdataStream)=#
#=close(zgdataStream)=#
#=close(gammadataStream)=#
