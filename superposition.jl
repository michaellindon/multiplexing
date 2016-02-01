importall Base.Random
importall Distributions
using PyPlot

include("poissonpointprocess.jl")
include("statespace.jl")
include("cleandataGeneration.jl")
#=include("dataGeneration.jl")=#
#=plot(grid,alpha,c="green")=#

niter=2000
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
Ξₚ=PPProcess(Λ)
σ²=1.0
for iter=1:niter
	println(iter)
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,λ,q)
	ξ₀₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(gₚ[t][1]))*λ₀(t)/Ξₚ.λ(t)))==1)
			ξ₀₁₀[t]=Dict{UTF8String,Float64}("y"=>rand(Truncated(Normal(gₚ[t][1],1),-Inf,0)), "γ"=>0)
		end
	end
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,λ,q)
	ξ₁₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli(Φ(gₚ[t][1])*λ₁(t)/Ξₚ.λ(t)))==1)
			ξ₁₁₀[t]=Dict{UTF8String,Float64}("y"=>rand(Truncated(Normal(gₚ[t][1],1),0,Inf)),"γ"=>1)
		end
	end
	for t in keys(ξ)
		denominator=Φ(g[t][1])*λ₀(t)+(1-Φ(g[t][1]))*λ₁(t)
		ξ[t]["γ"]=rand(Bernoulli((1-Φ(g[t][1]))*λ₁(t)/denominator))
		if(ξ[t]["γ"]==1)
			ξ[t]["y"]=rand(Truncated(Normal(g[t][1],1),-Inf,0))
		else
			ξ[t]["y"]=rand(Truncated(Normal(g[t][1],1),0,Inf))
		end
	end
	y=Dict{Float64,Float64}()
	for t in keys(ξ)
		y[t]=ξ[t]["y"]
	end
	for t in keys(ξ₁₁₀)
		y[t]=ξ₁₁₀[t]["y"]
	end
	for t in keys(ξ₀₁₀)
		y[t]=ξ₀₁₀[t]["y"]
	end
	g=FFBS(y,collect(keys(y)),μ,1,σ²,λ,q)
	#=łₚ=ł+rand(Normal(0,0.03));	λₚ=sqrt(2*ν)/łₚ;	qₚ=2*ρ²*√π*λₚ^(2*p+1)*gamma(p+1)/gamma(p+0.5);=#
	#=if(log(rand(Uniform(0,1)))<glogdensity(g,gᵧ,λₚ,qₚ)-glogdensity(g,gᵧ,λ,q))=#
		#=ł=łₚ=#
		#=λ=λₚ=#
		#=q=qₚ=#
		#=accepted=accepted+1=#
	#=end=#
	#=ρ²,q=rho(ρ²,g,gᵧ,λ,q)=#
	#=n1=length(g)=#
	#=μ=rand(Normal((n1/σ²)*mean([y[key]-g[key][1] for key in keys(g)])*(1/((n1/σ²)+(1/σ²ₘ))),sqrt(1/((n1/σ²)+(1/σ²ₘ)))))=#
	ρ²Trace[iter]=ρ²
	łTrace[iter]=ł
	μTrace[iter]=μ
	logoddsTrace[iter]=logodds
	gᵧTrace[iter]=gᵧ
	plot(sort(collect(keys(g))),Φ([g[key][1] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);
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
