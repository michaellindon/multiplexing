using Distributions
using PyPlot

include("functions.jl")
include("cleandataGeneration.jl")
#=include("dataGeneration.jl")=#
#=plot(grid,alpha,c="green")=#

niter=10
ZgTrace=Array(Float64,length(ξ),niter)
gTrace=Array(Float64,length(ξ),niter)
γTrace=Array(Float64,length(ξ),niter)
Ξₚ=PPProcess(Λ)
for iter=1:niter
	Tₚ=realization(Ξₚ,0,Tobs)
	gₚ=PredictGP(Tₚ,g,1.0,ρ²,ψ²,"function")
	ξ₀₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(gₚ[t]))*λ₀(t)/Ξₚ.λ(t)))==1)
			ξ₀₁₀[t]=Dict{UTF8String,Float64}("Zg"=>rand(Truncated(Normal(gₚ[t],1),-Inf,0)), "γ"=>0)
		end
	end
	Tₚ=realization(Ξₚ,0,Tobs)
	gₚ=PredictGP(Tₚ,g,1.0,ρ²,ψ²,"function")
	ξ₁₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli(Φ(gₚ[t])*λ₁(t)/Ξₚ.λ(t)))==1)
			ξ₁₁₀[t]=Dict{UTF8String,Float64}("Zg"=>rand(Truncated(Normal(gₚ[t],1),0,Inf)),"γ"=>1)
		end
	end
	for t in keys(ξ)
		denominator=Φ(g[t])*λ₀(t)+(1-Φ(g[t]))*λ₁(t)
		ξ[t]["γ"]=rand(Bernoulli((1-Φ(g[t]))*λ₁(t)/denominator))
		if(ξ[t]["γ"]==1)
			ξ[t]["Zg"]=rand(Truncated(Normal(g[t],1),-Inf,0))
		else
			ξ[t]["Zg"]=rand(Truncated(Normal(g[t],1),0,Inf))
		end
	end
	Zgc=Dict{Float64,Float64}()
	for t in keys(ξ)
		Zgc[t]=ξ[t]["Zg"]
	end
	for t in keys(ξ₁₁₀)
		Zgc[t]=ξ₁₁₀[t]["Zg"]
	end
	for t in keys(ξ₀₁₀)
		Zgc[t]=ξ₀₁₀[t]["Zg"]
	end
	g=PosteriorGP(Zgc,1.0,ρ²,ψ²)
	gTrace[:,iter]=[g[t] for t in keys(ξ)]
	ZgTrace[:,iter]=[ξ[key]["Zg"] for key in keys(ξ)]
	γTrace[:,iter]=[ξ[key]["γ"] for key in keys(ξ)]
	plot(sort(collect(keys(g))),Φ([g[key] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);
end

tdataStream=open("tjuliaHigh.txt","w")
gdataStream=open("gjuliaHigh.txt","w")
gammadataStream=open("gammajuliaHigh.txt","w")
zgdataStream=open("zgjuliaHigh.txt","w")
writedlm(tdataStream,t)
writedlm(gdataStream,gTrace)
writedlm(zgdataStream,zgTrace)
writedlm(gammadataStream,γTrace)
close(tdataStream)
close(gdataStream)
close(zgdataStream)
close(gammadataStream)
