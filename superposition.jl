importall Base.Random
importall Distributions
using PDMats
using PyPlot

include("spikeandslab.jl")
include("poissonpointprocess.jl")
include("gaussianprocess.jl")
include("cleandataGeneration.jl")
#=include("dataGeneration.jl")=#
#=plot(grid,alpha,c="green")=#

niter=50
ZgTrace=Array(Float64,length(ξ),niter)
gTrace=Array(Float64,length(ξ),niter)
γTrace=Array(Float64,length(ξ),niter)
ψ²Trace=Array(Float64,niter)
gConditional=GP(g,σ²,ρ²,ψ²,"function")
Ξₚ=PPProcess(Λ)
σ²=1.0
prior=Dict(:ψ²=>SpikeAndSlab([0.0,1.0],Gamma(5,10)))
proposal=SpikeAndSlab([0.0,1.0],Normal(0,10))
for iter=1:niter
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=rand(gConditional,Tₚ)
	ξ₀₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(gₚ[t]))*λ₀(t)/Ξₚ.λ(t)))==1)
			ξ₀₁₀[t]=Dict{UTF8String,Float64}("Zg"=>rand(Truncated(Normal(gₚ[t],1),-Inf,0)), "γ"=>0)
		end
	end
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=rand(gConditional,Tₚ)
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
	gProc=GP(Zgc,σ²,ρ²,ψ²,"response")
	ψ²ₚ=rand(SpikeAndSlab([0.0,1.0],Normal(ψ²,10)))
	if(ψ²ₚ>=0.0)
		gProcProp=GP(Zgc,σ²,ρ²,ψ²ₚ,"response")
		#=if(log(rand(Uniform(0,1)))<procprop.logDensity+log(pdf(prior[:ψ²],ψ²ₚ)[1])+log(q(ψ²,ψ²ₚ,proposal))-procprop.logDensity-log(pdf(prior[:ψ²],ψ²)[1])-log(q(ψ²ₚ,ψ²,proposal)))=#
		if(log(rand(Uniform(0,1)))<gProcProp.logDensity+log(pdf(prior[:ψ²],ψ²ₚ))+log(q(ψ²,ψ²ₚ,proposal))-gProc.logDensity-log(pdf(prior[:ψ²],ψ²))-log(q(ψ²ₚ,ψ²,proposal)))
			ψ²=ψ²ₚ
			gProc=GP(Zgc,σ²,ρ²,ψ²,"response")
		end
	end
	g=rand(gProc)
	gConditional=GP(g,σ²,ρ²,ψ²,"function")
	ψ²Trace[iter]=ψ²
	#=gTrace[:,iter]=[g[t] for t in keys(ξ)]=#
	#=ZgTrace[:,iter]=[ξ[key]["Zg"] for key in keys(ξ)]=#
	#=γTrace[:,iter]=[ξ[key]["γ"] for key in keys(ξ)]=#
	plot(sort(collect(keys(g))),Φ([g[key] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);
end

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
