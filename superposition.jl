importall Base.Random
importall Distributions
using PDMats
using PyPlot
using IterativeSolvers

include("spikeandslab.jl")
include("poissonpointprocess.jl")
include("gaussianprocess.jl")
include("statespace.jl")
include("cleandataGeneration.jl")
#=include("dataGeneration.jl")=#
#=plot(grid,alpha,c="green")=#

niter=20
Ξₚ=PPProcess(Λ)
σ²=1.0
for iter=1:niter
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,λ,q)
	ξ₀₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(gₚ[t][1]))*λ₀(t)/Ξₚ.λ(t)))==1)
			ξ₀₁₀[t]=Dict{UTF8String,Float64}("Zg"=>rand(Truncated(Normal(gₚ[t][1],1),-Inf,0)), "γ"=>0)
		end
	end
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,λ,q)
	ξ₁₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli(Φ(gₚ[t][1])*λ₁(t)/Ξₚ.λ(t)))==1)
			ξ₁₁₀[t]=Dict{UTF8String,Float64}("Zg"=>rand(Truncated(Normal(gₚ[t][1],1),0,Inf)),"γ"=>1)
		end
	end
	for t in keys(ξ)
		denominator=Φ(g[t][1])*λ₀(t)+(1-Φ(g[t][1]))*λ₁(t)
		ξ[t]["γ"]=rand(Bernoulli((1-Φ(g[t][1]))*λ₁(t)/denominator))
		if(ξ[t]["γ"]==1)
			ξ[t]["Zg"]=rand(Truncated(Normal(g[t][1],1),-Inf,0))
		else
			ξ[t]["Zg"]=rand(Truncated(Normal(g[t][1],1),0,Inf))
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
	g=FFBS(Zgc,collect(keys(Zgc)),σ²,λ,q)
	#=ψ²ₚ=rand(SpikeAndSlab([0.0,1.0],Normal(ψ²,10)))=#
	#=if(ψ²ₚ>=0.0)=#
		#=gProcProp=GP(Zgc,σ²,ρ²,ψ²ₚ,"response")=#
		#=[>if(log(rand(Uniform(0,1)))<procprop.logDensity+log(pdf(prior[:ψ²],ψ²ₚ)[1])+log(q(ψ²,ψ²ₚ,proposal))-procprop.logDensity-log(pdf(prior[:ψ²],ψ²)[1])-log(q(ψ²ₚ,ψ²,proposal)))<]=#
		#=if(log(rand(Uniform(0,1)))<gProcProp.logDensity+log(pdf(prior[:ψ²],ψ²ₚ))+log(q(ψ²,ψ²ₚ,proposal))-gProc.logDensity-log(pdf(prior[:ψ²],ψ²))-log(q(ψ²ₚ,ψ²,proposal)))=#
			#=ψ²=ψ²ₚ=#
			#=gProc=GP(Zgc,σ²,ρ²,ψ²,"response")=#
		#=end=#
	#=end=#
	#=gTrace[:,iter]=[g[t] for t in keys(ξ)]=#
	#=ZgTrace[:,iter]=[ξ[key]["Zg"] for key in keys(ξ)]=#
	#=γTrace[:,iter]=[ξ[key]["γ"] for key in keys(ξ)]=#
	plot(sort(collect(keys(g))),Φ([g[key][1] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);
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
