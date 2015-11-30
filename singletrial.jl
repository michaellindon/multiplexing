importall Base.Random
importall Distributions
using PDMats
using PyPlot
using IterativeSolvers

include("statespace.jl")
include("spikeandslab.jl")
include("poissonpointprocess.jl")
include("gaussianprocess.jl")


srand(30)
disc=0.01
Tobs=1
grid=collect(0:disc:Tobs)
σ²=1.0
ł=0.1
p=3
ν=p+0.5
λ=sqrt(2*ν)/ł
g=chol(Kernel(grid,grid,ρ²,ψ²)+0.000000001*eye(length(grid)))'*rand(Normal(0,1),length(grid))
gfun(arg)=g[maximum(find( y->(y <= arg), grid))]
Λ=1000
intensity(x)=Λ*Φ(gfun(x))
plot(grid,map(intensity,grid),c="green")
plot(grid,Λ*ones(length(grid)),c="black")
nₚ=rand(Poisson(Λ*(Tobs)))
Tₚ=sort(rand(Uniform(0,Tobs),nₚ))
accepted=Array(Float64,0)
rejected=Array(Float64,0)
for t in Tₚ
	if(rand(Bernoulli(intensity(t)/Λ))==1)
		push!(accepted,t)
	else
		push!(rejected,t)
	end
end
plot(Tₚ,-1*ones(length(Tₚ)),c="black",marker="|",linestyle="None",markersize=20)
plot(accepted,-0*ones(length(accepted)),c="green",marker="|",linestyle="None",markersize=20)
plot(rejected,-140*ones(length(rejected)),c="red",marker="|",linestyle="None",markersize=20)

niter=20
ξ=Dict{Float64,Dict{UTF8String,Float64}}()
for t in accepted
	ξ[t]=Dict("Zg"=>rand(Truncated(Normal(0,1),0,Inf)),"γ"=>0)
end
g=Dict()
for foo in collect(keys(ξ))
	g[foo]=rand(4,1)
end
Ξₚ=PPProcess(Λ)
σ²=1.0
for iter=1:niter
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,λ,q)
	ξ₀₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(gₚ[t][1]))))==1)
			ξ₀₁₀[t]=Dict{UTF8String,Float64}("Zg"=>rand(Truncated(Normal(gₚ[t][1],1),-Inf,0)), "γ"=>0)
		end
	end
	for t in keys(ξ)
		ξ[t]["Zg"]=rand(Truncated(Normal(g[t][1],1),0,Inf))
	end
	Zgc=Dict{Float64,Float64}()
	for t in keys(ξ)
		Zgc[t]=ξ[t]["Zg"]
	end
	for t in keys(ξ₀₁₀)
		Zgc[t]=ξ₀₁₀[t]["Zg"]
	end
	g=FFBS(Zgc,collect(keys(Zgc)),σ²,λ,q)
	plot(sort(collect(keys(g))),Λ*Φ([g[key][1] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);
end
