importall Base.Random
importall Distributions
using PDMats
using PyPlot
using IterativeSolvers

include("statespace.jl")
include("poissonpointprocess.jl")
include("gaussianprocess.jl")


srand(1)
disc=0.01
Tobs=1
grid=collect(0:disc:Tobs)
σ²=1.0
ł=0.1
p=3
d=p+1
ν=p+0.5
λ=sqrt(2*ν)/ł
q=2*σ²*√π*λ^(2*p+1)*gamma(p+1)/gamma(p+0.5)
g=chol(Kernel(grid,grid,1.0,1.0)+0.000000001*eye(length(grid)))'*rand(Normal(0,1),length(grid))
gfun(arg)=g[maximum(find( y->(y <= arg), grid))]
Λ=1500
μ=1.0
σ²ₘ=1
β=0.4
σ²ᵦ=1
intensity(x)=Λ*Φ(μ+β*gfun(x))
plot(grid,map(intensity,grid),c="green")
#=plot(grid,Λ*ones(length(grid)),c="black")=#
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

niter=500
łTrace=Array{Float64}(niter)
μTrace=Array{Float64}(niter)
βTrace=Array{Float64}(niter)
βᵧTrace=Array{Float64}(niter)
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
		if(rand(Bernoulli((1-Φ(μ+β*gₚ[t][1]))))==1)
			ξ₀₁₀[t]=Dict{UTF8String,Float64}("Zg"=>rand(Truncated(Normal(μ+β*gₚ[t][1],1),-Inf,0)), "γ"=>0)
		end
	end
	for t in keys(ξ)
		ξ[t]["Zg"]=rand(Truncated(Normal(μ+β*g[t][1],1),0,Inf))
	end
	Zgc=Dict{Float64,Float64}()
	for t in keys(ξ)
		Zgc[t]=ξ[t]["Zg"]
	end
	for t in keys(ξ₀₁₀)
		Zgc[t]=ξ₀₁₀[t]["Zg"]
	end
	łₚ=ł+rand(Normal(0,0.05))
	λₚ=sqrt(2*ν)/łₚ
	qₚ=2*σ²*√π*λₚ^(2*p+1)*gamma(p+1)/gamma(p+0.5)
	if(log(rand(Uniform(0,1)))<sslogdensity(Zgc,μ,β,λₚ,qₚ)-sslogdensity(Zgc,μ,β,λ,q))
		ł=łₚ
		λ=λₚ
		q=qₚ
	end
	g=FFBS(Zgc,collect(keys(Zgc)),μ,β,σ²,λ,q)
	plot(sort(collect(keys(g))),Λ*Φ([μ+β*g[key][1] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);
	g1=[g[i][1] for i in sort(collect(keys(g)))]
	Zgc1=[Zgc[i] for i in sort(collect(keys(Zgc)))]
	n1=length(g1)
	μ=rand(Normal(sum(Zgc1-β*g1)*(1/((n1/σ²)+(1/σ²ₘ))),sqrt(1/((n1/σ²)+(1/σ²ₘ)))))
	odds=sqrt(1/(dot(g1,g1)*σ²ᵦ+1))*exp(0.5*(σ²ᵦ/(σ²ᵦ*dot(g1,g1)+1))*(dot(Zgc1-μ,g1)^2))*Φ(dot(Zgc1-μ,g1)*σ²ᵦ/(σ²ᵦ*dot(g1,g1)+1))
	if(odds==Inf)
		βᵧ=1
	elseif(odds==-Inf)
		βᵧ=0
	else
		βᵧ=rand(Bernoulli(odds/(1+odds)))
	end
	if(βᵧ==1)
		β=rand(Truncated(Normal(dot(g1,Zgc1-μ)/(dot(g1,g1)+(1/σ²ᵦ)),1/sqrt(dot(g1,g1)+(1/σ²ᵦ))),0,Inf))
	else
		β=0.0
	end
	βᵧTrace[iter]=βᵧ
	βTrace[iter]=β
	łTrace[iter]=ł
	μTrace[iter]=μ
end
