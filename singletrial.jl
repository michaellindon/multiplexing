importall Base.Random
importall Distributions
using PDMats
using PyPlot
using IterativeSolvers

include("statespace.jl")
include("poissonpointprocess.jl")

srand(1)
ρ²=0.9; ł=0.1; p=3; d=p+1; ν=p+0.5; λ=sqrt(2*ν)/ł; q=2*ρ²*√π*λ^(2*p+1)*gamma(p+1)/gamma(p+0.5);
Tobs=1
t=Dict(zip(1:101,collect(0:0.01:Tobs)))
g=Dict();A=Dict();Q=Dict();Δ=Dict()
g[t[1]]=rand(MvNormal(statcov(λ,q)))
for i=2:length(t)
	Δ[i]=t[i]-t[i-1]
	A[t[i-1]]=transition(Δ[i],λ)
	Q[t[i-1]]=innovation(Δ[i],λ,q)
	g[t[i]]=A[t[i-1]]*g[t[i-1]]+rand(MvNormal(Q[t[i-1]]+0.00000001*eye(d)))
end
grid=convert(Array{Float64},sort(collect(keys(g))))
g=[g[key][1] for key in sort(collect(keys(g)))]
gfun(arg)=g[maximum(find( y->(y <= arg), grid))]
Λ=1000; μ=1.0; intensity(x)=Λ*Φ(μ+gfun(x));
plot(grid,map(intensity,grid),c="green")
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
plot(accepted,-0*ones(length(accepted)),c="green",marker="|",linestyle="None",markersize=20)
plot(rejected,-140*ones(length(rejected)),c="red",marker="|",linestyle="None",markersize=20)

niter=800
gᵧTrace=Array{Float64}(niter)
łTrace=Array{Float64}(niter)
μTrace=Array{Float64}(niter)
ρ²Trace=Array{Float64}(niter)
logoddsTrace=Array{Float64}(niter)
ξ=Dict{Float64,Dict{UTF8String,Float64}}()
for t in accepted
	ξ[t]=Dict("y"=>rand(Truncated(Normal(0,1),0,Inf)),"γ"=>0)
end
g=Dict()
for foo in collect(keys(ξ))
	g[foo]=zeros(4,1)
end
Ξₚ=PPProcess(Λ)
σ²=1.0
gᵧ=1
r=1
σ²ₘ=5
μ=0.0
ρ²=0.1
q=2*ρ²*√π*λ^(2*p+1)*gamma(p+1)/gamma(p+0.5);
for iter=1:niter
	if(iter%100==0)
		println(iter)
	end
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,λ,q)
	ξ₀₁₀=Dict{Float64,Dict{UTF8String,Float64}}()
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(μ+gₚ[t][1]))))==1)
			ξ₀₁₀[t]=Dict{UTF8String,Float64}("y"=>rand(Truncated(Normal(μ+gₚ[t][1],1),-Inf,0)), "γ"=>0)
		end
	end
	for t in keys(ξ)
		ξ[t]["y"]=rand(Truncated(Normal(μ+g[t][1],1),0,Inf))
	end
	y=Dict{Float64,Float64}()
	for t in keys(ξ)
		y[t]=ξ[t]["y"]
	end
	for t in keys(ξ₀₁₀)
		y[t]=ξ₀₁₀[t]["y"]
	end
	logodds=(sslogdensity2(y,1,μ,1,σ²,λ,q)-sslogdensity2(y,0,μ,1,σ²,λ,q))
	odds=exp(logodds)
	if(odds==Inf)
		gᵧ=1
	elseif(odds==-Inf)
		gᵧ=0
	else
		gᵧ=rand(Bernoulli(odds/(1+odds)))
	end
	if(gᵧ==1)
		g=FFBS(y,collect(keys(y)),μ,1,σ²,λ,q)
	else
		g=Dict()
		for key in keys(y)
			g[key]=zeros(4,1)
		end
	end
	#=for key in keys(g)=#
		#=g[key][2]=gtrue[key][2]=#
		#=g[key][3]=gtrue[key][3]=#
		#=g[key][4]=gtrue[key][4]=#
	#=end=#
	#=plot([g[key][4] for key in sort(collect(keys(g)))],c="red",alpha=0.1)=#
	łₚ=ł+rand(Normal(0,0.03));	λₚ=sqrt(2*ν)/łₚ;	qₚ=2*ρ²*√π*λₚ^(2*p+1)*gamma(p+1)/gamma(p+0.5);
	if(log(rand(Uniform(0,1)))<glogdensity(g,gᵧ,λₚ,qₚ)-glogdensity(g,gᵧ,λ,q))
		ł=łₚ
		λ=λₚ
		q=qₚ
		accepted=accepted+1
	end
	ρ²,q=rho(ρ²,g,gᵧ,λ,q)
	n1=length(g)
	μ=rand(Normal((n1/σ²)*mean([y[key]-g[key][1] for key in keys(g)])*(1/((n1/σ²)+(1/σ²ₘ))),sqrt(1/((n1/σ²)+(1/σ²ₘ)))))
	ρ²Trace[iter]=ρ²
	łTrace[iter]=ł
	μTrace[iter]=μ
	logoddsTrace[iter]=logodds
	gᵧTrace[iter]=gᵧ
	plot(sort(collect(keys(g))),Λ*Φ([μ+g[key][1] for key in sort(collect(keys(g)))]),c="blue",alpha=0.1);
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



srand(1)
r=2
σ²=1.0
gᵧ=1
niter=1000
gᵧTrace=zeros(niter)
ρ²Trace=zeros(niter)
łTrace=zeros(niter)
μTrace=zeros(niter)
logoddsTrace=zeros(niter)
ρ²=0.1; ł=2; p=3; d=p+1; ν=p+0.5; λ=sqrt(2*ν)/ł; q=2*ρ²*√π*λ^(2*p+1)*gamma(p+1)/gamma(p+0.5); μ=5; σ²ₘ=10; gᵧ=1;
Tobs=100
disc=collect(0:0.1:Tobs)
t=Dict(zip(1:length(disc),disc))
g=Dict();A=Dict();Q=Dict();Δ=Dict()
g[t[1]]=rand(MvNormal(statcov(λ,q)))
y=Dict()
y[t[1]]=μ+g[t[1]][1]+rand(Normal(0,sqrt(σ²)))
for i=2:length(t)
	Δ[i]=t[i]-t[i-1]
	A[t[i-1]]=transition(Δ[i],λ)
	Q[t[i-1]]=innovation(Δ[i],λ,q)
	g[t[i]]=A[t[i-1]]*g[t[i-1]]+rand(MvNormal(Q[t[i-1]]+0.00000001*eye(d)))
	y[t[i]]=μ+g[t[i]][1]+rand(Normal(0,sqrt(σ²)))
end
#=plot([y[key] for key in sort(collect(keys(y)))])=#
#=gtrue=copy(g)=#
#=μ=-3=#
#=ρ²=10=#
#=q=2*ρ²*√π*λ^(2*p+1)*gamma(p+1)/gamma(p+0.5); =#
iter=1
accepted=0
for iter=1:niter
	if(iter%100==0)
		println(iter)
	end
	logodds=(sslogdensity2(y,1,μ,1,σ²,λ,q)-sslogdensity2(y,0,μ,1,σ²,λ,q))
	odds=exp(logodds)
	if(odds==Inf)
		gᵧ=1
	elseif(odds==-Inf)
		gᵧ=0
	else
		gᵧ=rand(Bernoulli(odds/(1+odds)))
	end
	if(gᵧ==1)
		g=FFBS(y,collect(keys(y)),μ,1,σ²,λ,q)
	else
		g=Dict()
		for key in keys(y)
			g[key]=zeros(4,1)
		end
	end
	#=for key in keys(g)=#
		#=g[key][2]=gtrue[key][2]=#
		#=g[key][3]=gtrue[key][3]=#
		#=g[key][4]=gtrue[key][4]=#
	#=end=#
	#=plot([g[key][4] for key in sort(collect(keys(g)))],c="red",alpha=0.1)=#
	łₚ=ł+rand(Normal(0,0.03));	λₚ=sqrt(2*ν)/łₚ;	qₚ=2*ρ²*√π*λₚ^(2*p+1)*gamma(p+1)/gamma(p+0.5);
	if(log(rand(Uniform(0,1)))<glogdensity(g,gᵧ,λₚ,qₚ)-glogdensity(g,gᵧ,λ,q))
		ł=łₚ
		λ=λₚ
		q=qₚ
		accepted=accepted+1
	end
	ρ²,q=rho(ρ²,g,gᵧ,λ,q)
	n1=length(g)
	μ=rand(Normal((n1/σ²)*mean([y[key]-g[key][1] for key in keys(g)])*(1/((n1/σ²)+(1/σ²ₘ))),sqrt(1/((n1/σ²)+(1/σ²ₘ)))))
	ρ²Trace[iter]=ρ²
	łTrace[iter]=ł
	μTrace[iter]=μ
	logoddsTrace[iter]=logodds
	gᵧTrace[iter]=gᵧ
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





E,V=eig(Q[t[i]])
logpdf(MvNormal(vec(A[t[i]]*g[t[i]]),Q[t[i]]),g[t[i+1]])
E.>0.001
subind=zeros(Int64,4)
for w=1:length(E)
	if(E[w]<sqrt(eps(real(float(one(eltype(Q[t[i]])))))))
		E[w]=Inf
	else
		subind[w]=1
	end
end
subindices=find(x->x==1,subind)
-0.5*length(subindices)*log(2*π) - 0.5*logdet(Q[t[i]][subindices,subindices]) - (g[t[i+1]]-A[t[i]]*g[t[i]])'*V*pinv(Diagonal(E))*V'*(g[t[i+1]]-A[t[i]]*g[t[i]])
