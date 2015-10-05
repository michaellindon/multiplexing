using Distributions
#=using PyPlot=#
using DataFrames

include("dataGeneration.jl")
include("functions.jl")
#=plot(x,alpha,c="green")=#

ξ011=convert(Array,mppp[(mppp[:label].=="011"),[:t,:Zg]])
ξ011=hcat(ξ011,zeros(size(ξ011,1)))
ξ111=convert(Array,mppp[(mppp[:label].=="111"),[:t,:Zg]])
ξ111=hcat(ξ111,ones(size(ξ111,1)))
ξ=vcat(ξ011,ξ111)
indices=sortperm(ξ[:,1])
ξ=ξ[indices,:]
ξ010=convert(Array,mppp[(mppp[:label].=="010"),[:t,:Zg]])
ξ110=convert(Array,mppp[(mppp[:label].=="110"),[:t,:Zg]])
t=convert(Array,mppp[(mppp[:label].=="011")|(mppp[:label].=="111"),:t])
g=convert(Array,mppp[(mppp[:label].=="011")|(mppp[:label].=="111"),:g])
tc=convert(Array,mppp[(mppp[:label].=="011")|(mppp[:label].=="111")|(mppp[:label].=="010")|(mppp[:label].=="110"),:t])
gc=convert(Array,mppp[(mppp[:label].=="011")|(mppp[:label].=="111")|(mppp[:label].=="010")|(mppp[:label].=="110"),:g])
ξ=hcat(t,convert(Array,mppp[(mppp[:label].=="011")|(mppp[:label].=="111"),:Zg]),ones(length(t)))
for i=1:size(ξ,1)
	if((mppp[(mppp[:label].=="011")|(mppp[:label].=="111"),:label])[i]=="011")
		ξ[i,3]=0
	else
		ξ[i,3]=1
	end
end


niter=40
zgTrace=Array(Float64,length(t),niter)
gTrace=Array(Float64,length(t),niter)
γTrace=Array(Float64,length(t),niter)
for iter=1:niter
	np=rand(Poisson(Λ*T))
	tp=sort(rand(Uniform(0,T),np))
	ξ010=Array(Float64,0,2)
	yp,gp=PredictGP(tp,gc,tc,1.0,ρ²,ψ²,"function")
	for i=1:np
		if(rand(Bernoulli((1-Φ(gp[i]))*λ₀(tp[i])/Λ),1)[1]==1)
			zg=rand(Truncated(Normal(gp[i],1),-Inf,0),1)[1]
			ξ010=vcat(ξ010,[tp[i],zg]')
		end
	end
	np=rand(Poisson(Λ*T))
	tp=sort(rand(Uniform(0,T),np))
	ξ110=Array(Float64,0,2)
	yp,gp=PredictGP(tp,gc,tc,1.0,ρ²,ψ²,"function")
	for i=1:np
		if(rand(Bernoulli(Φ(gp[i])*λ₁(tp[i])/Λ),1)[1]==1)
			zg=rand(Truncated(Normal(gp[i],1),0,Inf),1)[1]
			ξ110=vcat(ξ110,[tp[i],zg]')
		end
	end
	for i=1:size(ξ,1)
		denominator=Φ(g[i])*λ₀(t[i])+(1-Φ(g[i]))*λ₁(t[i])
		ξ[i,3]=rand(Bernoulli((1-Φ(g[i]))*λ₁(t[i])/denominator))[1]
		if(ξ[i,3]==1)
			ξ[i,2]=rand(Truncated(Normal(g[i],1),-Inf,0),1)[1]
		else
			ξ[i,2]=rand(Truncated(Normal(g[i],1),0,Inf),1)[1]
		end
	end
	tc=[ξ[:,1],ξ010[:,1],ξ110[:,1]]
	zgc=[ξ[:,2],ξ010[:,2],ξ110[:,2]]
	gc=PosteriorGP(zgc,tc,1.0,ρ²,ψ²)
	g=gc[1:size(ξ,1)]
	gTrace[:,iter]=g
	zgTrace[:,iter]=ξ[:,2]
	γTrace[:,iter]=ξ[:,3]
	indices=sortperm(tc)
	tcPlot=tc[indices]
	gcPlot=gc[indices]
	plot(tcPlot,Φ(gcPlot),c="grey",alpha=0.1);
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
