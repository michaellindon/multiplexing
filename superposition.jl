using Distributions
using PyPlot
using DataFrames

include("dataGeneration.jl")
plot(x,alpha,c="green")


nc=size(mppp)[1]
nα=length(find(mppp[:mf].==1))
for iter=1:20
	ixα=find(mppp[:mf].==1)
	nα=length(ixα)
	K=Array(Float64,nα,nα);
	for i=1:nα
		for j=1:nα
			K[i,j]=rho2*exp(-psi2*(mppp[ixα[i],:t]-mppp[ixα[j],:t])^2)
		end
	end
	V=\(K+eye(nα),K)
	L=svd(V)
	L=L[1]*Diagonal(sqrt(L[2]))
	np=rand(Poisson(λc*T))
	tp=sort(rand(Uniform(0,T),np))
	mppp010=DataFrame(t=tp,γ=zeros(Int64,np),mf=zeros(Int64,np),mg=Array(Int64,np),Zf=Array(Float64,np),Zg=Array(Float64,np),label=Array(String,np))
	Kpp=Array(Float64,np,np);
	for i=1:np
		for j=1:np
			Kpp[i,j]=rho2*exp(-psi2*(mppp010[i,:t]-mppp010[j,:t])^2)
		end
	end
	Kpα=Array(Float64,np,nα);
	for i=1:np
		for j=1:nα
			Kpα[i,j]=rho2*exp(-psi2*(mppp010[i,:t]-mppp[ixα[j],:t])^2)
		end
	end
	Lp=svd(Kpp-Kpα*\(K,Kpα'))
	Lp=Lp[1]*Diagonal(sqrt(Lp[2]))
	mppp010[:g]=Kpα*\(K,convert(Array,mppp[ixα,:g]))+Lp*rand(Normal(0,1),np);
	for i=1:np
		mppp010[i,:Zf]=rand(Normal(f₀(mppp[i,:t]),1))
		if(mppp010[i,:Zf]>0)
			mppp010[i,:mf]=1
			mppp010[i,:Zg]=rand(Normal(mppp010[i,:g],1))
			if(mppp010[i,:Zg]>0)
				mppp010[i,:mg]=1
				mppp010[i,:label]=string(mppp010[i,:γ],mppp010[i,:mf],mppp010[i,:mg])
			else
				mppp010[i,:mg]=0
				mppp010[i,:label]=string(mppp010[i,:γ],mppp010[i,:mf],mppp010[i,:mg])
			end
		else 
			mppp010[i,:mf]=0
			mppp010[i,:mg]=NA
			mppp010[i,:Zg]=NA
			mppp010[i,:label]=string(mppp010[i,:γ],mppp010[i,:mf],mppp010[i,:mg])
		end
	end
	mppp010=mppp010[mppp010[:label].=="010",:]
	np=rand(Poisson(λc*T))
	tp=sort(rand(Uniform(0,T),np))
	mppp110=DataFrame(t=tp,γ=ones(Int64,np),mf=zeros(Int64,np),mg=Array(Int64,np),Zf=Array(Float64,np),Zg=Array(Float64,np),label=Array(String,np))
	Kpp=Array(Float64,np,np);
	for i=1:np
		for j=1:np
			Kpp[i,j]=rho2*exp(-psi2*(mppp110[i,:t]-mppp110[j,:t])^2)
		end
	end
	Kpα=Array(Float64,np,nα);
	for i=1:np
		for j=1:nα
			Kpα[i,j]=rho2*exp(-psi2*(mppp110[i,:t]-mppp[ixα[j],:t])^2)
		end
	end
	Lp=svd(Kpp-Kpα*\(K,Kpα'))
	Lp=Lp[1]*Diagonal(sqrt(Lp[2]))
	mppp110[:g]=Kpα*\(K,convert(Array,mppp[ixα,:g]))+Lp*rand(Normal(0,1),np);
	for i=1:np
		mppp110[i,:Zf]=rand(Normal(f₁(mppp[i,:t]),1))
		if(mppp110[i,:Zf]>0)
			mppp110[i,:mf]=1
			mppp110[i,:Zg]=rand(Normal(mppp110[i,:g],1))
			if(mppp110[i,:Zg]<0)
				mppp110[i,:mg]=1
				mppp110[i,:label]=string(mppp110[i,:γ],mppp110[i,:mf],mppp110[i,:mg])
			else
				mppp110[i,:mg]=0
				mppp110[i,:label]=string(mppp110[i,:γ],mppp110[i,:mf],mppp110[i,:mg])
			end
		else 
			mppp110[i,:mf]=0
			mppp110[i,:mg]=NA
			mppp110[i,:Zg]=NA
			mppp110[i,:label]=string(mppp110[i,:γ],mppp110[i,:mf],mppp110[i,:mg])
		end
	end
	mppp110=mppp110[mppp110[:label].=="110",:]
	ixos=find((mppp[:label].=="011")|(mppp[:label].=="111"))
	for ix in ixos
		numerator=cdf(Normal(0,1),mppp[ix,:g])*λ₀(mppp[ix,:t])+(1-cdf(Normal(0,1),mppp[ix,:g]))*λ₁(mppp[ix,:t])
		mppp[ix,:γ]=rand(Bernoulli((1-cdf(Normal(0,1),mppp[ix,:g]))*λ₁(mppp[ix,:t])/numerator))
		if(mppp[ix,:γ]==0)
			mppp[ix,:Zf]=rand(Truncated(Normal(f₀(mppp[ix,:t]),1), 0, Inf),1)[1]
			mppp[ix,:Zg]=rand(Truncated(Normal(mppp[ix,:g],1), 0, Inf),1)[1]
		else
			mppp[ix,:Zf]=rand(Truncated(Normal(f₁(mppp[ix,:t]),1), 0, Inf),1)[1]
			mppp[ix,:Zg]=rand(Truncated(Normal(mppp[ix,:g],1), -Inf, 0),1)[1]
		end
	end
	#=mppp=vcat(mppp[ixos,:],mppp010,mppp110)=#
	mppp=vcat(mppp[ixos,:],mppp[mppp[:label].=="010",:],mppp[mppp[:label].=="110",:])
	indices=sortperm(mppp[:,:t])
	mppp=mppp[indices,:]
	ixα=find(mppp[:mf].==1)
	nα=length(ixα)
	K=Array(Float64,nα,nα);
	for i=1:nα
		for j=1:nα
			K[i,j]=rho2*exp(-psi2*(mppp[ixα[i],:t]-mppp[ixα[j],:t])^2)
		end
	end
	V=\(K+eye(nα),K)
	L=svd(V)
	L=L[1]*Diagonal(sqrt(L[2]))
	mppp[ixα,:g]=V*mppp[ixα,:Zg]+L*rand(Normal(0,1),nα);
	plot(mppp[:t],cdf(Normal(0,1),mppp[:g]),c="grey",alpha=0.1);
	#=λ=rand(Gamma(nc,1/(T)))=#
end


