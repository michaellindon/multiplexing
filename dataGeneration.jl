f₀(x)=sin(2*(4*x-2))+2*exp(-(16^2)*(x-0.5).^2)
f₁(x)=cos(2*(4*x-2))+2*exp(-(3^2)*(x-0.5).^2)
λc=2000
λ₀(x)=λc*cdf(Normal(0,1),f₀(x))
λ₁(x)=λc*cdf(Normal(0,1),f₁(x))
δ=[0:0.01:1]

srand(1)
x=[0:0.01:1]
nx=length(x)
rho2=1
psi2=80
K=Array(Float64,nx,nx)
for i=1:nx
	for j=1:nx
		K[i,j]=rho2*exp(-psi2*(x[i]-x[j])^2)
	end
end
Ksvd=svd(K)
L=Ksvd[1]*Diagonal(sqrt(Ksvd[2]))
s2=1
z=rand(Normal(0,1),nx)
g=sqrt(s2)*L*z
gfun(arg)=g[maximum(find( y->(y <= arg), x))]
alpha=cdf(Normal(0,1),g)
#=plot(x,500*alpha)=#
alphafun(arg)=alpha[maximum(find( y->(y <= arg), x))]
δ=[0:0.01:1]
λs=Array(Float64,length(δ))
for t=1:length(δ)
	λs[t]=alphafun(δ[t])*λ₀(δ[t])+(1-alphafun(δ[t]))*λ₁(δ[t])
end
figure()
subplot(211)
plot(δ,λ₀(δ),c="blue",linestyle="--")
plot(δ,λ₁(δ),c="red",linestyle="--")
plot(δ,λs,c="purple")


T=1
nc=rand(Poisson(λc*T))
tc=sort(rand(Uniform(0,T),nc))
#Times#Component#λ-thinned#α-thinned#λ-Z#α-Z
mppp0=DataFrame(t=tc,γ=zeros(Int64,nc),mf=zeros(Int64,nc),mg=Array(Int64,nc),Zf=Array(Float64,nc),Zg=Array(Float64,nc),label=Array(String,nc),g=Array(Float64,nc))
for i=1:nc
	mppp0[i,:g]=gfun(tc[i])
	mppp0[i,:Zf]=rand(Normal(f₀(tc[i]),1))
	if(mppp0[i,:Zf]>0)
		#Accept
		mppp0[i,:mf]=1
		mppp0[i,:Zg]=rand(Normal(gfun(tc[i]),1))
		if(mppp0[i,:Zg]>0)
			mppp0[i,:mg]=1
			mppp0[i,:label]=string(mppp0[i,:γ],mppp0[i,:mf],mppp0[i,:mg])
		else
			mppp0[i,:mg]=0
			mppp0[i,:label]=string(mppp0[i,:γ],mppp0[i,:mf],mppp0[i,:mg])
		end
	else 
		#Thin
		mppp0[i,:mf]=0
		mppp0[i,:mg]=NA
		mppp0[i,:Zg]=NA
		mppp0[i,:label]=string(mppp0[i,:γ],mppp0[i,:mf],mppp0[i,:mg])
	end
end
###Sanity Checks###
if(!all(mppp0[mppp0[:label].=="00NA",:Zf].<0))
	println("Error")
end
if(!all(mppp0[mppp0[:label].=="011",:Zf].>0))
	println("Error")
end
if(!all(mppp0[mppp0[:label].=="011",:Zg].>0))
	println("Error")
end
if(!all(mppp0[mppp0[:label].=="010",:Zf].>0))
	println("Error")
end
if(!all(mppp0[mppp0[:label].=="010",:Zg].<0))
	println("Error")
end
##################

nc=rand(Poisson(λc*T))
tc=sort(rand(Uniform(0,T),nc))
#Times#Component#λ-thinned#α-thinned#λ-Z#α-Z
mppp1=DataFrame(t=tc,γ=ones(Int64,nc),mf=Array(Int64,nc),mg=Array(Int64,nc),Zf=Array(Float64,nc),Zg=Array(Float64,nc),label=Array(String,nc),g=Array(Float64,nc))
for i=1:nc
	mppp1[i,:g]=gfun(tc[i])
	mppp1[i,:Zf]=rand(Normal(f₁(tc[i]),1))
	if(mppp1[i,:Zf]>0)
		#Accept
		mppp1[i,:mf]=1
		mppp1[i,:Zg]=rand(Normal(gfun(tc[i]),1))
		if(mppp1[i,:Zg]<0)
			mppp1[i,:mg]=1
			mppp1[i,:label]=string(mppp1[i,:γ],mppp1[i,:mf],mppp1[i,:mg])
		else
			mppp1[i,:mg]=0
			mppp1[i,:label]=string(mppp1[i,:γ],mppp1[i,:mf],mppp1[i,:mg])
		end
	else 
		#Thin
		mppp1[i,:mf]=0
		mppp1[i,:mg]=NA
		mppp1[i,:Zg]=NA
		mppp1[i,:label]=string(mppp1[i,:γ],mppp1[i,:mf],mppp1[i,:mg])
	end
end
###Sanity Checks###
if(!all(mppp1[mppp1[:label].=="10NA",:Zf].<0))
	println("Error")
end
if(!all(mppp1[mppp1[:label].=="111",:Zf].>0))
	println("Error")
end
if(!all(mppp1[mppp1[:label].=="111",:Zg].<0))
	println("Error")
end
if(!all(mppp1[mppp1[:label].=="110",:Zf].>0))
	println("Error")
end
if(!all(mppp1[mppp1[:label].=="110",:Zg].>0))
	println("Error")
end
##################


mppp=vcat(mppp0,mppp1)
indices=sortperm(mppp[:t])
mppp=mppp[indices,:]
ix0=find((mppp[:,:γ].==0.0)&(mppp[:,:mg].==1.0))
ix1=find((mppp[:,:γ].==1.0)&(mppp[:,:mg].==1.0))
plot(mppp[ix0,1],-35*ones(length(ix0)),c="blue",marker="|",linestyle="None",markersize=10)
plot(mppp[ix1,1],-35*ones(length(ix1)),c="red",marker="|",linestyle="None",markersize=10)
subplot(212)
plot(x,alpha,c="green")
ρ²=rho2
ψ²=Array(Float64,1)
ψ²[1]=psi2
