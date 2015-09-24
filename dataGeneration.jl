f₀(x)=sin(2*(4*x-2))+2*exp(-(16^2)*(x-0.5).^2)
f₁(x)=cos(2*(4*x-2))+2*exp(-(3^2)*(x-0.5).^2)
λc=1000
λ₀(x)=λc*cdf(Normal(0,1),f₀(x))
λ₁(x)=λc*cdf(Normal(0,1),f₁(x))
δ=[0:0.01:1]
#=plot(δ,f₁(δ),c="blue")=#
#=plot(δ,f₂(δ),c="blue")=#
#=plot(δ,λ₁(δ),c="blue")=#
#=plot(δ,λ₂(δ),c="blue")=#
#=GENERATE DATA=#

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
mppp0=Array(Float64,nc,7)
#Times#Component#λ-thinned#α-thinned#λ-Z#α-Z
mppp0[:,1]=tc
mppp0[:,2]=zeros(nc)
for i=1:nc
	mppp0[i,5]=rand(Normal(f₀(tc[i]),1))
	if(mppp0[i,5]>0)
		#Accept
		mppp0[i,3]=1.0
		mppp0[i,6]=rand(Normal(gfun(tc[i]),1))
		if(mppp0[i,6]>0)
			mppp0[i,4]=1.0
		else
			mppp0[i,4]=0.0
		end
	else 
		#Thin
		mppp0[i,3]=0.0
		mppp0[i,4]=0.0
	end
end
nc=rand(Poisson(λc*T))
tc=sort(rand(Uniform(0,T),nc))
mppp1=Array(Float64,nc,7)
#Times#Component#λ-thinned#α-thinned#λ-Z#α-Z
mppp1[:,1]=tc
mppp1[:,2]=ones(nc)
for i=1:nc
	mppp1[i,5]=rand(Normal(f₁(tc[i]),1))
	if(mppp1[i,5]>0)
		#Accept
		mppp1[i,3]=1.0
		mppp1[i,6]=rand(Normal(gfun(tc[i]),1))
		if(mppp1[i,6]<0)
			mppp1[i,4]=1.0
		else
			mppp1[i,4]=0.0
		end
	else 
		#Thin
		mppp1[i,3]=0.0
		mppp1[i,4]=0.0
	end
end


mppp=vcat(mppp0,mppp1)
indices=sortperm(mppp[:,1])
mppp=mppp[indices,:]
ix0=find((mppp[:,2].==0.0)&(mppp[:,4].==1.0))
ix1=find((mppp[:,2].==1.0)&(mppp[:,4].==1.0))
plot(mppp[ix0,1],-35*ones(length(ix0)),c="blue",marker="|",linestyle="None",markersize=10)
plot(mppp[ix1,1],-35*ones(length(ix1)),c="red",marker="|",linestyle="None",markersize=10)
subplot(212)
plot(x,alpha,c="green")
