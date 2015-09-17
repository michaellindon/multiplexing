using Distributions
using PyPlot

f₀(x)=sin(2*(4*x-2))+2*exp(-(16^2)*(x-0.5).^2)-1
f₁(x)=cos(2*(4*x-2))+2*exp(-(3^2)*(x-0.5).^2)-1
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
mppp1=Array(Float64,nc,6)
#Times#Component#λ-thinned#α-thinned#λ-Z#α-Z
mppp1[:,1]=tc
mppp1[:,2]=zeros(nc)
for i=1:nc
	mppp1[i,5]=rand(Normal(f₀(tc[i]),1))
	if(mppp1[i,5]>0)
		#Accept
		mppp1[i,3]=1.0
		mppp1[i,6]=rand(Normal(gfun(tc[i]),1))
		if(mppp1[i,6]>0)
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
nc=rand(Poisson(λc*T))
tc=sort(rand(Uniform(0,T),nc))
mppp2=Array(Float64,nc,6)
#Times#Component#λ-thinned#α-thinned
mppp2[:,1]=tc
mppp2[:,2]=ones(nc)
for i=1:nc
	mppp2[i,5]=rand(Normal(f₁(tc[i]),1))
	if(mppp2[i,5]>0)
		#Accept
		mppp2[i,3]=1.0
		mppp2[i,6]=rand(Normal(gfun(tc[i]),1))
		if(mppp2[i,6]<0)
			mppp2[i,4]=1.0
		else
			mppp2[i,4]=0.0
		end
	else 
		#Thin
		mppp2[i,3]=0.0
		mppp2[i,4]=0.0
	end
end


mppp=vcat(mppp1,mppp2)
indices=sortperm(mppp[:,1])
mppp=mppp[indices,:]
ix0=find((mppp[:,2].==0.0)&(mppp[:,4].==1.0))
ix1=find((mppp[:,2].==1.0)&(mppp[:,4].==1.0))
plot(mppp[ix0,1],-35*ones(length(ix0)),c="blue",marker="|",linestyle="None",markersize=10)
plot(mppp[ix1,1],-35*ones(length(ix1)),c="red",marker="|",linestyle="None",markersize=10)
subplot(212)
plot(x,alpha,c="green")

#=Prior Specification=#
ν=1
λ₀=1

γ=rand(Bernoulli(0.5),nc)
nc=size(mppp)[1]
g=zeros(Float64,nc)
for iter=1:200
	Z=Array(Float64,nc);
	for i=1:nc  #m switches depending on the function it came from!
		if(mppp[i,2]==0 && mppp[i,4]==1)
			Z[i]=rand(Truncated(Normal(g[i],1), 0, Inf),1)[1]
		elseif(mppp[i,2]==0 && mppp[i,4]==0)
			Z[i]=rand(Truncated(Normal(g[i],1), -Inf, 0),1)[1]
		elseif(mppp[i,2]==1 && mppp[i,4]==1)
			Z[i]=rand(Truncated(Normal(g[i],1), -Inf, 0),1)[1]
		else
			Z[i]=rand(Truncated(Normal(g[i],1), 0, Inf),1)[1]
		end
	end
	K=Array(Float64,nc,nc);
	for i=1:nc
		for j=1:nc
			K[i,j]=rho2*exp(-psi2*(mppp[i,1]-mppp[j,1])^2)
		end
	end
	#=M=K*inv(eye(nc)+K)=#
	#=Msvd=svd(M)=#
	#=L=Msvd[1]*Diagonal(sqrt(Msvd[2]))=#
	V=\(K+eye(nc),K)
	L=svd(V)
	L=L[1]*Diagonal(sqrt(L[2]))
	g=V*Z+L*rand(Normal(0,1),nc);
	plot(mppp[:,1],cdf(Normal(0,1),g),c="grey",alpha=0.1);
	#=np1=rand(Poisson(λc*T));=#
	#=tp1=sort(rand(Uniform(0,T),np1));=#
	#=Kpp1=Array(Float64,np1,np1);=#
	#=for i=1:np1=#
		#=for j=1:np1=#
			#=Kpp1[i,j]=rho2*exp(-psi2*(tp1[i]-tp1[j])^2)=#
		#=end=#
	#=end=#
	#=Kpc1=Array(Float64,np1,nc);=#
	#=for i=1:np1=#
		#=for j=1:nc=#
			#=Kpc1[i,j]=rho2*exp(-psi2*(tp1[i]-tc[j])^2)=#
		#=end=#
	#=end=#
	#=αp1=Array(Float64,np1);=#
	#=#gp=Kpc*\(K,g)+chol(Kpp-Kpc*\(K,Kpc')+0.0001*eye(np))'*rand(Normal(0,1),np);=#
	#=L=svd(Kpp1-Kpc1*\(K,Kpc1'))=#
	#=L=L[1]*Diagonal(sqrt(L[2]))=#
	#=αp1=Kpc1*\(K,α)+L*rand(Normal(0,1),np1);=#
	#=αu1=Array(Float64,0);=#
	#=tu1=Array(Float64,0);=#
	#=γu1=Array(Int64,0)=#
	#=mu1=Array(Int64,0)=#
	#=for i=1:np1=#
		#=u=rand(Uniform(0,1))=#
		#=if(u<=cdf(Normal(0,1),αp1[i]))=#
			#=#Accept=#
		#=else =#
			#=#Thin=#
			#=push!(tu1,tp1[i])=#
			#=push!(αu1,αp1[i])=#
			#=push!(γu1,1)=#
			#=push!(mu1,0)=#
		#=end=#
	#=end=#
	#=nu1=length(tu1);=#
	#=np2=rand(Poisson(λc*T));=#
	#=tp2=sort(rand(Uniform(0,T),np2));=#
	#=Kpp2=Array(Float64,np2,np2);=#
	#=for i=1:np2=#
		#=for j=1:np2=#
			#=Kpp2[i,j]=rho2*exp(-psi2*(tp2[i]-tp2[j])^2)=#
		#=end=#
	#=end=#
	#=Kpc2=Array(Float64,np2,nc);=#
	#=for i=1:np2=#
		#=for j=1:nc=#
			#=Kpc2[i,j]=rho2*exp(-psi2*(tp2[i]-tc[j])^2)=#
		#=end=#
	#=end=#
	#=αp2=Array(Float64,np2);=#
	#=#gp=Kpc*\(K,g)+chol(Kpp-Kpc*\(K,Kpc')+0.0001*eye(np))'*rand(Normal(0,1),np);=#
	#=L=svd(Kpp2-Kpc2*\(K,Kpc2'))=#
	#=L=L[1]*Diagonal(sqrt(L[2]))=#
	#=αp2=Kpc2*\(K,α)+L*rand(Normal(0,1),np2);=#
	#=αu2=Array(Float64,0);=#
	#=tu2=Array(Float64,0);=#
	#=γu2=Array(Int64,0)=#
	#=mu2=Array(Int64,0)=#
	#=for i=1:np2=#
		#=u=rand(Uniform(0,1))=#
		#=if(u<=cdf(Normal(0,1),αp2[i]))=#
			#=#Accept=#
			#=push!(tu2,tp2[i])=#
			#=push!(αu2,αp2[i])=#
			#=push!(γu2,2)=#
			#=push!(mu2,1)=#
		#=else =#
			#=#Thin=#
		#=end=#
	#=end=#
	#=nu2=length(tu2);=#

	#=γo=Array(Int64,no)=#
	#=mo=Array(Int64,no)=#
	#=αo=α[find(observed)]=#
	#=for i=1:no=#
		#=numerator=cdf(Normal(0,1),αo[i])*λ₁(to[i])+(1-cdf(Normal(0,1),αo[i]))*λ₂(to[i])=#
		#=γo[i]=rand(Bernoulli((1-cdf(Normal(0,1),αo[i]))*λ₂(to[i])/numerator))+1=#
		#=if(γo[i]==1)=#
			#=mo[i]=1=#
		#=else=#
			#=mo[i]=0=#
		#=end=#
	#=end=#
	#=observed=[ones(no),zeros(nu1),zeros(nu2)]=#
	#=nc=no+nu1+nu2;=#
	#=tc=[to,tu1,tu2];=#
	#=α=[αo,αu1,αu2];=#
	#=m=[mo,mu1,mu2]=#
	#=γ=[γo,γu1,γu2]=#
	#=indices=sortperm(tc)=#
	#=observed=observed[indices]=#
	#=tc=tc[indices]=#
	#=α=α[indices]=#
	#=m=m[indices]=#
	#=γ=γ[indices]=#
	#=λ=rand(Gamma(nc,1/(T)))=#
end
