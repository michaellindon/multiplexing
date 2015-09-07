using Distributions
using PyPlot

#=GENERATE DATA=#
srand(7)
srand(10)
srand(11)

x= [0:0.1:8]
nx=length(x)
rho2=1
psi2=0.1
K=Array(Float64,nx,nx)
for i=1:nx
	for j=1:nx
		K[i,j]=rho2*exp(-psi2*(x[i]-x[j])^2)
	end
end

#=L=chol(K+0.0001*eye(nx))'=#
Ksvd=svd(K)
L=Ksvd[1]*Diagonal(sqrt(Ksvd[2]))
s2=1
z=rand(Normal(0,1),nx)
g=sqrt(s2)*L*z
plot(x,g)
lam_star=100
intensity=lam_star*cdf(Normal(0,1),g)
plot(x,intensity)
intensity1(arg)=intensity[maximum(find( y->(y <= arg), x))]

T=8
λ=lam_star
nc=rand(Poisson(λ*T))
tc=sort(rand(Uniform(0,T),nc))
m=Array(Bool,nc)
for i=1:nc
	u=rand(Uniform(0,1))
	if(u<=intensity1(tc[i])/λ)
		#Accept
		m[i]=true
	else 
		#Thin
		m[i]=false
	end
end
to=tc[m]
tu=tc[!m]
no=length(to)
nu=length(tu)

plot(to,0.5*ones(no),c="blue",marker="|",linestyle="None",markersize=30)

K=Array(Float64,nc,nc)
for i=1:nc
	for j=1:nc
		K[i,j]=rho2*exp(-psi2*(tc[i]-tc[j])^2)
	end
end

K=K+0.0001*eye(nc)
L=chol(K)'
z=rand(Normal(0,1),nc)
g=sqrt(s2)*L*z

#=Prior Specification=#
ν=1
λ₀=1

for iter=1:10
	Z=Array(Float64,nc);
	for i=1:nc
		if(m[i]==1)
			Z[i]=rand(Truncated(Normal(g[i],1), 0, Inf),1)[1]
		else
			Z[i]=rand(Truncated(Normal(g[i],1), -Inf, 0),1)[1]
		end
	end
	K=Array(Float64,nc,nc);
	for i=1:nc
		for j=1:nc
			K[i,j]=rho2*exp(-psi2*(tc[i]-tc[j])^2)
		end
	end
	#=M=K*inv(eye(nc)+K)=#
	#=Msvd=svd(M)=#
	#=L=Msvd[1]*Diagonal(sqrt(Msvd[2]))=#
	V=\(K+eye(nc),K)
	L=svd(V)
	L=L[1]*Diagonal(sqrt(L[2]))
	g=V*Z+L*rand(Normal(0,1),nc);
	plot(tc,λ*cdf(Normal(0,1),g),c="grey",alpha=0.1);
	np=rand(Poisson(λ*T));
	tp=sort(rand(Uniform(0,T),np));
	Kpp=Array(Float64,np,np);
	for i=1:np
		for j=1:np
			Kpp[i,j]=rho2*exp(-psi2*(tp[i]-tp[j])^2)
		end
	end
	Kpc=Array(Float64,np,nc);
	for i=1:np
		for j=1:nc
			Kpc[i,j]=rho2*exp(-psi2*(tp[i]-tc[j])^2)
		end
	end
	gp=Array(Float64,np);
	#gp=Kpc*\(K,g)+chol(Kpp-Kpc*\(K,Kpc')+0.0001*eye(np))'*rand(Normal(0,1),np);
	L=svd(Kpp-Kpc*\(K,Kpc'))
	L=L[1]*Diagonal(sqrt(L[2]))
	gp=Kpc*\(K,g)+L*rand(Normal(0,1),np);
	gu=Array(Float64,0);
	tu=Array(Float64,0);
	for i=1:np
		u=rand(Uniform(0,1))
		if(u<=cdf(Normal(0,1),gp[i]))
			#Accept
		else 
			#Thin
			push!(tu,tp[i])
			push!(gu,gp[i])
		end
	end
	nu=length(tu);
	nc=no+nu;
	go=g[find(m)];
	tc=[to,tu];
	g=[go,gu];
	m=[ones(no),zeros(nu)]
	indices=sortperm([to,tu])
	tc=tc[indices]
	g=g[indices]
	m=m[indices]
	λ=rand(Gamma(nc,1/(T)))
end
