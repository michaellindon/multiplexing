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
C=Array(Float64,nx,nx)
for i=1:nx
	for j=1:nx
		C[i,j]=rho2*exp(-psi2*(x[i]-x[j])^2)
	end
end

#=L=chol(C+0.0001*eye(nx))'=#
Csvd=svd(C)
L=Csvd[1]*Diagonal(sqrt(Csvd[2]))
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
no=rand(Poisson(λ*T))
o=sort(rand(Uniform(0,T),no))
oacc=Array(Float64,0)
othin=Array(Float64,0)
binary=Array(Int64,0)
for i=1:no
	u=rand(Uniform(0,1))
	if(u<=intensity1(o[i])/λ)
		#Accept
		push!(oacc,o[i])
		push!(binary,1)
	else 
		#Thin
		push!(othin,o[i])
		push!(binary,0)
	end
end
na=length(oacc)

plot(oacc,0.5*ones(length(oacc)),c="blue",marker="|",linestyle="None",markersize=30)
#=plot(othin,-0.5*ones(length(othin)),marker="|",linestyle="None",c="green")=#

C=Array(Float64,no,no)
for i=1:no
	for j=1:no
		C[i,j]=rho2*exp(-psi2*(o[i]-o[j])^2)
	end
end

C=C+0.0001*eye(no)
L=chol(C)'
z=rand(Normal(0,1),no)
g=sqrt(s2)*L*z

M=C*inv(eye(no)+C)
L=chol(M)'


for iter=1:10
	Z=Array(Float64,no);
	for i=1:no
		if(binary[i]==1)
			Z[i]=rand(Truncated(Normal(g[i],1), 0, Inf),1)[1]
		else
			Z[i]=rand(Truncated(Normal(g[i],1), -Inf, 0),1)[1]
		end
	end
	C=Array(Float64,no,no);
	for i=1:no
		for j=1:no
			C[i,j]=rho2*exp(-psi2*(o[i]-o[j])^2)
		end
	end
	#=C=C+0.0001*eye(no);=#
	M=C*inv(eye(no)+C);
	Msvd=svd(M)
	L=Msvd[1]*Diagonal(sqrt(Msvd[2]))
	#=L=chol(M)';=#
	g=M*Z+L*rand(Normal(0,1),no);
	plot(o,lam_star*cdf(Normal(0,1),g),c="grey",alpha=0.1);
	n_prop=rand(Poisson(λ*T));
	oprop=sort(rand(Uniform(0,T),n_prop));
	Cpp=Array(Float64,n_prop,n_prop);
	for i=1:n_prop
		for j=1:n_prop
			Cpp[i,j]=rho2*exp(-psi2*(oprop[i]-oprop[j])^2)
		end
	end
	Cpo=Array(Float64,n_prop,no);
	for i=1:n_prop
		for j=1:no
			Cpo[i,j]=rho2*exp(-psi2*(oprop[i]-o[j])^2)
		end
	end
	gprop=Array(Float64,n_prop);
	gprop=Cpo*\(C,g)+chol(Cpp-Cpo*\(C,Cpo')+0.0001*eye(n_prop))'*rand(Normal(0,1),n_prop);
	gthin=Array(Float64,0);
	othin=Array(Float64,0);
	for i=1:n_prop
		u=rand(Uniform(0,1))
		if(u<=lam_star*cdf(Normal(0,1),gprop[i])/λ)
			#Accept
		else 
			#Thin
			push!(othin,oprop[i])
			push!(gthin,gprop[i])
		end
	end
	nt=length(othin);
	no=nt+na;
	o=Array(Float64,0);
	gacc=g[find(x->x==1,binary)];
	g=Array(Float64,0);
	binary=Array(Int64,0);
	ai=1;
	ti=1;
	while ti<=nt && ai<=na
		if(othin[ti]<oacc[ai])
			push!(o,othin[ti])
			push!(binary,0)
			push!(g,gthin[ti])
			ti=ti+1
		else
			push!(o,oacc[ai])
			push!(binary,1)
			push!(g,gacc[ai])
			ai=ai+1
		end
	end
	while ti<=nt
		push!(o,othin[ti])
		push!(g,gthin[ti])
		push!(binary,0)
		ti=ti+1
	end
	while ai<=na
		push!(o,oacc[ai])
		push!(binary,1)
		push!(g,gacc[ai])
		ai=ai+1
	end
end
