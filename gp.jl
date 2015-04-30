using Distributions
using PyPlot

srand(4)

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

L=chol(C+0.0001*eye(nx))'

s2=1
z=rand(Normal(0,1),nx)
g=sqrt(s2)*L*z

lam_star=100
intensity=lam_star*cdf(Normal(0,1),g)
plot(x,intensity)






intensity1(arg)=intensity[maximum(find( y->(y <= arg), x))]

tobs=8
lam_dom=lam_star
n_dom=rand(Poisson(lam_dom*tobs))
o_dom=sort(rand(Uniform(0,tobs),n_dom))
o_acc=Array(Float64,0)
o_thin=Array(Float64,0)
binary=Array(Int64,0)
for i=1:n_dom
	u=rand(Uniform(0,1))
	if(u<=intensity1(o_dom[i])/lam_dom)
		#Accept
		push!(o_acc,o_dom[i])
		push!(binary,1)
	else 
		#Thin
		push!(o_thin,o_dom[i])
		push!(binary,0)
	end
end

PyPlot.plot(o,0.5*ones(length(o)),c="red",marker="|",linestyle="None")
PyPlot.plot(o_thin,-0.5*ones(length(o_thin)),marker="|",linestyle="None",c="green")

C=Array(Float64,n_dom,n_dom)
for i=1:n_dom
	for j=1:n_dom
		C[i,j]=rho2*exp(-psi2*(o_dom[i]-o_dom[j])^2)
	end
end

C=C+0.0001*eye(n_dom)
L=chol(C)'
z=rand(Normal(0,1),n_dom)
g=sqrt(s2)*L*z

M=C*inv(eye(n_dom)+C)
L=chol(M)'


for iter=1:1000

	Z=Array(Float64,n_dom)
	for i=1:n_dom
		if(binary[i]==1)
			Z[i]=rand(Truncated(Normal(g[i],1), 0, Inf),1)[1]
		else
			Z[i]=rand(Truncated(Normal(g[i],1), -Inf, 0),1)[1]
		end
	end

	z=rand(Normal(0,1),n_dom)
	g=M*Z+L*z
	plot(o_dom,lam_star*cdf(Normal(0,1),g),c="red",alpha=0.1)
end
