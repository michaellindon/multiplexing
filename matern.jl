using Distributions
using PyPlot
include("transition.jl")
n=1000
ł=0.1
p=3
ν=p+0.5
λ=sqrt(2*ν)/ł
σ²=1.0
q=2*σ²*√π*λ^(2*p+1)*gamma(p+1)/gamma(p+0.5)
bline=reverse([λ^i for i=1:(p+1)])
for i=1:(p+1)
	bline[i]=bline[i]*binomial(p+1,i-1)
end
F=vcat(hcat(zeros(p),eye(p)),-bline')
F=convert(Array{Float64,2},F)
d=p+1
L=zeros(d,1)
L[d,1]=1
H=zeros(1,d)
H[1,1]=1
# sylvester(A, B, C) Computes the solution X to the Sylvester equation AX + XB + C = 0
#=M=Dict(0=>sylvester(F,F',L*q*L'))=#
# lyap(A, C) Computes the solution X to the continuous Lyapunov equation AX + XA' + C = 0
m=Dict(0=>zeros(d,1))
M=Dict(0=>lyap(F,L*q*L'))
t=Dict(zip(collect(0:n),[0;sort(rand(n))]))
Δ=Dict(zip(collect(1:n),[t[k]-t[k-1] for k=1:n]))
Q=Dict()
A=Dict()
AMAQ=Dict()
x=Dict(0=>rand(MvNormal(M[0]+0.0000000001(eye(d))),1))
y=Dict()
@time for i=1:n
	#=A=expm(F*Δ[i])=#
	H=expm(Δ[i]*vcat(hcat(-F,L*q*L'),hcat(zeros(d,d),F')))
	A[i-1]=H[d+1:end,d+1:end]'
	Q[i-1]=A[i-1]*H[1:d,d+1:end]
	#=x[i]=A[i-1]*x[i-1]+rand(MvNormal(Q[i-1]+0.0000000001(eye(d))))=#
end
@time for i=1:n
	A[i-1]=transition(Δ[i],λ)
	Q[i-1]=innovation(Δ[i],λ,q)
	x[i]=A[i-1]*x[i-1]+rand(MvNormal(Q[i-1]+0.0000000001(eye(d))),1)
	y[i]=H*x[i]+sqrt(σ²)*rand(Normal(0,1))
end
plot([x[i][1] for i in sort(collect(keys(x)))],linestyle="None",marker="o")
plot([y[i][1] for i in sort(collect(keys(y)))],linestyle="None",marker="o")
plot([x[i][1] for i in sort(collect(keys(x)))])
plot([y[i][1] for i in sort(collect(keys(y)))])


#Forward Filtering
@time for i=1:n
	AMAQ[i-1]=A[i-1]*M[i-1]*A[i-1]'+Q[i-1]
	#=m[i]=A[i-1]*m[i-1]+AMAQ[i-1]*H'*\(σ²+H*AMAQ[i-1]*H',y[i]-H*A[i-1]*m[i-1])=#
	#=M[i]=AMAQ[i-1]-AMAQ[i-1]*H'*\(σ²+H*AMAQ[i-1]*H',H*AMAQ[i-1])=#
	m[i]=A[i-1]*m[i-1]+AMAQ[i-1][:,1]*(y[i]-A[i-1][1,:]*m[i-1])/(σ²+AMAQ[i-1][1,1])
	M[i]=AMAQ[i-1]-(AMAQ[i-1][:,1]*AMAQ[i-1][1,:])/(σ²+AMAQ[i-1][1,1])
end

#Backward Sampling
x[n]=m[n]+rand(MvNormal(M[n]+0.00000000001*eye(d)),1)
@time for i in reverse(0:n-1)
	x[i]=m[i]+M[i]*A[i]'*\(AMAQ[i],x[i+1]-A[i]*m[i])+rand(MvNormal(M[i]-M[i]*A[i]'*\(AMAQ[i],A[i]*M[i])+0.000000001(eye(d))))
end
plot([x[i][1] for i in sort(collect(keys(x)))])



V=zeros(n*d,n*d)
V[1:d,1:d]=A[0]*M[0]*A[0]'+Q[0]
for i=2:n
	V[(i-1)*d+1:i*d,(i-1)*d+1:i*d]=A[i-1]*V[(i-1-1)*d+1:(i-1)*d,(i-1-1)*d+1:(i-1)*d]*A[i-1]'+Q[i-1]
end
for i=2:n
	for j=1:(i-1)
		V[(i-1)*d+1:i*d,(j-1)*d+1:j*d]=V[(j-1)*d+1:j*d,(j-1)*d+1:j*d]
		for k=j:(i-1)
			V[(i-1)*d+1:i*d,(j-1)*d+1:j*d]=A[k]*V[(i-1)*d+1:i*d,(j-1)*d+1:j*d]
		end
		V[(j-1)*d+1:j*d,(i-1)*d+1:i*d]=V[(i-1)*d+1:i*d,(j-1)*d+1:j*d]'
	end
end

H=zeros(n*d,1)
for i=1:400
	if((i-1)%4==0)
		H[i,1]=1
	end
end
V[H.==1,H.==1]

matshow(inv(V[H.==1,H.==1]))


g=chol(Kernel(sort(collect(values(t))),sort(collect(values(t))),1.0,1.0)+0.000000001*eye(length(sort(collect(values(t))))))'*rand(Normal(0,1),length(sort(collect(values(t)))))
plot(g)
