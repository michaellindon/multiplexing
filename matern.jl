using Distributions
using PyPlot
include("statespace.jl")
n=5
ł=0.1
p=3
ν=p+0.5
λ=sqrt(2*ν)/ł
σ²=1.0
q=2*σ²*√π*λ^(2*p+1)*gamma(p+1)/gamma(p+0.5)
d=p+1
L=zeros(d,1)
L[d,1]=1
H=zeros(1,d)
H[1,1]=1
m=Dict(0.0=>zeros(d,1))
M=Dict(0.0=>statcov(λ,q))
t=Dict(zip(collect(0:n),[0;sort(rand(n))]))
Δ=Dict(zip(collect(1:n),[t[k]-t[k-1] for k=1:n]))
Q=Dict()
A=Dict()
x=Dict(0.0=>rand(MvNormal(M[0]+0.0000000001(eye(d))),1))
y=Dict()
@time for i=1:n
	A[t[i-1]]=transition(Δ[i],λ)
	Q[t[i-1]]=innovation(Δ[i],λ,q)
	x[t[i]]=A[t[i-1]]*x[t[i-1]]+rand(MvNormal(Q[t[i-1]]+0.0000000001(eye(d))),1)
	y[t[i]]=H*x[t[i]]+sqrt(σ²)*rand(Normal(0,1))
end
delete!(x,0.0)
plot(sort(collect(keys(x))),[x[i][1] for i in sort(collect(keys(x)))],linestyle="None",marker="o")
plot(sort(collect(keys(y))),[y[i][1] for i in sort(collect(keys(y)))],linestyle="None",marker="o")
plot(sort(collect(keys(x))),[x[i][1] for i in sort(collect(keys(x)))])
plot(sort(collect(keys(y))),[y[i][1] for i in sort(collect(keys(y)))])

@time for i=1:100
	foo=predictFunction(x,rand(Uniform(0,1),1000),λ,q)
	plot(sort(collect(keys(foo))),[foo[i][1] for i in sort(collect(keys(foo)))],c="blue",alpha=0.1)
end

for i=1:100
	foo=predictFunction(x,rand(Uniform(0,1),1000),λ,q)
	plot(sort(collect(keys(foo))),[foo[i][1] for i in sort(collect(keys(foo)))],c="blue",alpha=0.1)
end

#Forward Filtering
delete!(t,0.0)
t=Dict(zip(collect(1:1000),sort([collect(values(t)),rand(995)])))
t[0]=t[1]-(t[2]-t[1])
m=Dict(t[0]=>zeros(d,1))
M=Dict(t[0]=>lyap(F,L*q*L'))
AMAQ=Dict()
@time for i=1:(length(t)-1)
	Δ[i]=t[i]-t[i-1]
	A[t[i-1]]=transition(Δ[i],λ)
	Q[t[i-1]]=innovation(Δ[i],λ,q)
	if(haskey(y,t[i]))
		AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
		#=m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]]*H'*\(σ²+H*AMAQ[t[i-1]]*H',y[t[i]]-H*A[t[i-1]]*m[t[i-1]])=#
		#=M[t[i]]=AMAQ[t[i-1]]-AMAQ[t[i-1]]*H'*\(σ²+H*AMAQ[t[i-1]]*H',H*AMAQ[t[i-1]])=#
		m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
		M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
	else
		AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
		m[t[i]]=A[t[i-1]]*m[t[i-1]]
		M[t[i]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
	end
end

#Backward Sampling
x[t[n]]=m[t[n]]+rand(MvNormal(M[t[n]]+0.00000000001*eye(d)),1)
for i=(n-1):-1:0
	x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+rand(MvNormal(M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],A[t[i]]*M[t[i]])+0.000000001(eye(d))))
end
plot(sort(collect(keys(x))),[x[i][1] for i in sort(collect(keys(x)))])



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



for iiiiii=1:100
	x=FFBS(y,rand(1000),σ²,λ,q)
	plot(sort(collect(keys(x))),[x[i][1] for i in sort(collect(keys(x)))],c="red",alpha=0.1)
end
@time for iiiiii=1:100
	foo=copy(xcopy)
	x=FFBS2(foo,rand(1000),λ,q)
	plot(sort(collect(keys(x))),[x[i][1] for i in sort(collect(keys(x)))],c="red",alpha=0.1)
end

xₒ=Array{Float64}(collect(keys(y)))
output=Array{Float64}([y[key][1] for key in collect(keys(y))])
xₚ=sort(rand(1000))
ρ²=1.0
ψ²=0.5
Kₒₒ=Kernel(xₒ,xₒ,ρ²,ψ²)
Kₚₚ=Kernel(xₚ,xₚ,ρ²,ψ²)
Kₚₒ=Kernel(xₚ,xₒ,ρ²,ψ²)
L=chol(Kₚₚ-Kₚₒ*\(Kₒₒ+0.0000000000001*eye(length(output)),Kₚₒ')+0.00000000000001*eye(length(xₚ)))'
for iiii=1:100
	onehalf=Kₚₒ*\(Kₒₒ+0.00000000001*eye(length(output)),output)+L*rand(Normal(0,1),length(xₚ))
	plot(xₚ,onehalf,c="blue",alpha=0.1)
end
plot(xₒ,output,linestyle="None",marker="o")
