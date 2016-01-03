f₀(x)=sin(2*(4*x-2))+2*exp(-(16^2)*(x-0.5).^2)
f₁(x)=1-sin(2*(4*x-2))+2*exp(-(3^2)*(x-0.5).^2)-2
Λ=4500
λ₀(x)=Λ*cdf(Normal(0,1),f₀(x))
λ₁(x)=Λ*cdf(Normal(0,1),f₁(x))

srand(4)
n=1500
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
end
g=copy(x)
g1=Dict()
[g1[key]=g[key][1] for key in keys(g)]
function α(time)
	i=0
	while((t[i]<time))
		i=i+1
		if(i==length(t))
			break
		end
	end
	i=i-1
	if(i<0)
		i=0
	end
	return(Φ(g[t[i]][1]))
end
gfun(arg)=g[maximum(find( y->(y <= arg), grid))]
#=plot(x,500*alpha)=#
λs(t)=α(t)*λ₀(t)+(1-α(t))*λ₁(t)
λ₀₁₁(t)=α(t)*λ₀(t)
Ξ₀₁₁=PPProcess(λ₀₁₁,Λ)
Tobs=1
ξ₀₁₁=rand(Ξ₀₁₁,0,Tobs)
λ₁₁₁(t)=(1-α(t))*λ₁(t)
Ξ₁₁₁=PPProcess(λ₁₁₁,Λ)
ξ₁₁₁=rand(Ξ₁₁₁,0,Tobs)
figure()
subplot(211)
grid=collect(0:0.01:1)
plot(grid,λ₀(grid),c="blue",linestyle="--")
plot(grid,λ₁(grid),c="red",linestyle="--")
plot(grid,[λs(t) for t in grid],c="purple")
plot(ξ₀₁₁,-1*ones(length(ξ₀₁₁)),c="blue",marker="|",linestyle="None",markersize=10)
plot(ξ₁₁₁,-1*ones(length(ξ₁₁₁)),c="red",marker="|",linestyle="None",markersize=10)
subplot(212)
plot(grid,[α(t) for t in grid],c="green")
ylim(0,1)


ξ=Dict{Float64,Dict{UTF8String,Float64}}()
for t in ξ₀₁₁
	println(t)
	ξ[t]=Dict("Zg"=>rand(Truncated(Normal(0,1),0,Inf)),"γ"=>0)
end
for t in ξ₁₁₁
	ξ[t]=Dict("Zg"=>rand(Truncated(Normal(0,1),-Inf,0)),"γ"=>1)
end

n=length(ξ)
t=Dict(zip(collect(0:n),[0;sort(collect(keys(ξ)))]))
Δ=Dict(zip(collect(1:n),[t[k]-t[k-1] for k=1:n]))
Q=Dict()
A=Dict()
x=Dict(0.0=>rand(MvNormal(M[0]+0.0000000001(eye(d))),1))
y=Dict()
@time for i=1:n
	A[t[i-1]]=transition(Δ[i],λ)
	Q[t[i-1]]=innovation(Δ[i],λ,q)
	x[t[i]]=A[t[i-1]]*x[t[i-1]]+rand(MvNormal(Q[t[i-1]]+0.0000000001(eye(d))),1)
end
delete!(x,0.0)
g=copy(x)
