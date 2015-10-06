f₀(x)=sin(2*(4*x-2))+2*exp(-(16^2)*(x-0.5).^2)+2
f₁(x)=-sin(2*(4*x-2))+2*exp(-(3^2)*(x-0.5).^2)-2
Λ=1000
λ₀(x)=Λ*cdf(Normal(0,1),f₀(x))
λ₁(x)=Λ*cdf(Normal(0,1),f₁(x))

srand(2)
disc=0.01
Tobs=1
grid=collect(0:disc:Tobs)
σ²=1
ρ²=1
ψ²=40.0*ones(1)
g=realizationGP(grid,σ²,ρ²,ψ²)
gfun(arg)=g[maximum(find( y->(y <= arg), grid))]
alpha=Φ(g)
#=plot(x,500*alpha)=#
α(arg::Float64)=alpha[maximum(find( y->(y <= arg), grid))]
λs(t)=α(t)*λ₀(t)+(1-α(t))*λ₁(t)
λ₀₁₁(t)=α(t)*λ₀(t)
Ξ₀₁₁=PPProcess(λ₀₁₁,Λ)
ξ₀₁₁=realization(Ξ₀₁₁,0,Tobs)
λ₁₁₁(t)=(1-α(t))*λ₁(t)
Ξ₁₁₁=PPProcess(λ₁₁₁,Λ)
ξ₁₁₁=realization(Ξ₁₁₁,0,Tobs)
figure()
subplot(211)
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
	ξ[t]=Dict("Zg"=>rand(Truncated(Normal(gfun(t),1),0,Inf)),"γ"=>0)
end
for t in ξ₁₁₁
	ξ[t]=Dict("Zg"=>rand(Truncated(Normal(gfun(t),1),-Inf,0)),"γ"=>1)
end
g=Dict{Float64,Float64}(zip(collect(keys(ξ)),[gfun(t)::Float64 for t in collect(keys(ξ))]))
