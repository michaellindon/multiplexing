Φ(x)=cdf(Normal(0,1),x)

type PPProcess
	λ::Function
	Λ::Any
end

function PPProcess(Λ::Any)
	λ(x)=Λ
	PPProcess(λ,Λ)
end

function PPProcess(Λ::Any,Φ::Function,f::Function)
	λ(x)=Λ*Φ(f(x))
	PPProcess(λ,Λ)
end

function rand(Ξ::PPProcess, a,b)
	nₚ=rand(Poisson(Ξ.Λ*(b-a)))
	Tₚ=sort(rand(Uniform(a,b),nₚ))
	times=Array(Float64,0)
	for t in Tₚ
		if(rand(Bernoulli(Ξ.λ(t)/Ξ.Λ))==1)
			push!(times,t)
		end
	end
	return times
end
