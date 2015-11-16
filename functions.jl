Φ(x)=cdf(Normal(0,1),x)

function realizationGP(x,σ²,ρ²,ψ²)
	ϵ=0.0000001
	n=size(x,1)
	K=Kernel(x,x,ρ²,ψ²)
	return rand(MvNormal(σ²*K+ϵ*eye(n)))
end

function logDensityGP(y,x,σ²,ρ²,ψ²)
	n=size(x,1)
	Iₒₒ=eye(n)
	K=Kernel(x,x,ρ²,ψ²)
	return 	-0.5*n*log(σ²)-0.5*logdet(K+Iₒₒ)-0.5*(1/σ²)*dot(y,\(K+Iₒₒ,y))
end
function logDensityGP(xy::Dict,σ²,ρ²,ψ²)
	x=collect(keys(xy))
	y=collect(values(xy))
	return logDensityGP(y,x,σ²,ρ²,ψ²)
end

function Kernel(x::Vector{Float64},ρ²::Float64,ψ²::Float64)
	K=Array(Float64,length(x),length(x))
	for r=1:length(x)
		for c=1:length(x)
			K[r,c]=ρ²*exp(-ψ²*(x[r]-x[c])*(x[r]-x[c]))
		end
	end
	ϵ=0.00000000001
	K=PDMat(K+ScalMat(length(x),ϵ))
	return(K)
end
function Kernel(x₁::Vector{Float64},x₂::Vector{Float64},ρ²,ψ²)
	K=Array(Float64,length(x₁),length(x₂))
	for r=1:length(x₁)
		for c=1:length(x₂)
			K[r,c]=ρ²*exp(-ψ²*(x₁[r]-x₂[c])*(x₁[r]-x₂[c]))
		end
	end
	return(K)
end

function PredictGP(xₚ,zₒ,xₒ,σ²,ρ²,ψ²,inputType)
	ϵ=0.0000001
	Kₒₒ=Kernel(xₒ,xₒ,ρ²,ψ²)
	Kₚₚ=Kernel(xₚ,xₚ,ρ²,ψ²)
	Kₚₒ=Kernel(xₚ,xₒ,ρ²,ψ²)
	nₚ=size(xₚ,1)
	Iₚₚ=eye(nₚ)
	nₒ=size(xₒ,1)
	Iₒₒ=eye(nₒ)
	if(inputType=="response")
		yₒ=zₒ
		μ=Kₚₒ*\(Kₒₒ+Iₒₒ,yₒ)
		Σ=Kₚₚ-Kₚₒ*\(Kₒₒ+Iₒₒ,Kₚₒ')
	elseif(inputType=="function")
		fₒ=zₒ
		μ=Kₚₒ*\(Kₒₒ+ϵ*Iₒₒ,fₒ)
		Σ=Kₚₚ-Kₚₒ*\(Kₒₒ+ϵ*Iₒₒ,Kₚₒ')
	else
		println("Error: Type of input must be either response or function")
	end
	fₚ=rand(MvNormal(μ,σ²*Σ+ϵ*Iₚₚ),1)[:,1]
	#=yₚ=fₚ+rand(Normal(0,sqrt(σ²)),nₚ)=#
	return Dict(zip(xₚ,fₚ))
end	
function PredictGP(xₚ,xₒzₒ::Dict,σ²,ρ²,ψ²,inputType)
	xₒ=collect(keys(xₒzₒ))
	zₒ=collect(values(xₒzₒ))
	return PredictGP(xₚ,zₒ,xₒ,σ²,ρ²,ψ²,inputType)
end	
function PosteriorGP(yₒ,xₒ,σ²,ρ²,ψ²)
	ϵ=0.00000001
	Kₒₒ=Kernel(xₒ,xₒ,ρ²,ψ²)
	nₒ=size(xₒ,1)
	Iₒₒ=eye(nₒ)
	μ=Kₒₒ*\(Kₒₒ+Iₒₒ,yₒ)
	Σ=Kₒₒ-Kₒₒ*\(Kₒₒ+Iₒₒ,Kₒₒ')
	fₒ=rand(MvNormal(μ,σ²*Σ+ϵ*Iₒₒ),1)[:,1]
	return Dict(zip(xₒ,fₒ))
end	

function PosteriorGP(xₒyₒ::Dict,σ²,ρ²,ψ²)
	xₒ=collect(keys(xₒyₒ))
	yₒ=collect(values(xₒyₒ))
	return PosteriorGP(yₒ,xₒ,σ²,ρ²,ψ²)
end	

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

function realization(Ξ::PPProcess, a,b)
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
