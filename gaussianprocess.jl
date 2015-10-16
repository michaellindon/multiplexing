function Kernel(x₁::Vector{Float64},x₂::Vector{Float64},ρ²,ψ²)
	K=Array(Float64,length(x₁),length(x₂))
	for r=1:length(x₁)
		for c=1:length(x₂)
			K[r,c]=ρ²*exp(-ψ²*(x₁[r]-x₂[c])*(x₁[r]-x₂[c]))
		end
	end
	return(K)
end
type GP
	xₒ::Vector{Float64}
	yₒ::Vector{Float64}
	nₒ::Int
	σ²::Float64
	ρ²::Float64
	ψ²::Float64
	K::Matrix{Float64}
	KIY::Vector{Float64}
	KI::PDMat
	inputType::ASCIIString
	logDensity::Float64
	function GP(xₒyₒ::Dict{Float64,Float64},  σ²::Float64, ρ²::Float64, ψ²::Float64,inputType::ASCIIString)
		xₒ=collect(keys(xₒyₒ))
		yₒ=collect(values(xₒyₒ))
		nₒ=length(yₒ)
		K=Kernel(xₒ,xₒ,ρ²,ψ²)
		if(inputType=="response")
			KI=PDMat(K+ScalMat(nₒ, 1.0)) #Acquire Idnetity from Response
		elseif(inputType=="function")
			KI=PDMat(K+ScalMat(nₒ, 0.000001))
		else
			println("Type of input must be either response or function")
		end
		KIY=\(KI,yₒ)
		logDensity=-0.5*nₒ*log(σ²)-0.5*logdet(KI)-0.5*(1/σ²)*dot(yₒ,KIY)
		#=logDensity=logpdf(MvNormal(σ²*KI),zeros(nₒ))=#
		gp = new(xₒ, yₒ, nₒ, σ²,ρ²,ψ²,K,KIY,KI,inputType,logDensity)
		return gp
	end
	function GP(xₒ::Vector{Float64}, yₒ::Vector{Float64},  σ²::Float64, ρ²::Float64, ψ²::Float64,inputType::ASCIIString)
		nₒ=length(yₒ)
		K=Kernel(xₒ,xₒ,ρ²,ψ²)
		if(inputType=="response")
			KI=PDMat(K+ScalMat(nₒ, 1.0)) #Acquire Idnetity from Response
		elseif(inputType=="function")
			KI=PDMat(K+ScalMat(nₒ, 0.000001))
		else
			println("Type of input must be either response or function")
		end
		KIY=\(KI,yₒ)
		#=logDensity=logpdf(MvNormal(σ²*KI),zeros(nₒ))=#
		logDensity=-0.5*nₒ*log(σ²)-0.5*logdet(KI)-0.5*(1/σ²)*dot(yₒ,KIY)
		gp = new(xₒ, yₒ, nₒ, σ²,ρ²,ψ²,K,KIY,KI,inputType,logDensity)
		return gp
	end
end
function rand(d::GP)
	μ=d.K*d.KIY
	Σ=d.σ²*(d.K-X_invA_Xt(d.KI,d.K))
	ϵ=0.0000001
	ϵIₒₒ=ScalMat(d.nₒ,ϵ)
	return Dict(zip(d.xₒ,rand(MvNormal(μ,Σ+ϵIₒₒ))))
end
function rand(d::GP,xₚ::Vector{Float64})
	Kₚₚ=Kernel(xₚ,xₚ,d.ρ²,d.ψ²)
	Kₚₒ=Kernel(xₚ,d.xₒ,d.ρ²,d.ψ²)
	μ=Kₚₒ*d.KIY
	Σ=d.σ²*(Kₚₚ-X_invA_Xt(d.KI,Kₚₒ))
	ϵ=0.0000001
	nₚ=length(xₚ)
	ϵIₚₚ=ScalMat(nₚ,ϵ)
	fₚ=rand(MvNormal(μ,Σ+ϵIₚₚ))
	return Dict(zip(xₚ,fₚ))
end	

