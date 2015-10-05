function Kernel(x₁,x₂,ρ²,ψ²)
	K=Array(Float64,size(x₁,1),size(x₂,1))
	Ψ²=Diagonal(ψ²)
	for r=1:size(x₁,1)
		for c=1:size(x₂,1)
			K[r,c]=ρ²*exp(-((x₁[r,:]-x₂[c,:])*Ψ²*(x₁[r,:]-x₂[c,:])')[1])
		end
	end
	return(K)
end

function PredictGP(xₚ,zₒ,xₒ,σ²,ρ²,ψ²,inputType)
	ϵ=0.000000001
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
	yₚ=fₚ+rand(Normal(0,sqrt(σ²)),nₚ)
	return yₚ,fₚ
end	
function PosteriorGP(yₒ,xₒ,σ²,ρ²,ψ²)
	ϵ=0.000000000001
	Kₒₒ=Kernel(xₒ,xₒ,ρ²,ψ²)
	nₒ=size(xₒ,1)
	Iₒₒ=eye(nₒ)
	μ=Kₒₒ*\(Kₒₒ+Iₒₒ,yₒ)
	Σ=Kₒₒ-Kₒₒ*\(Kₒₒ+Iₒₒ,Kₒₒ')
	fₒ=rand(MvNormal(μ,σ²*Σ+ϵ*Iₒₒ),1)[:,1]
	return fₒ
end	
