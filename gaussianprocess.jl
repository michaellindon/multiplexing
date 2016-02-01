function Kernel(x₁::Vector{Float64},x₂::Vector{Float64},ł,ρ²)
	λ=sqrt(2*ν)/ł
	K=Array(Float64,length(x₁),length(x₂))
	for r=1:length(x₁)
		for c=1:length(x₂)
			τ=abs(x₁[r]-x₂[c])
			if(τ!=0.0)
				K[r,c]=ρ²*(2.0^(1.0-ν))*((λ*τ)^ν)*besselk(ν,λ*τ)/gamma(ν)
			else
				K[r,c]=ρ²
			end
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
	ł::Float64
	Kₒₒ::Matrix{Float64}
	Kₒₒsvd
	inputType::ASCIIString
	logDensity::Float64
	function GP(xₒyₒ::Dict{Float64,Float64},  σ²::Float64, ρ²::Float64, ł::Float64,inputType::ASCIIString)
		xₒ=collect(keys(xₒyₒ))
		yₒ=collect(values(xₒyₒ))
		nₒ=length(yₒ)
		Kₒₒ=Kernel(xₒ,xₒ,ł,ρ²)
		Kₒₒsvd=rsvd(Kₒₒ,20)
		#=logDensity=-0.5*nₒ*log(σ²)-0.5*logdet(KI)-0.5*(1/σ²)*dot(yₒ,KIY)=#
		logDensity=0
		gp = new(xₒ, yₒ, nₒ, σ²,ł,ρ²,Kₒₒ,Kₒₒsvd,inputType,logDensity)
		return gp
	end
	function GP(xₒ::Vector{Float64}, yₒ::Vector{Float64},  σ²::Float64, ρ²::Float64, ł::Float64,inputType::ASCIIString)
		nₒ=length(yₒ)
		Kₒₒ=Kernel(xₒ,xₒ,ł,ρ²)
		Kₒₒsvd=rsvd(Kₒₒ,20)
		#=logDensity=-0.5*nₒ*log(σ²)-0.5*logdet(KI)-0.5*(1/σ²)*dot(yₒ,KIY)=#
		logDensity=0
		gp = new(xₒ, yₒ, nₒ, σ²,ł,ρ²,Kₒₒ,Kₒₒsvd,inputType,logDensity)
		return gp
	end
end
function rand(d::GP)
	U=d.Kₒₒsvd.U
	S=d.Kₒₒsvd.S
	nₒ=d.nₒ
	nₐ=length(S)
	yₒ=d.yₒ
	xₒ=d.xₒ
	L=U*Diagonal(sqrt(S))
	LtY=L'*yₒ
	LtL=L'*L;
	W=PDMat(I+LtL)
	fₒ=L*LtY-L*(LtL*\(W,LtY))+L*chol(eye(nₐ)-LtL+(X_invA_Xt(W,LtL)))*rand(Normal(0,1),nₐ)
	return Dict(zip(xₒ,fₒ))
end
function rand(d::GP,xₚ::Vector{Float64})
	xₒ=d.xₒ
	yₒ=d.yₒ
	ρ²=d.ρ²
	ł=d.ł
	Kₒₒ=d.Kₒₒ
	nₒ=d.nₒ
	nₚ=length(xₚ)
	if(d.inputType=="function")
		ϵ=10000;
	else
		ϵ=1;
	end
	Kₚₚ=Kernel(xₚ,xₚ,ł,ρ²)
	Kₚₒ=Kernel(xₚ,xₒ,ł,ρ²)
	Ksvd=rsvd(vcat(hcat(Kₒₒ,Kₚₒ'),hcat(Kₚₒ,Kₚₚ)),20)
	nₐ=length(Ksvd.S)
	L=Ksvd.U*Diagonal(sqrt(Ksvd.S))
	J=L[1:nₒ,:]
	L=L[(nₒ+1):(nₒ+nₚ),:]
	JtY=J'*yₒ
	JtJ=J'*J
	W=PDMat(I+ϵ*JtJ)
	fₚ=ϵ*L*JtY-ϵ*L*(JtJ*\(W,JtY))*ϵ+L*chol(eye(nₐ)-ϵ*JtJ+ϵ*(X_invA_Xt(W,JtJ)*ϵ))*rand(Normal(0,1),nₐ)
	return Dict(zip(xₚ,fₚ))
end	

