
type Atrial
	id::Int64
	Tobs::Float64
	y₀::SortedDict
	ξ₀ₐ::Array{Float64,1}
	ξ₀ᵣ::Array{Float64,1}
end

function Atrial(id,ξ₀ₐ,Tobs) 
	y₀=SortedDict(Dict{Float64,Float64}())
	for t in ξ₀ₐ
		y₀[t]=0.0
	end
	return Atrial(id,Tobs,y₀,ξ₀ₐ,[])
end

#=function params(trial::Atrial)=#
	#=return (trial.id,trial.Tobs,trial.y₀,trial.ξ₀ₐ,trial.ξ₀ᵣ)=#
#=end=#


type Btrial
	id::Int64
	Tobs::Float64
	y₁::SortedDict
	ξ₁ₐ::Array{Float64,1}
	ξ₁ᵣ::Array{Float64,1}
end

#=function params(trial::Btrial)=#
	#=return (trial.id,trial.Tobs,trial.y₁,trial.ξ₁ₐ,trial.ξ₁ᵣ)=#
#=end=#

function Btrial(id,ξ₁ₐ,Tobs) 
	y₁=SortedDict(Dict{Float64,Float64}())
	for t in ξ₁ₐ
		y₁[t]=0.0
	end
	return Btrial(id,Tobs,y₁,ξ₁ₐ,[])
end

type ABtrial
	id::Int64
	Tobs::Float64
	μg::Float64
	y₀::SortedDict
	y₁::SortedDict
	yg::SortedDict
	ξ₀ₐᵣ::Array{Float64,1}
	ξ₀ᵣᵣ::Array{Float64,1}
	ξ₁ₐᵣ::Array{Float64,1}
	ξ₁ᵣᵣ::Array{Float64,1}
	ξ₀ₐₐ::Array{Float64,1}
	ξ₁ₐₐ::Array{Float64,1}
	𝑇::Array{Float64,1}
	g::Function 
	mₐ::Function
	α::Function
	gᵧ::Int64
	μprior::Float64
end

#=function params(trial::ABtrial)=#
	#=return (trial.id,trial.Tobs,trial.μg,trial.y₀,trial.y₁,trial.yg,trial.ξ₀ₐᵣ,trial.ξ₀ᵣᵣ,trial.ξ₁ₐᵣ,trial.ξ₁ᵣᵣ,trial.ξ₀ₐₐ,trial.ξ₁ₐₐ,trial.𝑇,trial.g,trial.m,trial.α,trial.gᵧ,trial.μprior)=#
#=end=#

function ABtrial(id,𝑇,Tobs) 
	y₀=SortedDict(Dict{Float64,Float64}())
	y₁=SortedDict(Dict{Float64,Float64}())
	yg=SortedDict(Dict{Float64,Float64}())
	for t in 𝑇
		y₀[t]=0.0
		y₁[t]=0.0
		yg[t]=0.0
	end
	return ABtrial(id,Tobs,0.0,y₀,y₁,yg,[],[],[],[],[],[],𝑇,x->0,x->0,x->0.5,1,0)
end

function ABtrial(id,𝑇,Tobs,g) 
	y₀=SortedDict(Dict{Float64,Float64}())
	y₁=SortedDict(Dict{Float64,Float64}())
	yg=SortedDict(Dict{Float64,Float64}())
	for t in 𝑇
		y₀[t]=0.0
		y₁[t]=0.0
		yg[t]=0.0
	end
	return ABtrial(id,Tobs,0.0,y₀,y₁,yg,[],[],[],[],[],[],𝑇,g,x->0+g(x),x->Φ(0+g(x)),1,0)
end
