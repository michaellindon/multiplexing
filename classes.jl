
type Atrial
	id::Int64
	Tobs::Float64
	y₀::SortedDict
	ξ₀ₐ::Dict
	ξ₀ᵣ::Dict
end

function Atrial(id,ξ₀ₐ,Tobs) 
	y₀=SortedDict(Dict{Float64,Float64}())
	for key in collect(keys(ξ₀ₐ))
		y₀[key]=0.0
	end
	return Atrial(id,Tobs,y₀,ξ₀ₐ,Dict())
end

function params(trial::Atrial)
	return (trial.id,trial.Tobs,trial.y₀,trial.ξ₀ₐ,trial.ξ₀ᵣ)
end


type Btrial
	id::Int64
	Tobs::Float64
	y₁::SortedDict
	ξ₁ₐ::Dict
	ξ₁ᵣ::Dict
end

function params(trial::Btrial)
	return (trial.id,trial.Tobs,trial.y₁,trial.ξ₁ₐ,trial.ξ₁ᵣ)
end

function Btrial(id,ξ₁ₐ,Tobs) 
	y₁=SortedDict(Dict{Float64,Float64}())
	for key in collect(keys(ξ₁ₐ))
		y₁[key]=0.0
	end
	return Btrial(id,Tobs,y₁,ξ₁ₐ,Dict())
end

type ABtrial
	id::Int64
	Tobs::Float64
	μg::Float64
	y₀::SortedDict
	y₁::SortedDict
	yg::SortedDict
	ξ₀ₐᵣ::Dict
	ξ₀ᵣᵣ::Dict
	ξ₁ₐᵣ::Dict
	ξ₁ᵣᵣ::Dict
	ξ₀ₐₐ::Dict
	ξ₁ₐₐ::Dict
	𝑇::Array{Float64,1}
	g::SortedDict
	gᵧ::Int64
	μprior::Float64
end

function params(trial::ABtrial)
	return (trial.id,trial.Tobs,trial.μg,trial.y₀,trial.y₁,trial.yg,trial.ξ₀ₐᵣ,trial.ξ₀ᵣᵣ,trial.ξ₁ₐᵣ,trial.ξ₁ᵣᵣ,trial.ξ₀ₐₐ,trial.ξ₁ₐₐ,trial.𝑇,trial.g,trial.gᵧ,trial.μprior)
end

function ABtrial(id,𝑇,Tobs) 
	g=SortedDict(Dict{Float64,Array{Float64,1}}())
	y₀=SortedDict(Dict{Float64,Float64}())
	y₁=SortedDict(Dict{Float64,Float64}())
	yg=SortedDict(Dict{Float64,Float64}())
	for t in 𝑇
		g[t]=zeros(Float64,3)
		y₀[t]=0.0
		y₁[t]=0.0
		yg[t]=0.0
	end
	return ABtrial(id,Tobs,0.0,y₀,y₁,yg,Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),𝑇,g,1,0)
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
	return ABtrial(id,Tobs,0.0,y₀,y₁,yg,Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),𝑇,g,1,0)
end
