immutable SpikeAndSlab <: ContinuousUnivariateDistribution
	Probability::Array{Float64}
	Slab::Distribution

	function SpikeAndSlab(Probability::Array{Float64}, Slab::Distribution)
		new(Probability, Slab)
	end
end

params(d::SpikeAndSlab) = (d.Probability, d.Slab)

function rand(d::SpikeAndSlab)
	(Probability,Slab)=params(d)
	if(rand(Uniform(0,1))<Probability[1])
		return 0.0
	else
		return rand(Slab)
	end
end

function pdf(d::SpikeAndSlab,x::Float64)
	(Probability,Slab)=params(d)
	if(x==0.0)
		return Probability[1]
	else
		return Probability[2]*pdf(Slab,x)
	end
end


function q(θ₁,θ₂,d::SpikeAndSlab)
	(Probability,Slab)=params(d)
	if(θ₁==0.0)
		return Probability[1]
	else
		return Probability[2]*pdf(Slab,θ₁-θ₂)
	end
end
