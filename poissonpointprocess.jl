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
	Tₚ=rand(Uniform(a,b),nₚ)
	times=Array(Float64,0)
	for t in Tₚ
		if(rand(Bernoulli(Ξ.λ(t)/Ξ.Λ))[1]==1.0)
			push!(times,t)
		end
	end
	return sort(times)
end

function Ξ₀ₐ(ξ₀ₐ,μ₀,f₀)
	for t in keys(ξ₀ₐ)
		ξ₀ₐ[t]["yf"]=rand(Truncated(Normal(μ₀+f₀[t][1],1),0,Inf))
	end
end

function Ξ₁ₐ(ξ₁ₐ,μf₁,f₁)
	for t in keys(ξ₁ₐ)
		ξ₁ₐ[t]["yf"]=rand(Truncated(Normal(μf₁+f₁[t][1],1),0,Inf))
	end
end

function Ξ₀ᵣ(μ₀,f₀,ł₀,ρ²₀,Ξₚ,Tobs)
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₀,Tₚ,ł₀,ρ²₀)
	ξ₀ᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	sizehint!(ξ₀ᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(μ₀+fₚ[t][1])))[1]==1.0)
			ξ₀ᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(μ₀+fₚ[t][1],1),-Inf,0)))
		end
	end
	return ξ₀ᵣ
end

function Ξ₁ᵣ(μf₁,f₁,ł₁,ρ²₁,Ξₚ,Tobs)
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₁,Tₚ,ł₁,ρ²₁)
	ξ₁ᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	sizehint!(ξ₁ᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(μf₁+fₚ[t][1])))[1]==1.0)
			ξ₁ᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(μf₁+fₚ[t][1],1),-Inf,0)))
		end
	end
	return ξ₁ᵣ
end


function Ξ₀ₐᵣ(μg,g,łg,ρ²g,μ₀,f₀,ł₀,ρ²₀,Ξₚ,Tobs)
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,łg,ρ²g)
	fₚ=FFBS2(f₀,Tₚ,ł₀,ρ²₀)
	ξ₀ₐᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	sizehint!(ξ₀ₐᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(μg+gₚ[t][1]))*Φ(μ₀+fₚ[t][1])))[1]==1.0)
			ξ₀ₐᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(μ₀+fₚ[t][1],1),0,Inf)),"yg"=>rand(Truncated(Normal(μg+gₚ[t][1],1),-Inf,0)))
		end
	end
	return ξ₀ₐᵣ
end

function Ξ₀ᵣᵣ(μ₀,f₀,ł₀,ρ²₀,Ξₚ,Tobs)
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₀,Tₚ,ł₀,ρ²₀)
	ξ₀ᵣᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	sizehint!(ξ₀ᵣᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(μ₀+fₚ[t][1])))[1]==1.0)
			ξ₀ᵣᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(μ₀+fₚ[t][1],1),-Inf,0)))
		end
	end
	return ξ₀ᵣᵣ
end

function Ξ₁ₐᵣ(μg,g,łg,ρ²g,μ₁,f₁,ł₁,ρ²₁,Ξₚ,Tobs)
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,łg,ρ²g)
	fₚ=FFBS2(f₁,Tₚ,ł₁,ρ²₁)
	ξ₁ₐᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	sizehint!(ξ₁ₐᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(Φ(μg+gₚ[t][1])*Φ(μ₁+fₚ[t][1])))[1]==1.0)
			ξ₁ₐᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(μ₁+fₚ[t][1],1),0,Inf)),"yg"=>rand(Truncated(Normal(μg+gₚ[t][1],1),0,Inf)))
		end
	end
	return ξ₁ₐᵣ
end

function Ξ₁ᵣᵣ(μ₁,f₁,ł₁,ρ²₁,Ξₚ,Tobs)
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₁,Tₚ,ł₁,ρ²₁)
	ξ₁ᵣᵣ=Dict{Float64,Dict{UTF8String,Float64}}()
	sizehint!(ξ₁ᵣᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(μ₁+fₚ[t][1])))[1]==1.0)
			ξ₁ᵣᵣ[t]=Dict{UTF8String,Float64}("yf"=>rand(Truncated(Normal(μ₁+fₚ[t][1],1),-Inf,0)))
		end
	end
	return ξ₁ᵣᵣ
end


function Ξ₀ₐₐ₁ₐₐ(μ₀,f₀,ł₀,ρ²₀,μ₁,f₁,ł₁,ρ²₁,μg,g,łg,ρ²g,Td)
	ξ₀ₐₐ=Dict{Float64,Dict{UTF8String,Float64}}()
	sizehint!(ξ₀ₐₐ,length(Td))
	ξ₁ₐₐ=Dict{Float64,Dict{UTF8String,Float64}}()
	sizehint!(ξ₁ₐₐ,length(Td))
	gₚ=FFBS2(g,Td,łg,ρ²g)
	f₀ₚ=FFBS2(f₀,Td,ł₀,ρ²₀)
	f₁ₚ=FFBS2(f₁,Td,ł₁,ρ²₁)
	for t in Td
		denominator=Φ(μg+gₚ[t][1])*Φ(μ₀+f₀ₚ[t][1])+(1-Φ(μg+gₚ[t][1]))*Φ(μ₁+f₁ₚ[t][1])
		if(rand(Bernoulli((1-Φ(μg+gₚ[t][1]))*Φ(μ₁+f₁ₚ[t][1])/denominator))[1]==1.0)
			ξ₁ₐₐ[t]=Dict("yf"=>rand(Truncated(Normal(μ₁+f₁ₚ[t][1],1),0,Inf)), "yg"=>rand(Truncated(Normal(μg+gₚ[t][1],1),-Inf,0)))
		else
			ξ₀ₐₐ[t]=Dict("yf"=>rand(Truncated(Normal(μ₀+f₀ₚ[t][1],1),0,Inf)), "yg"=>rand(Truncated(Normal(μg+gₚ[t][1],1),0,Inf)))
		end
	end
	return ξ₀ₐₐ,ξ₁ₐₐ
end

function Y₀Y₁Yg(ξ₀ₐ,ξ₀ᵣ,ξ₁ₐ,ξ₁ᵣ,ξ₀ₐₐ,ξ₀ₐᵣ,ξ₁ₐₐ,ξ₁ₐᵣ,ξ₀ᵣᵣ,ξ₁ᵣᵣ)
	y₀=Dict{Float64,Float64}()
	sizehint!(y₀,length(ξ₀ₐ)+length(ξ₀ᵣ)+length(ξ₀ₐₐ)+length(ξ₀ₐᵣ)+length(ξ₀ᵣᵣ))
	y₁=Dict{Float64,Float64}()
	sizehint!(y₁,length(ξ₁ₐ)+length(ξ₁ᵣ)+length(ξ₁ₐₐ)+length(ξ₁ₐᵣ)+length(ξ₁ᵣᵣ))
	yg=Dict{Float64,Float64}()
	sizehint!(yg,length(ξ₀ₐₐ)+length(ξ₀ₐᵣ)+length(ξ₁ₐₐ)+length(ξ₁ₐᵣ)+length(ξ₀ᵣᵣ)+length(ξ₁ᵣᵣ))
	for t in keys(ξ₀ₐ)
		y₀[t]=ξ₀ₐ[t]["yf"]
	end
	for t in keys(ξ₀ᵣ)
		y₀[t]=ξ₀ᵣ[t]["yf"]
	end
	for t in keys(ξ₁ₐ)
		y₁[t]=ξ₁ₐ[t]["yf"]
	end
	for t in keys(ξ₁ᵣ)
		y₁[t]=ξ₁ᵣ[t]["yf"]
	end
	for t in keys(ξ₀ₐₐ)
		y₀[t]=ξ₀ₐₐ[t]["yf"]
		yg[t]=ξ₀ₐₐ[t]["yg"]
	end
	for t in keys(ξ₀ₐᵣ)
		y₀[t]=ξ₀ₐᵣ[t]["yf"]
		yg[t]=ξ₀ₐᵣ[t]["yg"]
	end
	for t in keys(ξ₁ₐₐ)
		y₁[t]=ξ₁ₐₐ[t]["yf"]
		yg[t]=ξ₁ₐₐ[t]["yg"]
	end
	for t in keys(ξ₁ₐᵣ)
		y₁[t]=ξ₁ₐᵣ[t]["yf"]
		yg[t]=ξ₁ₐᵣ[t]["yg"]
	end
	for t in keys(ξ₀ᵣᵣ)
		y₀[t]=ξ₀ᵣᵣ[t]["yf"]
	end
	for t in keys(ξ₁ᵣᵣ)
		y₁[t]=ξ₁ᵣᵣ[t]["yf"]
	end
	return y₀,y₁,yg
end


function Ξ!(trial::Atrial,μ₀,μ₀ₜ,f₀,ł₀,ρ²₀,Ξₚ)
	(id,Tobs,y₀,ξ₀ₐ,ξ₀ᵣ)=params(trial)
	empty!(y₀)

	#ξ₀ₐ
	for t in keys(ξ₀ₐ)
		y₀[t]=rand(Truncated(Normal(μ₀+μ₀ₜ(t)+f₀[t][1],1),0,Inf))
		ξ₀ₐ[t]["yf"]=y₀[t]
	end

	#ξ₀ᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₀,Tₚ,ł₀,ρ²₀)
	empty!(ξ₀ᵣ)
	sizehint!(ξ₀ᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(μ₀+μ₀ₜ(t)+fₚ[t][1])))[1]==1.0)
			y₀[t]=rand(Truncated(Normal(μ₀+μ₀ₜ(t)+fₚ[t][1],1),-Inf,0))
			ξ₀ᵣ[t]=Dict{UTF8String,Float64}("yf"=>y₀[t])
		end
	end

end

function Ξ!(trial::Btrial,μ₁,μ₁ₜ,f₁,ł₁,ρ²₁,Ξₚ)
	(id,Tobs,y₁,ξ₁ₐ,ξ₁ᵣ)=params(trial)
	empty!(y₁)

	#ξ₁ₐ
	for t in keys(ξ₁ₐ)
		y₁[t]=rand(Truncated(Normal(μ₁+μ₁ₜ(t)+f₁[t][1],1),0,Inf))
		ξ₁ₐ[t]["yf"]=y₁[t]
	end

	#ξ₁ᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₁,Tₚ,ł₁,ρ²₁)
	empty!(ξ₁ᵣ)
	sizehint!(ξ₁ᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(μ₁+μ₁ₜ(t)+fₚ[t][1])))[1]==1.0)
			y₁[t]=rand(Truncated(Normal(μ₁+μ₁ₜ(t)+fₚ[t][1],1),-Inf,0))
			ξ₁ᵣ[t]=Dict{UTF8String,Float64}("yf"=>y₁[t])
		end
	end
end

function Ξ!(trial::ABtrial,m₀,λ₀,m₁,λ₁,łg,ρ²g,Λ)
	(id,Tobs,μg,y₀,y₁,yg,ξ₀ₐᵣ,ξ₀ᵣᵣ,ξ₁ₐᵣ,ξ₁ᵣᵣ,ξ₀ₐₐ,ξ₁ₐₐ,𝑇,gg,gᵧ,μprior)=params(trial)
	empty!(y₀)
	empty!(y₁)
	empty!(yg)

	g=lazyGP(gg,łg)

	mₐ = x->μg+g(x)
	α = x->Φ(mₐ(x))
	
	ξ₀ₐᵣ=rand(PPProcess(x->(1-α(x))*λ₀(x),Λ),0,Tobs)
	for t in ξ₀ₐᵣ
			y₀[t]=rand(Truncated(Normal(m₀(t),1),0,Inf))
			yg[t]=rand(Truncated(Normal(mₐ(t),1),-Inf,0))
	end

	ξ₀ᵣᵣ=rand(PPProcess(x-> Λ-λ₀(x),Λ),0,Tobs)
	for t in ξ₀ᵣᵣ
			y₀[t]=rand(Truncated(Normal(m₀(t),1),-Inf,0))
	end

	ξ₁ₐᵣ=rand(PPProcess(x-> α(x)*λ₁(x),Λ),0,Tobs)
	for t in ξ₁ₐᵣ
			y₁[t]=rand(Truncated(Normal(m₁(t),1),0,Inf))
			yg[t]=rand(Truncated(Normal(mₐ(t),1),0,Inf))
	end

	ξ₁ᵣᵣ=rand(PPProcess(x-> Λ-λ₁(x),Λ),0,Tobs)
	for t in ξ₁ᵣᵣ
			y₁[t]=rand(Truncated(Normal(m₁(t),1),-Inf,0))
	end

	for t in 𝑇
		denominator=α(t)*λ₀(t)+(1-α(t))*λ₁(t)
		if(rand(Bernoulli((1-α(t))*λ₁(t)/(denominator)))[1]==1.0)
			y₁[t]=rand(Truncated(Normal(m₁(t),1),0,Inf))
			yg[t]=rand(Truncated(Normal(mₐ(t),1),-Inf,0))
		else
			y₀[t]=rand(Truncated(Normal(m₀(t),1),0,Inf))
			yg[t]=rand(Truncated(Normal(mₐ(t),1),0,Inf))
		end
	end
end

function 𝐺!(trial::ABtrial,σ²,σ²ₘ,łg,ρ²g,p)
	(id,Tobs,μg,y₀,y₁,yg,ξ₀ₐᵣ,ξ₀ᵣᵣ,ξ₁ₐᵣ,ξ₁ᵣᵣ,ξ₀ₐₐ,ξ₁ₐₐ,𝑇,g,gᵧ,μprior)=params(trial);

	#=trial.μg=rand(mu(trial.yg,trial.gᵧ,σ²,łg,ρ²g,σ²ₘ,trial.μprior))=#
	#=logodds=(sslogdensity(yg,1,μg,σ²,łg,ρ²g)-sslogdensity(yg,0,μg,σ²,łg,ρ²g))=#
	#=odds=exp(logodds)*p/(1-p)=#
	#=if(odds==Inf)=#
		#=trial.gᵧ=1=#
	#=elseif(odds==-Inf)=#
		#=trial.gᵧ=0=#
	#=else=#
		#=trial.gᵧ=rand(Bernoulli(odds/(1+odds)))=#
	#=end=#
	#=if(trial.gᵧ==1)=#
		trial.g=FFBS(yg,μg,σ²,łg,ρ²g)
	#=else=#
		#=empty!(g)=#
		#=for key in keys(yg)=#
			#=g[key]=zeros(Float64,3)=#
		#=end=#
	#=end=#
end

function Ξ(trial::ABtrial,μ₀,f₀,ł₀,ρ²₀,μ₁,f₁,ł₁,ρ²₁,łg,ρ²g,Ξₚ)
	(id,Tobs,μg,y₀,y₁,yg,ξ₀ₐᵣ,ξ₀ᵣᵣ,ξ₁ₐᵣ,ξ₁ᵣᵣ,ξ₀ₐₐ,ξ₁ₐₐ,𝑇,g,gᵧ)=params(trial)
	empty!(y₀)
	empty!(y₁)
	empty!(yg)

	#ξ₀ₐᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,łg,ρ²g)
	fₚ=FFBS2(f₀,Tₚ,ł₀,ρ²₀)
	empty!(ξ₀ₐᵣ)
	sizehint!(ξ₀ₐᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli((1-Φ(μg+gₚ[t][1]))*Φ(μ₀+fₚ[t][1])))[1]==1.0)
			y₀[t]=rand(Truncated(Normal(μ₀+fₚ[t][1],1),0,Inf))
			yg[t]=rand(Truncated(Normal(μg+gₚ[t][1],1),-Inf,0))
			ξ₀ₐᵣ[t]=Dict{UTF8String,Float64}("yf"=>y₀[t],"yg"=>yg[t])
		end
	end

	#ξ₀ᵣᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₀,Tₚ,ł₀,ρ²₀)
	empty!(ξ₀ᵣᵣ)
	sizehint!(ξ₀ᵣᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(μ₀+fₚ[t][1])))[1]==1.0)
			y₀[t]=rand(Truncated(Normal(μ₀+fₚ[t][1],1),-Inf,0))
			ξ₀ᵣᵣ[t]=Dict{UTF8String,Float64}("yf"=>y₀[t])
		end
	end

	#sample ξ₁ₐᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	gₚ=FFBS2(g,Tₚ,łg,ρ²g)
	fₚ=FFBS2(f₁,Tₚ,ł₁,ρ²₁)
	empty!(ξ₁ₐᵣ)
	sizehint!(ξ₁ₐᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(Φ(μg+gₚ[t][1])*Φ(μ₁+fₚ[t][1])))[1]==1.0)
			y₁[t]=rand(Truncated(Normal(μ₁+fₚ[t][1],1),0,Inf))
			yg[t]=rand(Truncated(Normal(μg+gₚ[t][1],1),0,Inf))
			ξ₁ₐᵣ[t]=Dict{UTF8String,Float64}("yf"=>y₁[t],"yg"=>yg[t])
		end
	end

	#sample ξ₁ᵣᵣ
	Tₚ=rand(Ξₚ,0,Tobs)
	fₚ=FFBS2(f₁,Tₚ,ł₁,ρ²₁)
	empty!(ξ₁ᵣᵣ)
	sizehint!(ξ₁ᵣᵣ,length(Tₚ))
	for t in Tₚ
		if(rand(Bernoulli(1-Φ(μ₁+fₚ[t][1])))[1]==1.0)
			y₁[t]=rand(Truncated(Normal(μ₁+fₚ[t][1],1),-Inf,0))
			ξ₁ᵣᵣ[t]=Dict{UTF8String,Float64}("yf"=>y₁[t])
		end
	end

	#sample ξ₀ₐₐ, ξ₁ₐₐ
	empty!(ξ₀ₐₐ)
	sizehint!(ξ₀ₐₐ,length(𝑇))
	empty!(ξ₁ₐₐ)
	sizehint!(ξ₁ₐₐ,length(𝑇))
	#=f₀ₚ=FFBS2(f₀,𝑇,ł₀,ρ²₀)=#
	#=f₁ₚ=FFBS2(f₁,𝑇,ł₁,ρ²₁)=#
	for t in 𝑇
		denominator=Φ(μg+g[t][1])*Φ(μ₀+f₀ₚ[t][1])+(1-Φ(μg+g[t][1]))*Φ(μ₁+f₁ₚ[t][1])
		if(rand(Bernoulli((1-Φ(μg+g[t][1]))*Φ(μ₁+f₁ₚ[t][1])/denominator))[1]==1.0)
			y₁[t]=rand(Truncated(Normal(μ₁+f₁ₚ[t][1],1),0,Inf))
			yg[t]=rand(Truncated(Normal(μg+g[t][1],1),-Inf,0))
			ξ₁ₐₐ[t]=Dict("yf"=>y₁[t], "yg"=>yg[t])
		else
			y₀[t]=rand(Truncated(Normal(μ₀+f₀ₚ[t][1],1),0,Inf))
			yg[t]=rand(Truncated(Normal(μg+g[t][1],1),0,Inf))
			ξ₀ₐₐ[t]=Dict("yf"=>y₀[t], "yg"=>yg[t])
		end
	end

	return ABtrial(id,Tobs,μg,y₀,y₁,yg,ξ₀ₐᵣ,ξ₀ᵣᵣ,ξ₁ₐᵣ,ξ₁ᵣᵣ,ξ₀ₐₐ,ξ₁ₐₐ,𝑇,g,gᵧ)
end



function lazyGP(f,łs)
	ł=sqrt(5.0)/łs;
	function predict(t)
		if(haskey(f,t))
			return f[t][1]
		end
		𝑍=rand(Normal(0,1),3)
		xout=zeros(Float64,3)
		fore=searchsortedlast(f,t); #This routine returns the semitoken of the last item in the container whose key is less than or equal to t. If no such key, then before-start semitoken is returned. 
		aft=searchsortedfirst(f,t); #This routine returns the semitoken of the first item in the container whose key is greater than or equal to t. If no such key, then past-end semitoken is returned. 
		if(fore==beforestartsemitoken(f)) #No key less than or equal to t
			t₁,x₁=deref((f,advance((f,fore))))
			ccall((:Backwards, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Float64,Float64,Ref{Cdouble},Float64),xout,𝑍,t,t₁,x₁,ł)
			f[t]=xout;
		elseif(aft==pastendsemitoken(f)) #No key greater than or equal to t
			tₙ,xₙ=deref((f,regress((f,aft))))
			ccall((:Forwards, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Float64,Float64,Ref{Cdouble},Float64),xout,𝑍,t,tₙ,xₙ,ł)
			f[t]=xout;
		else
			(tl,vl)=deref((f,fore))
			(tr,vr)=deref((f,aft))
			ccall((:ThreePoint, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Float64,Float64,Ref{Cdouble},Float64,Ref{Cdouble},Float64),xout,𝑍,t,tl,vl,tr,vr,ł)
			f[t]=xout;
		end
		return f[t][1]
	end
end

