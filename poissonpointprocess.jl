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



function Ξ!(trial::Atrial,m₀,λ₀,Λ)
	empty!(trial.y₀)

	#ξ₀ₐ
	for t in trial.ξ₀ₐ
		trial.y₀[t]=rand(Truncated(Normal(m₀(t),1),0,Inf))
	end

	#ξ₀ᵣ
	trial.ξ₀ᵣ=rand(PPProcess(x-> Λ-λ₀(x),Λ),0,trial.Tobs)
	for t in trial.ξ₀ᵣ
		trial.y₀[t]=rand(Truncated(Normal(m₀(t),1),-Inf,0))
	end

end

function Ξ!(trial::Btrial,m₁,λ₁,Λ)
	empty!(trial.y₁)

	#ξ₁ₐ
	for t in trial.ξ₁ₐ
		trial.y₁[t]=rand(Truncated(Normal(m₁(t),1),0,Inf))
	end

	#ξ₁ᵣ
	trial.ξ₁ᵣ=rand(PPProcess(x-> Λ-λ₁(x),Λ),0,trial.Tobs)
	for t in trial.ξ₁ᵣ
		trial.y₁[t]=rand(Truncated(Normal(m₁(t),1),-Inf,0))
	end
end

function Ξ!(trial::ABtrial,m₀,λ₀,m₁,λ₁,Λ)
	empty!(trial.y₀)
	empty!(trial.y₁)
	empty!(trial.yg)

	trial.ξ₀ₐᵣ=rand(PPProcess(x->(1-trial.α(x))*λ₀(x),Λ),0,trial.Tobs)
	for t in trial.ξ₀ₐᵣ
			trial.y₀[t]=rand(Truncated(Normal(m₀(t),1),0,Inf))
			trial.yg[t]=rand(Truncated(Normal(trial.mₐ(t),1),-Inf,0))
	end

	trial.ξ₀ᵣᵣ=rand(PPProcess(x-> Λ-λ₀(x),Λ),0,trial.Tobs)
	for t in trial.ξ₀ᵣᵣ
			trial.y₀[t]=rand(Truncated(Normal(m₀(t),1),-Inf,0))
	end

	trial.ξ₁ₐᵣ=rand(PPProcess(x-> trial.α(x)*λ₁(x),Λ),0,trial.Tobs)
	for t in trial.ξ₁ₐᵣ
			trial.y₁[t]=rand(Truncated(Normal(m₁(t),1),0,Inf))
			trial.yg[t]=rand(Truncated(Normal(trial.mₐ(t),1),0,Inf))
	end

	trial.ξ₁ᵣᵣ=rand(PPProcess(x-> Λ-λ₁(x),Λ),0,trial.Tobs)
	for t in trial.ξ₁ᵣᵣ
			trial.y₁[t]=rand(Truncated(Normal(m₁(t),1),-Inf,0))
	end

	for t in trial.𝑇
		denominator=trial.α(t)*λ₀(t)+(1-trial.α(t))*λ₁(t)
		if(rand(Bernoulli((1-trial.α(t))*λ₁(t)/(denominator)))[1]==1.0)
			trial.y₁[t]=rand(Truncated(Normal(m₁(t),1),0,Inf))
			trial.yg[t]=rand(Truncated(Normal(trial.mₐ(t),1),-Inf,0))
		else
			trial.y₀[t]=rand(Truncated(Normal(m₀(t),1),0,Inf))
			trial.yg[t]=rand(Truncated(Normal(trial.mₐ(t),1),0,Inf))
		end
	end
end

function 𝐺!(trial::ABtrial,σ²,σ²ₘ,łg,ρ²g,p)

	trial.μg=rand(mu(trial.yg,trial.gᵧ,σ²,łg,ρ²g,σ²ₘ,trial.μprior))
	logodds=(sslogdensity(trial.yg,1,trial.μg,σ²,łg,ρ²g)-sslogdensity(trial.yg,0,trial.μg,σ²,łg,ρ²g))
	odds=exp(logodds)*p/(1-p)
	if(odds==Inf)
		trial.gᵧ=1
	elseif(odds==-Inf)
		trial.gᵧ=0
	else
		trial.gᵧ=rand(Bernoulli(odds/(1+odds)))
	end
	trial.g = (trial.gᵧ==1) ? lazyGP(FFBS(trial.yg,trial.μg,σ²,łg,ρ²g),łg) : x->0
	trial.mₐ = x -> trial.μg + trial.g(x)
	trial.α = x -> Φ(trial.mₐ(x))
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


function generate(Tobs,m₀,λ₀,m₁,λ₁,Λ,p,łg,ρ²g,prior)
	μg=rand(prior[:μg])
	gᵧ=rand(Bernoulli(p))
	g = (gᵧ==1) ? lazyGP(realization(statcov(łg,ρ²g),transition(łg),innovation(łg,ρ²g),collect(0:Tobs/1000:Tobs)),łg) : x->0
	mₐ= x->μg+g(x)
	α= x->Φ(mₐ(x))
	y₀=SortedDict(Dict{Float64}{Float64}())
	y₁=SortedDict(Dict{Float64}{Float64}())
	yg=SortedDict(Dict{Float64}{Float64}())
	ξ₀ₐₐ=rand(PPProcess(x->α(x)*λ₀(x),Λ),0,Tobs)
	for t in ξ₀ₐₐ
		y₀[t]=rand(Truncated(Normal(m₀(t),1),0,Inf))
		yg[t]=rand(Truncated(Normal(mₐ(t),1),0,Inf))
	end
	ξ₁ₐₐ=rand(PPProcess(x->(1-α(x))*λ₁(x),Λ),0,Tobs)
	for t in ξ₁ₐₐ
		y₁[t]=rand(Truncated(Normal(m₁(t),1),0,Inf))
		yg[t]=rand(Truncated(Normal(mₐ(t),1),-Inf,0))
	end
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
	𝑇=sort(union(ξ₀ₐₐ,ξ₁ₐₐ))
	ABtrial(1,Tobs,μg,y₀,y₁,yg,ξ₀ₐᵣ,ξ₀ᵣᵣ,ξ₁ₐᵣ,ξ₁ᵣᵣ,ξ₀ₐₐ,ξ₁ₐₐ,𝑇,g,mₐ,α,gᵧ,0)
end
