function plot(data::Dict{Float64,Array{Float64,2}})
	plot(sort(collect(keys(data))),[data[key][1] for key in sort(collect(keys(data)))])
end

signM=[1.0 ; -1.0  ; 1.0 ; -1.0]
signM=signM*signM'

function SafeMvNormal(Σ)
	(U,S,V)=svd(Σ)
	S=map(x->max(0,x),S)
	U*Diagonal(sqrt(S))*rand(Normal(0,1),3)
end

function SafeInv(Σ)
	(U,S,V)=svd(Σ)
	tolerance=e-10
	S=map(x->( x>tolerance? 1/x: 0 ),S)
	U*Diagonal(S)*V'
end

function realization(Cₛ,transition,innovation,t)
	t=sort(t)
	x=SortedDict(Dict{Float64,Array{Float64,1}}())
	x[t[1]]=(SafeMvNormal(Cₛ))
	#=E,V=eig(Cₛ)=#
	#=x[t[1]]=V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d,1)=#
	for i=2:length(t)
		Δ=t[i]-t[i-1]
		A=transition(Δ)
		Q=innovation(Δ)
		x[t[i]]=A*x[t[i-1]]+(SafeMvNormal(Q))
		#=E,V=eig(ρ²*Q)=#
		#=x[t[i]]=A*x[t[i-1]]+V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d,1)=#
	end
	return(x)
end

function statcorr(ł)
	λ=sqrt(5)/ł
	q=(16*λ^5)/3
	M=[[(3*q)/(16*λ^5),0,-(q/(16*λ^3))]';
	[0,q/(16*λ^3),0]';
	[-(q/(16*λ^3)),0,(3*q)/(16*λ)]']
	return(M)
end

function statcov(ł,ρ²)
	return ρ²*statcorr(ł)
end

function innovation(ł,ρ²)
	λ=sqrt(5)/ł
	q=(16*λ^5)/3
	function innovation(Δ)
		em2Δλ=exp(-2*Δ*λ)
		Q=Array{Float64}(3,3)
		Q[1,1]=(q*(3+em2Δλ*(-3-2*Δ*λ*(3+Δ*λ*(3+Δ*λ*(2+Δ*λ))))))/(16*λ^5);
		Q[1,2]=1/8*em2Δλ*q*Δ^4;
		Q[1,3]=-((q*(-em2Δλ+1+2*Δ*λ*em2Δλ*(-1+Δ*λ*(-1+Δ*λ*(-2+Δ*λ)))))/(16*λ^3));
		Q[2,1]=1/8*em2Δλ*q*Δ^4
		Q[2,2]=(q*(1+em2Δλ*(-1-2*Δ*λ*(1+Δ*λ*(-1+Δ*λ)^2))))/(16*λ^3)
		Q[2,3]=1/8*em2Δλ*q*Δ^2*(-2+Δ*λ)^2
		Q[3,1]=-((q*(-em2Δλ+1+2*Δ*λ*em2Δλ*(-1+Δ*λ*(-1+Δ*λ*(-2+Δ*λ)))))/(16*λ^3))
		Q[3,2]=1/8*em2Δλ*q*Δ^2*(-2+Δ*λ)^2
		Q[3,3]=(q*(3+em2Δλ*(-3-2*Δ*λ*(-5+Δ*λ*(11+Δ*λ*(-6+Δ*λ))))))/(16*λ)
		return ρ²*Q
	end
end

function transition(ł)
	λ=sqrt(5)/ł
	function transition(Δ)
		emΔλ=exp(-Δ*λ)
		Φ=Array{Float64}(3,3)
		Φ[1,1]=1/2*emΔλ*(2+2*Δ*λ+Δ^2*λ^2)
		Φ[1,2]=emΔλ*Δ*(1+Δ*λ)
		Φ[1,3]=1/2*emΔλ*Δ^2
		Φ[2,1]=-(1/2)*emΔλ*Δ^2*λ^3
		Φ[2,2]=-emΔλ*(-1-Δ*λ+Δ^2*λ^2)
		Φ[2,3]=-(1/2)*emΔλ*Δ*(-2+Δ*λ)
		Φ[3,1]=1/2*emΔλ*Δ*λ^3*(-2+Δ*λ)
		Φ[3,2]=emΔλ*Δ*λ^2*(-3+Δ*λ)
		Φ[3,3]=1/2*emΔλ*(2-4*Δ*λ+Δ^2*λ^2)
		return Φ
	end
end

function FFBS(y,μ,σ²,t₀,Cₛ,transition,innovation)
	t=collect(keys(y))
	n=length(t)
	m=Dict{Float64,Array{Float64,2}}(); sizehint!(m,n);
	m[t[1]]=reshape(Cₛ[:,1]*(y[t[1]]-μ)/(σ²+Cₛ[1,1]),d,1)
	M=Dict{Float64,Array{Float64,2}}(); sizehint!(M,n);
	M[t[1]]=Cₛ-Cₛ[:,1]*Cₛ[1,:]/(σ²+Cₛ[1,1])
	AMAQ=Dict{Float64,Array{Float64,2}}(); sizehint!(AMAQ,n);
	A=Dict{Float64,Array{Float64,2}}(); sizehint!(A,n);
	for i=2:n
		Δ=t[i]-t[i-1]
		A[t[i-1]]=transition(Δ)
		Q=innovation(Δ)
		AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q
		m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
		M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
	end
	x=SortedDict(Dict{Float64,Array{Float64,2}}())
	#Backward Sampling
	Σ=M[t[n]]
	Σ=0.5*(Σ+Σ')
	E,V=eig(Σ)
	x[t[n]]=m[t[n]]+V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d)
	#=F=cholfact(Σ+eye(d),:L,Val{true}) =#
	#=x[t[n]]=m[t[n]]+F[:P]*F[:L]*rand(Normal(0,1),d)=#
	for i=(n-1):-1:1
		Σ=M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],A[t[i]]*M[t[i]])
		Σ=0.5*(Σ+Σ')
		E,V=eig(Σ)
		x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+V*Diagonal(sqrt(map(x->max(x,0.0),E)))*rand(Normal(0,1),d)
		#=F=cholfact(Σ+eye(d),:L,Val{true}) =#
		#=x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+F[:P]*F[:L]*rand(Normal(0,1),d)=#

	end
	Q=innovation(t[1]-t₀)
	A=transition(t[1]-t₀)
	if(t₀!=-Inf)
		x[t₀]=Cₛ*A'*\(A*Cₛ*A'+Q,x[t[1]])+(SafeMvNormal( Cₛ-Cₛ*A'*\(A*Cₛ*A'+Q,A*Cₛ) ))
	end
	return(x);
end


function FFBS2(y,tp,transition,innovation)
	z=SortedDict(Dict{Float64,Array{Float64,2}}())
	tunique=setdiff(tp,collect(keys(y)))
	for t in tunique
		fore=searchsortedlast(y,t); #This routine returns the semitoken of the last item in the container whose key is less than or equal to t. If no such key, then before-start semitoken is returned. 
		aft=searchsortedfirst(y,t); #This routine returns the semitoken of the first item in the container whose key is greater than or equal to t. If no such key, then past-end semitoken is returned. 
		if(fore==beforestartsemitoken(y)) #No key less than or equal to t
			t₁,x₁=deref((y,advance((y,fore))))
			println("Error in FFBS2")
		elseif(aft==pastendsemitoken(y)) #No key greater than or equal to t
			tₙ,xₙ=deref((y,regress((y,aft))))
			y[t]=transition(t-tₙ)*xₙ+(SafeMvNormal(innovation(t-tₙ)))
		else
			(tl,xl)=deref((y,fore))
			(tr,xr)=deref((y,aft))
			Al=transition(t-tl)
			Ar=transition(tr-t)
			Ql=innovation(t-tl)
			Qr=innovation(tr-t)
			AMAQ=Ar*Ql*Ar'+Qr
			AMAQi=SafeInv(AMAQ);
			y[t]=Al*xl+Ql*Ar'*AMAQi*(xr-Ar*Al*xl)+(SafeMvNormal(Ql-Ql*Ar'*AMAQi*(Ar*Ql)))
		end
	end
	for t in tp
		z[t]=y[t]
	end
	for t in tunique
		delete!(y,t)
	end
	return z;
end

function sslogdensity(y,gᵧ,μ,σ²,Cₛ,transition,innovation)
	n=length(y)
	t=collect(keys(y))
	if(gᵧ==1)
		m=Dict{Float64,Array{Float64,2}}()
		sizehint!(m,n)
		M=Dict{Float64,Array{Float64,2}}()
		sizehint!(M,n)
		AMAQ=Dict{Float64,Array{Float64,2}}()
		sizehint!(AMAQ,n)
		A=Dict{Float64,Array{Float64,2}}()
		sizehint!(A,n)
		m[t[1]]=reshape(Cₛ[:,1]*(y[t[1]]-μ)/(σ²+Cₛ[1,1]),d,1)
		M[t[1]]=Cₛ-Cₛ[:,1]*Cₛ[1,:]/(σ²+Cₛ[1,1])
		for i=2:n
			Δ=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ)
			Q=innovation(Δ)
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
			M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
		end
		logdensity=logpdf(Normal(μ,sqrt(σ²+Cₛ[1,1])),y[t[1]])
		for i=2:n
			logdensity=logdensity+logpdf(Normal(μ+(A[t[i-1]][1,:]*m[t[i-1]])[1],sqrt(σ²+AMAQ[t[i-1]][1,1])),y[t[i]])
		end
		return(logdensity)
	else
		logdensity=0
		for i=1:n
			logdensity=logdensity+logpdf(Normal(μ,sqrt(σ²)),y[t[i]])
		end
		return(logdensity)
	end
end

function rho(g,gᵧ,ł)
	n=length(g)
	t=Dict(zip(1:n,sort(collect(keys(g)))))
	if(gᵧ==1)
		Cₛ=statcorr(ł)
		ρ²shape=0
		ρ²rate=(g[t[1]]'*\(Cₛ,g[t[1]]))[1]
		for i=2:n
			Δ=t[i]-t[i-1]
			Q=innovation(Δ,ł)
			A=transition(Δ,ł)
			res=(g[t[i]]-A*g[t[i-1]])
			E,V=eig(Q)
			subind=zeros(Int64,d)
			for w=1:length(E)
				#=if(E[w]<sqrt(eps(real(float(one(eltype(Q)))))))=#
					if(E[w]<0)
						E[w]=Inf
					else
						subind[w]=1
					end
				end
				subindices=find(x->x==1,subind)
				ρ²shape=ρ²shape+length(subindices)
				ρ²rate=ρ²rate+(res'*V[:,subindices]*Diagonal(1./E[subindices])*V[:,subindices]'*res)[1]
			end
			newρ²=rand(InverseGamma(0.5*(ρ²shape+r),0.5*(ρ²rate+r)))
			return(newρ²)
		else
			newρ²=rand(InverseGamma(r/2,r/2))
			return(newρ²)
		end
	end


	function glogdensity(g,gᵧ,ł,ρ²)
		if(ł<0)
			return(-Inf)
		end
		if(ρ²<0)
			return(-Inf)
		end
		n=length(g)
		t=Dict{Float64,Float64}(zip(1:n,sort(collect(keys(g)))))
		if(gᵧ==1)
			Cₛ=statcorr(ł)
			logdensity=-0.5*logdet(Cₛ)[1]-0.5*(g[t[1]]'*\(Cₛ,g[t[1]]))[1]
			for i=2:n
				Δ=t[i]-t[i-1]
				Q=PDMat(innovation(Δ,ł))
				A=transition(Δ,ł)
				res=(g[t[i]]-A*g[t[i-1]])
				logdensity=logdensity-0.5*logdet(ρ²*Q)[1]-0.5*res'*\(ρ²*Q,res)
				#=E,V=eig(ρ²*Q)=#
				#=subind=zeros(Int64,d)=#
				#=for w=1:length(E)=#
					#=[>if(E[w]<sqrt(eps(real(float(one(eltype(ρ²*Q)))))))<]=#
						#=if(E[w]<0)=#
							#=E[w]=Inf=#
						#=else=#
							#=subind[w]=1=#
						#=end=#
					#=end=#
					#=subindices=find(x->x==1,subind)=#
					#=logdensity=logdensity-0.5*sum(log(E[subindices]))-0.5*(res'*V[:,subindices]*Diagonal(1./E[subindices])*V[:,subindices]'*res)[1]=#
				end
				return(logdensity+logpdf(Gamma(2,2),ł))
			else
				return(logpdf(Gamma(2,2),ł))
			end
		end

		#=function sslogdensity(trial::ABtrial,σ²,ł,ρ²)=#
			#=[>(id,Tobs,μg,y₀,y₁,yg,ξ₀ₐᵣ,ξ₀ᵣᵣ,ξ₁ₐᵣ,ξ₁ᵣᵣ,ξ₀ₐₐ,ξ₁ₐₐ,𝑇,g,gᵧ)=params(trial);<]=#
			#=println(trial.gᵧ)=#
			#=println((trial.gᵧ==1))=#
			#=println(trial.μg)=#
			#=println(length(trial.yg))=#
			#=if(trial.gᵧ==1)=#
				#=return sslogdensity(trial.yg,trial.gᵧ,trial.μg,σ²,ł,ρ²)=#
			#=else=#
				#=return 0=#
			#=end=#
		#=end=#
		function sslogdensity(trial::ABtrial,σ²,ł,ρ²)
			(id,Tobs,μg,y₀,y₁,yg,ξ₀ₐᵣ,ξ₀ᵣᵣ,ξ₁ₐᵣ,ξ₁ᵣᵣ,ξ₀ₐₐ,ξ₁ₐₐ,𝑇,g,gᵧ)=params(trial);
			if(gᵧ==1)
				return sslogdensity(yg,gᵧ,μg,σ²,ł,ρ²)
			else
				return 0
			end
		end

