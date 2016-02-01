function plot(data::Dict{Float64,Array{Float64,2}})
	plot(sort(collect(keys(data))),[data[key][1] for key in sort(collect(keys(data)))])
end

signM=[1.0 ; -1.0  ; 1.0 ; -1.0]
signM=signM*signM'

function realization(ł,ρ²,t)
	t=sort(t)
	x=Dict{Float64,Array{Float64,2}}()
	Cₛ=statcorr(ł)
	x[t[1]]=rand(MvNormal(ρ²*Cₛ),1)
	#=E,V=eig(ρ²*Cₛ)=#
	#=x[t[1]]=V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d,1)=#
	for i=2:length(t)
		Δ=t[i]-t[i-1]
		A=transition(Δ,ł)
		Q=innovation(Δ,ł)
		x[t[i]]=A*x[t[i-1]]+rand(MvNormal(ρ²*Q),1)
		#=E,V=eig(ρ²*Q)=#
		#=x[t[i]]=A*x[t[i-1]]+V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d,1)=#
	end
	return(x)
end

function statcorr(ł)
	λ=sqrt(2*ν)/ł
	if(d==4)
		q=(96*λ^7)/15
		M=[[(5*q)/(32*λ^7)	,0	,-(q/(32*λ^5))	,0]';
		[0,	q/(32*λ^5),	0,	-(q/(32*λ^3))]';
		[-(q/(32*λ^5)),	0,	q/(32*λ^3),	0]';
		[0,	-(q/(32*λ^3)),	0,	(5*q)/(32*λ)]']
		return(M)
	end
	if(d==3)
		q=(16*λ^5)/3
		M=[[(3*q)/(16*λ^5),0,-(q/(16*λ^3))]';
		[0,q/(16*λ^3),0]';
		[-(q/(16*λ^3)),0,(3*q)/(16*λ)]']
		return(M)
	end
	if(d==2)
		q=4*λ^3
		M=[[q/(4*λ^3),0]';[0,q/(4*λ)]';]
		return(M)
	end
end

function innovation(Δ,ł)
	λ=sqrt(2*ν)/ł
	em2Δλ=exp(-2*Δ*λ)
	e2Δλ=exp(2*Δ*λ)
	if(d==4)
		q=(96*λ^7)/15
		Q=[
		[(q*(45+em2Δλ*(-45-2*Δ*λ*(45+Δ*λ*(45+Δ*λ*(30+Δ*λ*(15+2*Δ*λ*(3+Δ*λ))))))))/(288*λ^7),
		1/72*em2Δλ*q*Δ^6,
		-((em2Δλ*q*(-9+9*e2Δλ+2*Δ*λ*(-9+Δ*λ*(-9+Δ*λ*(-6+Δ*λ*(-3+2*Δ*λ*(-3+Δ*λ)))))))/(288*λ^5)),
		1/72*em2Δλ*q*Δ^4*(3+Δ*λ*(-6+Δ*λ))]';
		[1/72*em2Δλ*q*Δ^6,
		(q*(9+em2Δλ*(-9-2*Δ*λ*(9+Δ*λ*(9+Δ*λ*(6+Δ*λ*(3+2*Δ*λ*(-3+Δ*λ))))))))/(288*λ^5),
		1/72*em2Δλ*q*Δ^4*(-3+Δ*λ)^2,
		(q*(-9+em2Δλ*(9+2*Δ*λ*(9+Δ*λ*(9+Δ*λ*(30+Δ*λ*(-45-2*Δ*λ*(-9+Δ*λ))))))))/(288*λ^3)]';
		[-((em2Δλ*q*(-9+9*e2Δλ+2*Δ*λ*(-9+Δ*λ*(-9+Δ*λ*(-6+Δ*λ*(-3+2*Δ*λ*(-3+Δ*λ)))))))/(288*λ^5)),
		1/72*em2Δλ*q*Δ^4*(-3+Δ*λ)^2,
		(q*(9+em2Δλ*(-9-2*Δ*λ*(9+Δ*λ*(9+Δ*λ*(-42+Δ*λ*(51+2*Δ*λ*(-9+Δ*λ))))))))/(288*λ^3),
		1/72*em2Δλ*q*Δ^2*(6+Δ*λ*(-6+Δ*λ))^2]';
		[1/72*em2Δλ*q*Δ^4*(3+Δ*λ*(-6+Δ*λ)),
		(q*(-9+em2Δλ*(9+2*Δ*λ*(9+Δ*λ*(9+Δ*λ*(30+Δ*λ*(-45-2*Δ*λ*(-9+Δ*λ))))))))/(288*λ^3),
		1/72*em2Δλ*q*Δ^2*(6+Δ*λ*(-6+Δ*λ))^2,
		(q*(45+em2Δλ*(-45+2*Δ*λ*(99+Δ*λ*(-333+Δ*λ*(354+Δ*λ*(-159-2*Δ*λ*(-15+Δ*λ))))))))/(288*λ)]';
		]
		return Q
	end
	if(d==3)
		q=(16*λ^5)/3
		Q=[
		[(q*(3+em2Δλ*(-3-2*Δ*λ*(3+Δ*λ*(3+Δ*λ*(2+Δ*λ))))))/(16*λ^5),
		1/8*em2Δλ*q*Δ^4,
		-((em2Δλ*q*(-1+e2Δλ+2*Δ*λ*(-1+Δ*λ*(-1+Δ*λ*(-2+Δ*λ)))))/(16*λ^3))]';
		[1/8*em2Δλ*q*Δ^4,
		(q*(1+em2Δλ*(-1-2*Δ*λ*(1+Δ*λ*(-1+Δ*λ)^2))))/(16*λ^3),
		1/8*em2Δλ*q*Δ^2*(-2+Δ*λ)^2]';
		[-((em2Δλ*q*(-1+e2Δλ+2*Δ*λ*(-1+Δ*λ*(-1+Δ*λ*(-2+Δ*λ)))))/(16*λ^3)),
		1/8*em2Δλ*q*Δ^2*(-2+Δ*λ)^2,
		(q*(3+em2Δλ*(-3-2*Δ*λ*(-5+Δ*λ*(11+Δ*λ*(-6+Δ*λ))))))/(16*λ)]']
		return Q
	end
	if(d==2)
		q=4*λ^3
		Q=[
		[(q*(1+em2Δλ*(-1-2*Δ*λ*(1+Δ*λ))))/(4*λ^3),
		1/2*em2Δλ*q*Δ^2]';
		[1/2*em2Δλ*q*Δ^2,
		(em2Δλ*q*(-1+e2Δλ-2*Δ*λ*(-1+Δ*λ)))/(4*λ)]';]
		return Q
	end
end

function transition(Δ,ł)
	λ=sqrt(2*ν)/ł
	emΔλ=exp(-Δ*λ)
	if(d==4)
		Φ=[
		[1/6*emΔλ*(6+6*Δ*λ+3*Δ^2*λ^2+Δ^3*λ^3),
		1/2*emΔλ*Δ*(2+2*Δ*λ+Δ^2*λ^2),
		1/2*emΔλ*Δ^2*(1+Δ*λ),
		1/6*emΔλ*Δ^3]';
		[-(1/6)*emΔλ*Δ^3*λ^4,
		-(1/2)*emΔλ*(-2-2*Δ*λ-Δ^2*λ^2+Δ^3*λ^3),
		-(1/2)*emΔλ*Δ*(-2-2*Δ*λ+Δ^2*λ^2),
		-(1/6)*emΔλ*Δ^2*(-3+Δ*λ)]';
		[1/6*emΔλ*Δ^2*λ^4*(-3+Δ*λ),
		1/2*emΔλ*Δ^2*λ^3*(-4+Δ*λ),
		1/2*emΔλ*(2+2*Δ*λ-5*Δ^2*λ^2+Δ^3*λ^3),
		1/6*emΔλ*Δ*(6-6*Δ*λ+Δ^2*λ^2)]';
		[-(1/6)*emΔλ*Δ*λ^4*(6-6*Δ*λ+Δ^2*λ^2),
		-(1/2)*emΔλ*Δ*λ^3*(8-7*Δ*λ+Δ^2*λ^2),
		-(1/2)*emΔλ*Δ*λ^2*(12-8*Δ*λ+Δ^2*λ^2),
		-(1/6)*emΔλ*(-6+18*Δ*λ-9*Δ^2*λ^2+Δ^3*λ^3)]';
		]
		return Φ
	end
	if(d==3)
		Φ=[
		[1/2*emΔλ*(2+2*Δ*λ+Δ^2*λ^2),
		emΔλ*Δ*(1+Δ*λ),
		1/2*emΔλ*Δ^2]';
		[-(1/2)*emΔλ*Δ^2*λ^3,-emΔλ*(-1-Δ*λ+Δ^2*λ^2),
		-(1/2)*emΔλ*Δ*(-2+Δ*λ)]';
		[1/2*emΔλ*Δ*λ^3*(-2+Δ*λ),
		emΔλ*Δ*λ^2*(-3+Δ*λ),
		1/2*emΔλ*(2-4*Δ*λ+Δ^2*λ^2)]';
		]
		return Φ
	end
	if(d==2)
		Φ=[
		[emΔλ*(1+λ*Δ),
		emΔλ*Δ]';
		[-emΔλ*λ^2*Δ,
		-emΔλ*(-1+λ*Δ)]';
		]
		return Φ
	end
end

function FFBS(y::Dict{Float64,Float64},tkeep,μ,σ²,ł,ρ²)
	t=union(keys(y),tkeep)
	n=length(t)
	t=Dict{Int64,Float64}(zip(collect(1:n),sort(t)))
	t[0]=t[1]-(t[2]-t[1])
	d=p+1
	L=zeros(d,1)
	L[d,1]=1
	m=Dict{Float64,Array{Float64,2}}(t[0]=>zeros(d,1))
	M=Dict{Float64,Array{Float64,2}}(t[0]=>ρ²*statcorr(ł))
	AMAQ=Dict{Float64,Array{Float64,2}}()
	Δ=Dict{Int64,Float64}()
	Q=Dict{Float64,Array{Float64,2}}()
	A=Dict{Float64,Array{Float64,2}}()
	for i=1:n
		Δ[i]=t[i]-t[i-1]
		A[t[i-1]]=transition(Δ[i],ł)
		Q[t[i-1]]=innovation(Δ[i],ł)
		if(haskey(y,t[i]))
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
			M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
		else
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]
			M[t[i]]=AMAQ[t[i-1]]
		end
	end
	x=Dict{Float64,Array{Float64,2}}()
	#Backward Sampling
	Σ=M[t[n]]
	Σ=0.5*(Σ+Σ')

	E,V=eig(Σ)
	x[t[n]]=m[t[n]]+V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d)

	#=L=chol(Σ,Val{:L})=#
	#=x[t[n]]=m[t[n]]+L*randn(d,1)=#
	for i=(n-1):-1:0

		Σ=M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],A[t[i]]*M[t[i]])
		Σ=0.5*(Σ+Σ')
		E,V=eig(Σ)
		x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+V*Diagonal(sqrt(map(x->max(x,0.0),E)))*rand(Normal(0,1),d)

		#=L=chol(Σ,Val{:L})=#
		#=x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+L*randn(d,1)=#

		#=x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+rand(MvNormal(M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],A[t[i]]*M[t[i]])))=#
	end
	#=for key in keys(x)=#slow
		#=if(!(key in tkeep))=#
			#=delete!(x,key)=#
		#=end=#
	#=end=#
	xreturn=Dict{Float64,Array{Float64,2}}()
	for time in tkeep
		xreturn[time]=x[time]
	end
	return(xreturn);
end

function FFBS2(y::Dict{Float64,Array{Float64,2}},tkeep,ł,ρ²)
	σ²=0.000000000000001
	t=union(keys(y),tkeep)
	n=length(t)
	t=Dict{Int64,Float64}(zip(collect(1:n),sort(t)))
	t[0]=t[1]-(t[2]-t[1])
	d=p+1
	L=zeros(d,1)
	L[d,1]=1
	m=Dict{Float64,Array{Float64,2}}(t[0]=>zeros(d,1))
	M=Dict{Float64,Array{Float64,2}}(t[0]=>ρ²*statcorr(ł))
	AMAQ=Dict{Float64,Array{Float64,2}}()
	Δ=Dict{Int64,Float64}()
	sizehint!(Δ,n)
	Q=Dict{Float64,Array{Float64,2}}()
	sizehint!(Q,n)
	A=Dict{Float64,Array{Float64,2}}()
	sizehint!(A,n)
	for i=1:n
		Δ[i]=t[i]-t[i-1]
		A[t[i-1]]=transition(Δ[i],ł)
		Q[t[i-1]]=innovation(Δ[i],ł)
		if(haskey(y,t[i]))
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q[t[i-1]]
			m[t[i]]=y[t[i]]
			M[t[i]]=zeros(d,d)
		else
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]
			M[t[i]]=AMAQ[t[i-1]]
		end
	end
	x=Dict{Float64,Array{Float64,2}}()
	#Backward Sampling
	if(!haskey(y,t[n]))
		#=L=chol(M[t[n]],Val{:L})=#
		#=x[t[n]]=m[t[n]]+L*randn(d,1)=#
		#=x[t[n]]=m[t[n]]+rand(MvNormal(M[t[n]]),1)=#
		Σ=M[t[n]]
		Σ=0.5*(Σ+Σ')
		#=L=chol(Σ,Val{:L})=#
		#=x[t[n]]=m[t[n]]+L*randn(d,1)=#
		E,V=eig(Σ)
		x[t[n]]=m[t[n]]+V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d)
	else
		x[t[n]]=y[t[n]]
	end
	for i=(n-1):-1:0
		if(!haskey(y,t[i]))
			#=x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+rand(MvNormal(M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],A[t[i]]*M[t[i]])))=#
			Σ=M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],A[t[i]]*M[t[i]])
			Σ=0.5*(Σ+Σ')
			#=L=chol(Σ,Val{:L})=#
			#=x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+L*randn(d,1)=#
			E,V=eig(Σ)
			x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+V*Diagonal(sqrt(map(x->max(x,0.0),E)))*rand(Normal(0,1),d)
		else
			x[t[i]]=y[t[i]]
		end
	end
	#=for key in keys(x)=#slow
		#=if(!(key in tkeep))=#
			#=delete!(x,key)=#
		#=end=#
	#=end=#
	xreturn=Dict{Float64,Array{Float64,2}}()
	for time in tkeep
		xreturn[time]=x[time]
	end
	return(xreturn);
end



function sslogdensity(y,gᵧ,μ,σ²,ł,ρ²)
	if(ł<0)
		return(-Inf)
	end
	if(ρ²<0)
		return(-Inf)
	end
	n=length(y)
	t=Dict(zip(1:n,sort(collect(keys(y)))))
	if(gᵧ==1)
		Cₛ=statcorr(ł)
		m=Dict(t[1]=>reshape(ρ²*Cₛ[:,1]*(y[t[1]]-μ)/(σ²+ρ²*Cₛ[1,1]),d,1))
		M=Dict(t[1]=>ρ²*Cₛ-ρ²*Cₛ[:,1]*Cₛ[1,:]*ρ²/(σ²+ρ²*Cₛ[1,1]))
		AMAQ=Dict{Float64,Array{Float64,2}}()
		Δ=Dict{Float64,Float64}()
		Q=Dict{Float64,Array{Float64,2}}()
		A=Dict{Float64,Array{Float64,2}}()
		for i=2:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],ł)
			Q[t[i-1]]=innovation(Δ[i],ł)
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
			M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
		end
		logdensity=logpdf(Normal(μ,sqrt(σ²+ρ²*Cₛ[1,1])),y[t[1]])
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
		Δ=Dict{Float64,Float64}()
		sizehint!(Δ,n)
		Q=Dict{Float64,Array{Float64,2}}()
		sizehint!(Q,n)
		A=Dict{Float64,Array{Float64,2}}()
		sizehint!(A,n)
		for i=2:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],ł)
			Q[t[i-1]]=innovation(Δ[i],ł)
		end
		ρ²shape=0
		ρ²rate=(g[t[1]]'*\(Cₛ,g[t[1]]))[1]
		for i=2:n
			res=(g[t[i]]-A[t[i-1]]*g[t[i-1]])
			E,V=eig(Q[t[i-1]])
			subind=zeros(Int64,d)
			for w=1:length(E)
				if(E[w]<sqrt(eps(real(float(one(eltype(Q[t[i-1]])))))))
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
		Δ=Dict{Float64,Float64}()
		sizehint!(Δ,n)
		Q=Dict{Float64,Array{Float64,2}}()
		sizehint!(Q,n)
		A=Dict{Float64,Array{Float64,2}}()
		sizehint!(A,n)
		for i=2:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],ł)
			Q[t[i-1]]=innovation(Δ[i],ł)
		end
		logdensity=-0.5*logdet(ρ²*Cₛ)-0.5*(g[t[1]]'*\(ρ²*Cₛ,g[t[1]]))[1]
		for i=2:n
			res=(g[t[i]]-A[t[i-1]]*g[t[i-1]])
			E,V=eig(ρ²*Q[t[i-1]])
			subind=zeros(Int64,d)
			for w=1:length(E)
				#=if(E[w]<sqrt(eps(real(float(one(eltype(ρ²*Q[t[i-1]])))))))=#
				if(E[w]<0)
					println(i)
					E[w]=Inf
				else
					subind[w]=1
				end
			end
			subindices=find(x->x==1,subind)
			logdensity=logdensity-0.5*sum(log(E[subindices]))-0.5*(res'*V[:,subindices]*Diagonal(1./E[subindices])*V[:,subindices]'*res)[1]
		end
		return(logdensity+logpdf(Gamma(2,2),ł))
	else
		return(logpdf(Gamma(2,2),ł))
	end
end

#Appears FFBS2 faster than this one
function predictFunction(xInput::Dict,tkeep,λ,q)
	x=copy(xInput)
	t=sort(union(keys(x),tkeep))
	n=length(t)
	if(!haskey(x,t[end])) #if x does not have key i.e. unobserved
		prevo=length(t) #Go backwards in time to find latest observed time
		while(!haskey(x,t[prevo])) #while x does not have key
			prevo=prevo-1 #go back
		end #exit when x has the key, this is the latest observed time before the current unobserved time
		Ap=transition(t[end]-t[prevo],ł)
		Qp=innovation(t[end]-t[prevo],ł)
		x[t[end]]=Ap*x[t[prevo]]+rand(MvNormal(Qp+0.0000000001(eye(d))),1)
	end
	for i=(n-1):-1:1
		if(!haskey(x,t[i])) #if x does not have key i.e. unobserved
			prevo=i #Go backwards in time to find latest observed time
			while(!haskey(x,t[prevo])) #while x does not have key
				prevo=prevo-1 #go back
				if(prevo==0)
					break
				end
			end #exit when x has the key, this is the latest observed time before the current unobserved time
			if(prevo!=0)
				Ap=transition(t[i]-t[prevo],ł)
				Qp=innovation(t[i]-t[prevo],ł)
				An=transition(t[i+1]-t[i],ł)
				Qn=innovation(t[i+1]-t[i],ł)
				x[t[i]]=Ap*x[t[prevo]]+Qp*An'*\(An*Qp*An'+Qn+0.0000001*eye(d),x[t[i+1]]-An*Ap*x[t[prevo]])+rand(MvNormal(Qp-Qp*An'*\(An*Qp*An'+Qn+0.00000001*eye(d),An*Qp)+0.000001*eye(d)))
			else
				An=transition(t[i+1]-t[i],ł)
				Qn=innovation(t[i+1]-t[i],ł)
				x[t[i]]=(signM.*An)*x[t[i+1]]+rand(MvNormal((signM.*Qn)+0.0001*eye(d)))
			end
		end
	end
	for key in keys(x)
		if(!(key in tkeep))
			delete!(x,key)
		end
	end
	return(x)
end

function sslogdensity2(y,gᵧ,μ,σ²,ł,ρ²)
	if(ł<0)
		return(-Inf)
	end
	if(ρ²<0)
		return(-Inf)
	end
	n=length(y)
	t=Dict(zip(1:n,sort(collect(keys(y)))))
	if(gᵧ==1)
		t[0]=t[1]-(t[2]-t[1])
		m=Dict(t[0]=>zeros(d,1))
		M=Dict(t[0]=>ρ²*statcorr(ł))
		AMAQ=Dict{Float64,Array{Float64,2}}()
		Δ=Dict{Float64,Float64}()
		Q=Dict{Float64,Array{Float64,2}}()
		A=Dict{Float64,Array{Float64,2}}()
		for i=1:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],ł)
			Q[t[i-1]]=innovation(Δ[i],ł)
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
			M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
		end
		logdensity=0
		for i=1:n
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



function mu(y,gᵧ,σ²ₘ,σ²,ł,ρ²)
	n=length(y)
	t=Dict(zip(1:n,sort(collect(keys(y)))))
	if(gᵧ==1)
		m=Dict{Float64,Array{Float64,2}}()
		sizehint!(m,n)
		M=Dict{Float64,Array{Float64,2}}()
		sizehint!(M,n)
		AMAQ=Dict{Float64,Array{Float64,2}}()
		sizehint!(AMAQ,n)
		Δ=Dict{Int64,Float64}()
		sizehint!(Δ,n)
		Q=Dict{Float64,Array{Float64,2}}()
		sizehint!(Q,n)
		A=Dict{Float64,Array{Float64,2}}()
		sizehint!(A,n)
		Mₛ=ρ²*statcorr(ł)
		m[t[1]]=reshape(Mₛ[:,1]*(y[t[1]])/(σ²+Mₛ[1,1]),d,1)
		M[t[1]]=Mₛ-Mₛ[:,1]*Mₛ[1,:]/(σ²+Mₛ[1,1])
		for i=2:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],ł)
			Q[t[i-1]]=innovation(Δ[i],ł)
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
			M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
		end
		muprec=1/σ²ₘ+1/(σ²+Mₛ[1,1])+reduce(+,[1/(σ²+AMAQ[t[i-1]][1,1]) for i=2:n])
		mumean=0/σ²ₘ+y[t[1]]/(σ²+Mₛ[1,1])+reduce(+,[(y[t[i]]-(A[t[i-1]]*m[t[i-1]])[1,1])/(σ²+AMAQ[t[i-1]][1,1]) for i=2:n])
		muvar=1/muprec
		return(Normal(muvar*mumean,√muvar))
	else
		return(Normal((n/σ²)*mean([y[key] for key in keys(y)])*(1/((n/σ²)+(1/σ²ₘ))),sqrt(1/((n/σ²)+(1/σ²ₘ)))))
	end
end

