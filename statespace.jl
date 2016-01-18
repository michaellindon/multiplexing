signM=[1.0 ; -1.0  ; 1.0 ; -1.0]
signM=signM*signM'

function statcov(λ,q)
	M=[[(5*q)/(32*λ^7)	,0	,-(q/(32*λ^5))	,0]';
	[0,	q/(32*λ^5),	0,	-(q/(32*λ^3))]';
	[-(q/(32*λ^5)),	0,	q/(32*λ^3),	0]';
	[0,	-(q/(32*λ^3)),	0,	(5*q)/(32*λ)]']
	return(M)
end

function innovation(Δ,λ,q)
	em2Δλ=exp(-2*Δ*λ)
	e2Δλ=exp(2*Δ*λ)
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

function transition(Δ,λ)
	emΔλ=exp(-Δ*λ)
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
		Ap=transition(t[end]-t[prevo],λ)
		Qp=innovation(t[end]-t[prevo],λ,q)
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
				Ap=transition(t[i]-t[prevo],λ)
				Qp=innovation(t[i]-t[prevo],λ,q)
				An=transition(t[i+1]-t[i],λ)
				Qn=innovation(t[i+1]-t[i],λ,q)
				x[t[i]]=Ap*x[t[prevo]]+Qp*An'*\(An*Qp*An'+Qn+0.0000001*eye(d),x[t[i+1]]-An*Ap*x[t[prevo]])+rand(MvNormal(Qp-Qp*An'*\(An*Qp*An'+Qn+0.00000001*eye(d),An*Qp)+0.000001*eye(d)))
			else
				An=transition(t[i+1]-t[i],λ)
				Qn=innovation(t[i+1]-t[i],λ,q)
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

function FFBS(y::Dict,tkeep,μ,β,σ²,λ,q)
	t=union(keys(y),tkeep)
	n=length(t)
	t=Dict(zip(collect(1:n),sort(t)))
	t[0]=t[1]-(t[2]-t[1])
	d=p+1
	L=zeros(d,1)
	L[d,1]=1
	m=Dict(t[0]=>zeros(d,1))
	M=Dict(t[0]=>statcov(λ,q))
	AMAQ=Dict()
	Δ=Dict()
	Q=Dict()
	A=Dict()
	for i=1:n
		Δ[i]=t[i]-t[i-1]
		A[t[i-1]]=transition(Δ[i],λ)
		Q[t[i-1]]=innovation(Δ[i],λ,q)
		if(haskey(y,t[i]))
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+β*AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-β*A[t[i-1]][1,:]*m[t[i-1]])/(σ²+β*AMAQ[t[i-1]][1,1]*β)
			M[t[i]]=AMAQ[t[i-1]]-β*(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])*β/(σ²+β*AMAQ[t[i-1]][1,1]*β)
		else
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]
			M[t[i]]=AMAQ[t[i-1]]
		end
	end
	x=Dict()
	#Backward Sampling
	x[t[n]]=m[t[n]]+rand(MvNormal(M[t[n]]+0.00000000001*eye(d)),1)
	for i=(n-1):-1:0
		x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+rand(MvNormal(M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],A[t[i]]*M[t[i]])+0.0001(eye(d))))
	end
	#=for key in keys(x)=#slow
		#=if(!(key in tkeep))=#
			#=delete!(x,key)=#
		#=end=#
	#=end=#
	xreturn=Dict()
	for time in tkeep
		xreturn[time]=x[time]
	end
	return(xreturn);
end

function FFBS2(y::Dict,tkeep,λ,q)
	σ²=0.000000000000001
	t=union(keys(y),tkeep)
	n=length(t)
	t=Dict(zip(collect(1:n),sort(t)))
	t[0]=t[1]-(t[2]-t[1])
	d=p+1
	L=zeros(d,1)
	L[d,1]=1
	m=Dict(t[0]=>zeros(d,1))
	M=Dict(t[0]=>statcov(λ,q))
	AMAQ=Dict()
	Δ=Dict()
	A=Dict()
	Q=Dict()
	for i=1:n
		Δ[i]=t[i]-t[i-1]
		A[t[i-1]]=transition(Δ[i],λ)
		Q[t[i-1]]=innovation(Δ[i],λ,q)
		if(haskey(y,t[i]))
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
			m[t[i]]=y[t[i]]
			M[t[i]]=zeros(d,d)
		else
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]
			M[t[i]]=AMAQ[t[i-1]]
		end
	end
	x=Dict()
	#Backward Sampling
	if(!haskey(y,t[n]))
		x[t[n]]=m[t[n]]+rand(MvNormal(M[t[n]]+0.000001*eye(d)),1)
	else
		x[t[n]]=y[t[n]]
	end
	for i=(n-1):-1:0
		if(!haskey(y,t[i]))
			x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]]+0.00001*eye(d),x[t[i+1]]-A[t[i]]*m[t[i]])+rand(MvNormal(M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]]+0.00001*eye(d),A[t[i]]*M[t[i]])+0.00001(eye(d))))
		else
			x[t[i]]=y[t[i]]
		end
	end
	#=for key in keys(x)=#slow
		#=if(!(key in tkeep))=#
			#=delete!(x,key)=#
		#=end=#
	#=end=#
	xreturn=Dict()
	for time in tkeep
		xreturn[time]=x[time]
	end
	return(xreturn);
end

function sslogdensity(y,gᵧ,μ,β,σ²,λ,q)
	if(λ<0)
		return(-Inf)
	end
	if(q<0)
		return(-Inf)
	end
	n=length(y)
	t=Dict(zip(1:n,sort(collect(keys(y)))))
	if(gᵧ==1)
		t[0]=t[1]-(t[2]-t[1])
		m=Dict(t[0]=>zeros(d,1))
		M=Dict(t[0]=>statcov(λ,q))
		AMAQ=Dict()
		Δ=Dict()
		Q=Dict()
		A=Dict()
		for i=1:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],λ)
			Q[t[i-1]]=innovation(Δ[i],λ,q)
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+β*AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-β*A[t[i-1]][1,:]*m[t[i-1]])/(σ²+β*AMAQ[t[i-1]][1,1]*β)
			M[t[i]]=AMAQ[t[i-1]]-β*(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])*β/(σ²+β*AMAQ[t[i-1]][1,1]*β)
		end
		logdensity=0
		for i=1:n
			logdensity=logdensity+logpdf(Normal(μ+β*(A[t[i-1]][1,:]*m[t[i-1]])[1],sqrt(σ²+β*AMAQ[t[i-1]][1,1]*β)),y[t[i]])
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

function sslogdensity2(y,gᵧ,μ,β,σ²,λ,q)
	if(λ<0)
		return(-Inf)
	end
	if(q<0)
		return(-Inf)
	end
	n=length(y)
	t=Dict(zip(1:n,sort(collect(keys(y)))))
	if(gᵧ==1)
		Mₛ=statcov(λ,q)
		m=Dict(t[1]=>reshape(Mₛ[:,1]*(y[t[1]]-μ)/(σ²+Mₛ[1,1]),4,1))
		M=Dict(t[1]=>Mₛ-Mₛ[:,1]*Mₛ[1,:]/(σ²+Mₛ[1,1]))
		AMAQ=Dict()
		Δ=Dict()
		Q=Dict()
		A=Dict()
		for i=2:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],λ)
			Q[t[i-1]]=innovation(Δ[i],λ,q)
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+β*AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-β*A[t[i-1]][1,:]*m[t[i-1]])/(σ²+β*AMAQ[t[i-1]][1,1]*β)
			M[t[i]]=AMAQ[t[i-1]]-β*(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])*β/(σ²+β*AMAQ[t[i-1]][1,1]*β)
		end
		logdensity=logpdf(Normal(μ,sqrt(σ²+Mₛ[1,1])),y[t[1]])
		for i=2:n
			logdensity=logdensity+logpdf(Normal(μ+β*(A[t[i-1]][1,:]*m[t[i-1]])[1],sqrt(σ²+β*AMAQ[t[i-1]][1,1]*β)),y[t[i]])
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


function mu(y,gᵧ,σ²ₘ,σ²,λ,q)
	n=length(y)
	t=Dict(zip(1:n,sort(collect(keys(y)))))
	if(gᵧ==1)
		Mₛ=statcov(λ,q)
		m=Dict(t[1]=>reshape(Mₛ[:,1]*(y[t[1]]-μ)/(σ²+Mₛ[1,1]),4,1))
		M=Dict(t[1]=>Mₛ-Mₛ[:,1]*Mₛ[1,:]/(σ²+Mₛ[1,1]))
		AMAQ=Dict()
		Δ=Dict()
		Q=Dict()
		A=Dict()
		for i=2:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],λ)
			Q[t[i-1]]=innovation(Δ[i],λ,q)
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+Q[t[i-1]]
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
			M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
		end
		muprec=1/σ²ₘ
		muprec=muprec+1/(σ²+Mₛ[1,1])
		for i=2:n
			muprec=muprec+1/(σ²+AMAQ[t[i-1]][1,1])
		end
		mumean=0/σ²ₘ
		mumean=mumean+y[t[1]]/(σ²+Mₛ[1,1])
		for i=2:n
			mumean=mumean+(y[t[i]]-(A[t[i-1]]*m[t[i-1]])[1,1])/(σ²+AMAQ[t[i-1]][1,1])
		end
		return(rand(Normal(mumean/muprec,sqrt(1/muprec))))
	else
		return(rand(Normal(n*mean([y[key] for key in keys(y)])/(n+1/σ²ₘ),sqrt(1/(n+1/σ²ₘ)))))
	end
end


function rho(ρ²,g,gᵧ,λ,q)
	n=length(g)
	t=Dict(zip(1:n,sort(collect(keys(g)))))
	if(gᵧ==1)
		Mₛ=statcov(λ,q)
		Δ=Dict()
		Q=Dict()
		A=Dict()
		for i=2:n
			Δ[i]=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ[i],λ)
			Q[t[i-1]]=innovation(Δ[i],λ,q)
		end
		ρ²shape=0
		ρ²rate=(g[t[1]]'*\(Mₛ/ρ²,g[t[1]]))[1]
		for i=2:n
			res=(g[t[i]]-A[t[i-1]]*g[t[i-1]])
			E,V=eig(Q[t[i-1]]/ρ²)
			subind=zeros(Int64,4)
			for w=1:length(E)
				if(E[w]<sqrt(eps(real(float(one(eltype(Q[t[i-1]]/ρ²)))))))
					E[w]=Inf
				else
					subind[w]=1
				end
			end
			subindices=find(x->x==1,subind)
			ρ²shape=ρ²shape+length(subindices)
			ρ²rate=ρ²rate+(res'*V[:,subindices]*Diagonal(1./E[subindices])*V[:,subindices]'*res)[1]
			#=ρ²rate=ρ²rate+(res'*\(Q[t[i-1]]/ρ²,res))[1]=#
			#=ρ²rate=ρ²rate+(res'*pinv(Q[t[i-1]]/ρ²)*res)[1]=#
			#=ρ²rate=ρ²rate+(res'*pinv(Q[t[i-1]]/ρ²,sqrt(eps(real(float(one(eltype(Q[t[i-1]])))))))*res)[1]=#
		end
		#=println((n*d+r)/2)=#
		#=println(2/(ρ²rate+r))=#
		return(rand(InverseGamma(0.5*(ρ²shape+r),0.5*(ρ²rate+r))))
	else
		return(rand(InverseGamma(r/2,r/2)))
	end
end

