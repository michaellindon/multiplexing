function plot(data::Dict{Float64,Array{Float64,2}})
	plot(sort(collect(keys(data))),[data[key][1] for key in sort(collect(keys(data)))])
end

signM=[1.0 ; -1.0  ; 1.0 ; -1.0]
signM=signM*signM'

function realization(ł,ρ²,t)
	t=sort(t)
	x=SortedDict(Dict{Float64,Array{Float64,1}}())
	Cₛ=statcorr(ł)
	x[t[1]]=rand(MvNormal(ρ²*Cₛ))
	#=E,V=eig(ρ²*Cₛ)=#
	#=x[t[1]]=V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d,1)=#
	for i=2:length(t)
		Δ=t[i]-t[i-1]
		A=transition(Δ,ł)
		Q=innovation(Δ,ł)
		x[t[i]]=A*x[t[i-1]]+rand(MvNormal(ρ²*Q))
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
		λ³=λ*λ*λ
		q=4*λ³
		M=Array{Float64}(2,2)
		M[1,1]=q/(4*λ³)
		M[1,2]=0
		M[2,1]=0
		M[2,2]=q/(4*λ)
		#=M=[[q/(4*λ^3),0]';[0,q/(4*λ)]';]=#
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
		Q=Array{Float64}(3,3)
		Q[1,1]=(q*(3+em2Δλ*(-3-2*Δ*λ*(3+Δ*λ*(3+Δ*λ*(2+Δ*λ))))))/(16*λ^5);
		Q[1,2]=1/8*em2Δλ*q*Δ^4;
		Q[1,3]=-((em2Δλ*q*(-1+e2Δλ+2*Δ*λ*(-1+Δ*λ*(-1+Δ*λ*(-2+Δ*λ)))))/(16*λ^3));
		Q[2,1]=1/8*em2Δλ*q*Δ^4
		Q[2,2]=(q*(1+em2Δλ*(-1-2*Δ*λ*(1+Δ*λ*(-1+Δ*λ)^2))))/(16*λ^3)
		Q[2,3]=1/8*em2Δλ*q*Δ^2*(-2+Δ*λ)^2
		Q[3,1]=-((em2Δλ*q*(-1+e2Δλ+2*Δ*λ*(-1+Δ*λ*(-1+Δ*λ*(-2+Δ*λ)))))/(16*λ^3))
		Q[3,2]=1/8*em2Δλ*q*Δ^2*(-2+Δ*λ)^2
		Q[3,3]=(q*(3+em2Δλ*(-3-2*Δ*λ*(-5+Δ*λ*(11+Δ*λ*(-6+Δ*λ))))))/(16*λ)
		return Q
		#=Q=[=#
		#=[(q*(3+em2Δλ*(-3-2*Δ*λ*(3+Δ*λ*(3+Δ*λ*(2+Δ*λ))))))/(16*λ^5),=#
		#=1/8*em2Δλ*q*Δ^4,=#
		#=-((em2Δλ*q*(-1+e2Δλ+2*Δ*λ*(-1+Δ*λ*(-1+Δ*λ*(-2+Δ*λ)))))/(16*λ^3))]';=#
		#=[1/8*em2Δλ*q*Δ^4,=#
		#=(q*(1+em2Δλ*(-1-2*Δ*λ*(1+Δ*λ*(-1+Δ*λ)^2))))/(16*λ^3),=#
		#=1/8*em2Δλ*q*Δ^2*(-2+Δ*λ)^2]';=#
		#=[-((em2Δλ*q*(-1+e2Δλ+2*Δ*λ*(-1+Δ*λ*(-1+Δ*λ*(-2+Δ*λ)))))/(16*λ^3)),=#
		#=1/8*em2Δλ*q*Δ^2*(-2+Δ*λ)^2,=#
		#=(q*(3+em2Δλ*(-3-2*Δ*λ*(-5+Δ*λ*(11+Δ*λ*(-6+Δ*λ))))))/(16*λ)]']=#
		#=return Q=#
	end
	if(d==2)
		λ³=λ*λ*λ
		q=4*λ³
		Q=Array{Float64}(2,2)
		Q[1,1]=(q*(1+em2Δλ*(-1-2*Δ*λ*(1+Δ*λ))))/(4*λ³)
		Q[1,2]=1/2*em2Δλ*q*Δ^2
		Q[2,1]=1/2*em2Δλ*q*Δ^2
		Q[2,2]=(em2Δλ*q*(-1+e2Δλ-2*Δ*λ*(-1+Δ*λ)))/(4*λ)
		#=Q=[=#
		#=[(q*(1+em2Δλ*(-1-2*Δ*λ*(1+Δ*λ))))/(4*λ^3),=#
		#=1/2*em2Δλ*q*Δ^2]';=#
		#=[1/2*em2Δλ*q*Δ^2,=#
		#=(em2Δλ*q*(-1+e2Δλ-2*Δ*λ*(-1+Δ*λ)))/(4*λ)]';]=#
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
		
		#=Φ=[=#
		#=[1/2*emΔλ*(2+2*Δ*λ+Δ^2*λ^2),=#
		#=emΔλ*Δ*(1+Δ*λ),=#
		#=1/2*emΔλ*Δ^2]';=#
		#=[-(1/2)*emΔλ*Δ^2*λ^3,-emΔλ*(-1-Δ*λ+Δ^2*λ^2),=#
		#=-(1/2)*emΔλ*Δ*(-2+Δ*λ)]';=#
		#=[1/2*emΔλ*Δ*λ^3*(-2+Δ*λ),=#
		#=emΔλ*Δ*λ^2*(-3+Δ*λ),=#
		#=1/2*emΔλ*(2-4*Δ*λ+Δ^2*λ^2)]';=#
		#=]=#
		#=return Φ=#
	end
	if(d==2)
		Φ=Array{Float64}(2,2)
		Φ[1,1]=emΔλ*(1+λ*Δ)
		Φ[1,2]=emΔλ*Δ
		Φ[2,1]=-emΔλ*λ^2*Δ
		Φ[2,2]=-emΔλ*(-1+λ*Δ)
		#=Φ=[=#
		#=[emΔλ*(1+λ*Δ),=#
		#=emΔλ*Δ]';=#
		#=[-emΔλ*λ^2*Δ,=#
		#=-emΔλ*(-1+λ*Δ)]';=#
		#=]=#
		return Φ
	end
end

function FFBS(y,μ,σ²,ł,ρ²)
	t=collect(keys(y))
	n=length(t)
	Cₛ=statcorr(ł)
	m=Dict{Float64,Array{Float64,2}}()
	sizehint!(m,n)
	m[t[1]]=reshape(ρ²*Cₛ[:,1]*(y[t[1]]-μ)/(σ²+ρ²*Cₛ[1,1]),d,1)
	M=Dict{Float64,Array{Float64,2}}()
	sizehint!(M,n)
	M[t[1]]=ρ²*Cₛ-ρ²*Cₛ[:,1]*Cₛ[1,:]*ρ²/(σ²+ρ²*Cₛ[1,1])
	AMAQ=Dict{Float64,Array{Float64,2}}()
	sizehint!(AMAQ,n)
	A=Dict{Float64,Array{Float64,2}}()
	sizehint!(A,n)
	for i=2:n
		Δ=t[i]-t[i-1]
		A[t[i-1]]=transition(Δ,ł)
		Q=innovation(Δ,ł)
		if(haskey(y,t[i]))
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q
			m[t[i]]=A[t[i-1]]*m[t[i-1]]+AMAQ[t[i-1]][:,1]*(y[t[i]]-μ-A[t[i-1]][1,:]*m[t[i-1]])/(σ²+AMAQ[t[i-1]][1,1])
			M[t[i]]=AMAQ[t[i-1]]-(AMAQ[t[i-1]][:,1]*AMAQ[t[i-1]][1,:])/(σ²+AMAQ[t[i-1]][1,1])
		else
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q
			m[t[i]]=A[t[i-1]]*m[t[i-1]]
			M[t[i]]=AMAQ[t[i-1]]
		end
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
	return(x);
end

function FFBS2(xc,tp,ł,ρ²)
	tc=collect(keys(xc))
	c=1
	p=1
	i=1
	n=length(tc)+length(tp)
	t=Array{Float64}(n)
	while(c<=length(tc) && p<=length(tp))
		if(tp[p]<tc[c])
			t[i]=tp[p]
			p=p+1;
			i=i+1;
		else
			t[i]=tc[c]
			c=c+1;
			i=i+1;
		end
	end
	while(p<=length(tp))
		t[i]=tp[p]
		i=i+1
		p=p+1
	end
	while(c<=length(tc))
		t[i]=tc[c]
		i=i+1
		c=c+1
	end

	m=Dict{Float64,Array{Float64,2}}()
	sizehint!(m,n)
	if(haskey(xc,t[1]))
		m[t[1]]=xc[t[1]]
	else
		m[t[1]]=zeros(Float64,d,1)
	end
	M=Dict{Float64,Array{Float64,2}}()
	sizehint!(M,n)
	if(haskey(xc,t[1]))
		M[t[1]]=zeros(Float64,d,d)
	else
		M[t[1]]=ρ²*statcorr(ł)
	end
	AMAQ=Dict{Float64,Array{Float64,2}}()
	sizehint!(AMAQ,n)
	A=Dict{Float64,Array{Float64,2}}()
	sizehint!(A,n)
	for i=2:n
		Δ=t[i]-t[i-1]
		A[t[i-1]]=transition(Δ,ł)
		Q=innovation(Δ,ł)
		if(haskey(xc,t[i]))
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q
			m[t[i]]=xc[t[i]]
			M[t[i]]=zeros(Float64,d,d)
		else
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q
			m[t[i]]=A[t[i-1]]*m[t[i-1]]
			M[t[i]]=AMAQ[t[i-1]]
		end
	end
	x=Dict{Float64,Array{Float64,2}}()
	#Backward Sampling
	if(!haskey(xc,t[n]))
		Σ=M[t[n]]
		Σ=0.5*(Σ+Σ')
		E,V=eig(Σ)
		x[t[n]]=m[t[n]]+V*Diagonal(sqrt(map(x->max(x,0),E)))*rand(Normal(0,1),d)
		#=F=cholfact(Σ+eye(d),:L,Val{true}) =#
		#=x[t[n]]=m[t[n]]+F[:P]*F[:L]*rand(Normal(0,1),d)=#
	else
		x[t[n]]=xc[t[n]]
	end
	for i=(n-1):-1:1
		if(!haskey(xc,t[i]))
			Σ=M[t[i]]-M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],A[t[i]]*M[t[i]])
			Σ=0.5*(Σ+Σ')
			E,V=eig(Σ)
			x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+V*Diagonal(sqrt(map(x->max(x,0.0),E)))*rand(Normal(0,1),d)
			#=F=cholfact(Σ+eye(d),:L,Val{true}) =#
			#=x[t[i]]=m[t[i]]+M[t[i]]*A[t[i]]'*\(AMAQ[t[i]],x[t[i+1]]-A[t[i]]*m[t[i]])+F[:P]*F[:L]*rand(Normal(0,1),d)=#
		else
			x[t[i]]=xc[t[i]]
		end
	end
	xp=SortedDict(Dict{Float64,Array{Float64,2}}())
	for Zeit in tp
		xp[Zeit]=x[Zeit]
	end
	return(xp);
end

function sslogdensity(y,gᵧ,μ,σ²,ł,ρ²)
	if(ł<0)
		return(-Inf)
	end
	if(ρ²<0)
		return(-Inf)
	end
	n=length(y)
	t=collect(keys(y))
	if(gᵧ==1)
		Cₛ=statcorr(ł)
		m=Dict{Float64,Array{Float64,2}}()
		sizehint!(m,n)
		M=Dict{Float64,Array{Float64,2}}()
		sizehint!(M,n)
		AMAQ=Dict{Float64,Array{Float64,2}}()
		sizehint!(AMAQ,n)
		A=Dict{Float64,Array{Float64,2}}()
		sizehint!(A,n)
		m[t[1]]=reshape(ρ²*Cₛ[:,1]*(y[t[1]]-μ)/(σ²+ρ²*Cₛ[1,1]),d,1)
		M[t[1]]=ρ²*Cₛ-ρ²*Cₛ[:,1]*Cₛ[1,:]*ρ²/(σ²+ρ²*Cₛ[1,1])
		for i=2:n
			Δ=t[i]-t[i-1]
			A[t[i-1]]=transition(Δ,ł)
			Q=innovation(Δ,ł)
			AMAQ[t[i-1]]=A[t[i-1]]*M[t[i-1]]*A[t[i-1]]'+ρ²*Q
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
		logdensity=-0.5*logdet(ρ²*Cₛ)[1]-0.5*(g[t[1]]'*\(ρ²*Cₛ,g[t[1]]))[1]
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

