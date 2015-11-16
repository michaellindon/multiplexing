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

function predictFunction(tInput::Dict,tOutput::Array{Float64},λ,q)
	x=copy(tInput)
	times=sort([collect(keys(x));tOutput])
	if(!haskey(x,times[end])) #if x does not have key i.e. unobserved
		prevo=length(times) #Go backwards in time to find latest observed time
		while(!haskey(x,times[prevo])) #while x does not have key
			prevo=prevo-1 #go back
		end #exit when x has the key, this is the latest observed time before the current unobserved time
		Ap=transition(times[end]-times[prevo],λ)
		Qp=innovation(times[end]-times[prevo],λ,q)
		x[times[end]]=Ap*x[times[prevo]]+rand(MvNormal(Qp+0.0000000001(eye(d))),1)
	end
	for i=length(times-1):-1:1
		if(!haskey(x,times[i])) #if x does not have key i.e. unobserved
			prevo=i #Go backwards in time to find latest observed time
			while(!haskey(x,times[prevo])) #while x does not have key
				prevo=prevo-1 #go back
				if(prevo==0)
					break
				end
			end #exit when x has the key, this is the latest observed time before the current unobserved time
			if(prevo!=0)
			Ap=transition(times[i]-times[prevo],λ)
			Qp=innovation(times[i]-times[prevo],λ,q)
			An=transition(times[i+1]-times[i],λ)
			Qn=innovation(times[i+1]-times[i],λ,q)
			x[times[i]]=Ap*x[times[prevo]]+Qp*An'*\(An*Qp*An'+Qn+0.0000001*eye(d),x[times[i+1]]-An*Ap*x[times[prevo]])+rand(MvNormal(Qp-Qp*An'*\(An*Qp*An'+Qn+0.00000001*eye(d),An*Qp)+0.000001*eye(d)))
		else
			An=transition(times[i+1]-times[i],λ)
			Qn=innovation(times[i+1]-times[i],λ,q)
			x[times[i]]=(signM.*An)*x[times[i+1]]+rand(MvNormal((signM.*Qn)+0.0000000001*eye(d)))
		end
		end
	end
	return(x)
end
