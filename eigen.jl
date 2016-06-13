

function ms(y,μ)
	x=deepcopy(y)
	for key in keys(x)
		x[key]=x[key]-μ(key)
	end
	return x
end

function FFBS(y,μ,σ²,ł,ρ²)
	t=collect(keys(y));
	y=collect(values(y))
	if(!all(isfinite(y)))
		println("inputs to FFBS are not finite")
	end
	n=length(t)
	x=zeros(Float64,3*n)
	normals=rand(Normal(0,1),3*n)
	ccall((:FFBS3d, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Int32,Float64,Float64,Float64,Float64,Ref{Cdouble}), y,x,t,n,ł,σ²,ρ²,μ,normals)
	if(!all(isfinite(x)))
		println("outputs to FFBS are not finite")
	end
	x=reshape(x,3,n)
	z=SortedDict(Dict{Float64,Array{Float64,1}}())
	for i=1:n
		z[t[i]]=x[:,i]
	end
	return z
end

function sslogdensity(y,gᵧ,μ,σ²,ł,ρ²)
	t=collect(keys(y));
	y=collect(values(y));
	n=length(t)
	return  ccall((:LogDensity3d, "./eigen.so"), Float64, (Ref{Cdouble},Ref{Cdouble},Int32,Int32,Float64,Float64,Float64,Float64), y,t,gᵧ,n,μ,σ²,ł,ρ²)
end

function mulogdensity(y,gᵧ,σ²,ł,ρ²,σ²ₘ)
	t=collect(keys(y));
	y=collect(values(y));
	n=length(t)
	if(gᵧ==1)
		return  ccall((:mulogdensity, "./eigen.so"), Float64, (Ref{Cdouble},Ref{Cdouble},Int32,Float64,Float64,Float64,Float64), y,t,n,σ²,ł,ρ²,σ²ₘ)
	else
		mm=sum(y)/σ²
		pp=(n/σ²+1/σ²ₘ)
		-0.5*n*log(2*pi*σ²)-0.5*log(2*pi*σ²ₘ)+0.5*log(2*pi/pp)-0.5*(dot(y,y)/σ²-mm*mm/pp)
	end
end

function oldFFBS2(y,tout,ł,ρ²)
	tin=collect(keys(y))
	xin=vcat(collect(values(y))...);	
	if(!all(isfinite(xin)))
		println("Inputs to FFBS2 are not finite")
	end
	xout=zeros(Float64,3*length(tout))
	#=tout=convert(Array{Float64},tout)=#
	ccall((:Predict3d, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Int32,Ref{Cdouble},Ref{Cdouble},Int32, Int32,Float64,Float64), xin,tin,length(tin),xout,tout,length(tout),1,ł,ρ²)
	if(!all(isfinite(xout)))
		println( "outputs to FFBS2 are not finite")
	end
	xout=reshape(xout,3,length(tout))
	z=SortedDict(Dict{Float64,Array{Float64,1}}())
	for i=1:length(tout)
		z[tout[i]]=xout[:,i]
	end
	return z
end

function innovation3d(Δ,ł)
	ccall((:Innovation3dPrint, "./eigen.so"), Void, (Float64,Float64), Δ,ł)
end

function FFBS2(y,tp,łs,ρ²)
	ł=sqrt(5.0)/łs;
	z=SortedDict(Dict{Float64,Array{Float64,1}}())
	for t in setdiff(tp,collect(keys(y)))
		𝑍=rand(Normal(0,1),3)
		xout=zeros(Float64,3)
		fore=searchsortedlast(y,t); #This routine returns the semitoken of the last item in the container whose key is less than or equal to t. If no such key, then before-start semitoken is returned. 
		aft=searchsortedfirst(y,t); #This routine returns the semitoken of the first item in the container whose key is greater than or equal to t. If no such key, then past-end semitoken is returned. 
		if(fore==beforestartsemitoken(y)) #No key less than or equal to t
			t₁,x₁=deref((y,advance((y,fore))))
			ccall((:Backwards, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Float64,Float64,Ref{Cdouble},Float64),xout,𝑍,t,t₁,x₁,ł)
			y[t]=xout;
		elseif(aft==pastendsemitoken(y)) #No key greater than or equal to t
			tₙ,xₙ=deref((y,regress((y,aft))))
			ccall((:Forwards, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Float64,Float64,Ref{Cdouble},Float64),xout,𝑍,t,tₙ,xₙ,ł)
			y[t]=xout;
		else
			(tl,vl)=deref((y,fore))
			(tr,vr)=deref((y,aft))
			ccall((:ThreePoint, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Float64,Float64,Ref{Cdouble},Float64,Ref{Cdouble},Float64),xout,𝑍,t,tl,vl,tr,vr,ł)
			y[t]=xout;
		end
	end
	for t in tp
		z[t]=y[t]
	end
	for t in setdiff(tp,collect(keys(y)))
		delete!(y,t)
	end
	return z;
end

function mu(y,gᵧ,σ²,ł,ρ²,σ²ₘ)
	t=collect(keys(y));
	y=collect(values(y));
	n=length(t)
	if(gᵧ==1)
		μmean=zeros(Float64,1)
		μprec=zeros(Float64,1)
		ccall((:mu3d, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Int32,Float64,Float64,Float64, Ref{Cdouble}, Ref{Cdouble},Float64), y,t,n,σ²,ł,ρ²,μmean,μprec,σ²ₘ)
		return Normal(μmean[1],sqrt(1.0/μprec[1]))
	else
		μprec=1/σ²ₘ+n/σ²
		μmean=(n/σ²)*mean(y)/μprec
		return Normal(μmean,sqrt(1.0/μprec))
	end
end
