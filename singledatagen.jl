srand(3)
const global p=2
const global ν=p+0.5
const global d=p+1
Λ=4000
#Create 0 branch
srand(2)
ρ²₀=9.0
ł₀=0.1

f₀=realization(statcov(ł₀,ρ²₀),transition(ł₀),innovation(ł₀,ρ²₀),collect(0:1/1000:1))
μ₀=0.0
Atrials=Array{Atrial}(1)
for i=1:length(Atrials)
	Tobs=rand(Beta(7,1))
	if(i==1) 
		Tobs=1.0
	end
	Tₚ=rand(PPProcess(Λ),0,Tobs)
	f₀ₚ=FFBS2(f₀,Tₚ,transition(ł₀),innovation(ł₀,ρ²₀))
	ξ₀ₐ=Dict()
	ξ₀ᵣ=Dict()
	y=SortedDict(Dict{Float64,Float64}())
	for t in Tₚ
		yf=rand(Normal(μ₀+f₀ₚ[t][1],1))
		y[t]=yf
		if(yf>0)
			ξ₀ₐ[t]=Dict("yf"=>yf)
		else
			ξ₀ᵣ[t]=Dict("yf"=>yf)
		end
	end
	Atrials[i]=Atrial(i,Tobs,y,ξ₀ₐ,ξ₀ᵣ)
end

function foof(x)
	if (x>=0 && x<0.2)
		return 0.3
	end
	if (x>=0.2 && x<0.3)
		return 0.7
	end
	if (x>=0.3 && x<0.5)
		return 0.9
	end
	if (x>=0.5 && x<0.7)
		return 0.2
	end
	if (x>=0.7 && x<1.0)
		return 0.5
	end
end

fooy=SortedDict(Dict{Float64,Float64}())
fooa=rand(PPProcess(x->Λ*foof(x),Λ),0,1)
foor=rand(PPProcess(x->Λ*(1-foof(x)),Λ),0,1)
fooξ₀ₐ=Dict()
fooξ₀ᵣ=Dict()
for foot in union(fooa,foor)
	fooy[foot]=0.0
end
for foot in fooa
	fooξ₀ₐ[foot]=Dict("yf"=>fooy[foot])
end
for foot in foor
	fooξ₀ᵣ[foot]=Dict("yf"=>fooy[foot])
end
Atrials[1]=Atrial(1,1,fooy,fooξ₀ₐ,fooξ₀ᵣ)




figure()
subplot(211)
#=plot((collect(keys(f₀))),[Λ*Φ(μ₀+f₀[key][1]) for key in (collect(keys(f₀)))])=#
foor=sort(rand(1000))
plot(foor,map(x->Λ*foof(x),foor))
xlim(0,1)
subplot(212)
for trial in Atrials
	plot(collect(keys(trial.ξ₀ₐ)),trial.id*ones(length(trial.ξ₀ₐ)),c="blue",marker="|",linestyle="None",markersize=1)
end
xlim(0,1)




T₀=Array{Float64}(0)
for trial in Atrials
	T₀=union(T₀,collect(keys(trial.ξ₀ₐ)))
end
T₀=sort(T₀)
T=union(T₀)

T=sort(T)
T=convert(Array{Float64},T)
f₀=FFBS2(f₀,T,transition(ł₀),innovation(ł₀,ρ²₀))
