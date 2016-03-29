importall Base.Random
importall Distributions
importall PyPlot

include("poissonpointprocess.jl")
include("analyticalEigen.jl")
include("statespace.jl")
include("fulldatageneration2.jl")


niter=400
iter=1
gᵧ=1
σ²ₘ=3
accepted=0
r=2
gᵧTrace=zeros(niter)
ρ²Trace=zeros(niter)
łTrace=zeros(niter)
μTrace=zeros(niter)
logodds=1
logoddsTrace=zeros(niter)
σ²=1.0
for key in keys(f₀)
	f₀[key]=zeros(d,1)
end
for key in keys(f₁)
	f₁[key]=zeros(d,1)
end
for key in keys(g)
	g[key]=zeros(d,1)
end

for iter=1:niter
	println(iter)
	#=println(iter)=#
	Ξₚ=PPProcess(Λ)

	#Sample Realizations
	Ξ₀ₐ(ξ₀ₐ,μ₀,f₀);
	ξ₀ᵣ=Ξ₀ᵣ(μ₀,f₀,ł₀,ρ²₀,Ξₚ,Tobs);
	Ξ₁ₐ(ξ₁ₐ,μ₁,f₁);
	ξ₁ᵣ=Ξ₁ᵣ(μ₁,f₁,ł₁,ρ²₁,Ξₚ,Tobs);
	ξ₀ₐᵣ=Ξ₀ₐᵣ(μg,g,łg,ρ²g,μ₀,f₀,ł₀,ρ²₀,Ξₚ,Tobs);
	ξ₀ᵣᵣ=Ξ₀ᵣᵣ(μ₀,f₀,ł₀,ρ²₀,Ξₚ,Tobs);
	ξ₁ₐᵣ=Ξ₁ₐᵣ(μg,g,łg,ρ²g,μ₁,f₁,ł₁,ρ²₁,Ξₚ,Tobs);
	ξ₁ᵣᵣ=Ξ₁ᵣᵣ(μ₁,f₁,ł₁,ρ²₁,Ξₚ,Tobs);
	ξ₀ₐₐ,ξ₁ₐₐ=Ξ₀ₐₐ₁ₐₐ(μ₀,f₀,ł₀,ρ²₀,μ₁,f₁,ł₁,ρ²₁,μg,g,łg,ρ²g,Td);

	#sample latent normals
	y₀,y₁,yg=Y₀Y₁Yg(ξ₀ₐ,ξ₀ᵣ,ξ₁ₐ,ξ₁ᵣ,ξ₀ₐₐ,ξ₀ₐᵣ,ξ₁ₐₐ,ξ₁ₐᵣ,ξ₀ᵣᵣ,ξ₁ᵣᵣ);

	logodds=(sslogdensity(yg,1,μg,σ²,łg,ρ²g)-sslogdensity(yg,0,μg,σ²,łg,ρ²g))
	odds=exp(logodds)
	if(odds==Inf)
		gᵧ=1
	elseif(odds==-Inf)
		gᵧ=0
	else
		gᵧ=rand(Bernoulli(odds/(1+odds)))
	end
	if(gᵧ==1)
		g=FFBS(yg,collect(keys(yg)),μg,σ²,łg,ρ²g)
	else
		g=Dict{Float64,Array{Float64,2}}()
		for key in keys(yg)
			g[key]=zeros(Float64,d,1)
		end
	end
	n1=length(g)
	μg=rand(Normal((n1/σ²)*mean([yg[key]-g[key][1] for key in keys(g)])*(1/((n1/σ²)+(1/σ²ₘ))),sqrt(1/((n1/σ²)+(1/σ²ₘ)))))
	f₀=FFBS(y₀,collect(keys(y₀)),μ₀,σ²,ł₀,ρ²₀)
	f₁=FFBS(y₁,collect(keys(y₁)),μ₁,σ²,ł₁,ρ²₁)
	#=łₚ=ł+rand(Normal(0,0.03));	λₚ=sqrt(2*ν)/łₚ;	qₚ=2*ρ²*√π*λₚ^(2*p+1)*gamma(p+1)/gamma(p+0.1);=#
	#=if(log(rand(Uniform(0,1)))<glogdensity(g,gᵧ,λₚ,qₚ)-glogdensity(g,gᵧ,ł,ρ²))=#
		#=ł=łₚ=#
		#=λ=λₚ=#
		#=q=qₚ=#
		#=accepted=accepted+1=#
	#=end=#
	#=ρ²,q=rho(ρ²,g,gᵧ,ł,ρ²)=#
	#=n1=length(g)=#
	#=μ=rand(Normal((n1/σ²)*mean([y[key]-g[key][1] for key in keys(g)])*(1/((n1/σ²)+(1/σ²ₘ))),sqrt(1/((n1/σ²)+(1/σ²ₘ)))))=#
	#=ρ²Trace[iter]=ρ²=#
	#=łTrace[iter]=ł=#
	μTrace[iter]=μg
	logoddsTrace[iter]=logodds
	gᵧTrace[iter]=gᵧ
	#=plot(sort(collect(keys(g))),Φ([g[key][1] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);=#
	subplot(411)
	plot(convert(Array{Float64,1},sort(collect(keys(f₀)))),convert(Array{Float64},[Λ*Φ(μ₀+f₀[key][1]) for key in sort(collect(keys(f₀)))]),c="grey",alpha=0.1)
	#=plot(sort(collect(keys(f₀))),([f₀[key][1] for key in sort(collect(keys(f₀)))]),c="grey",alpha=0.1);=#
	subplot(412)
	plot(convert(Array{Float64,1},sort(collect(keys(f₁)))),convert(Array{Float64},[Λ*Φ(μ₁+f₁[key][1]) for key in sort(collect(keys(f₁)))]),c="grey",alpha=0.1)
	#=plot(sort(collect(keys(f₁))),([f₁[key][1] for key in sort(collect(keys(f₁)))]),c="grey",alpha=0.1);=#
	subplot(413)
	plot(convert(Array{Float64,1},sort(collect(keys(g)))),convert(Array{Float64},[Φ(μg+g[key][1]) for key in sort(collect(keys(g)))]),c="grey",alpha=0.1)
	#=plot(sort(collect(keys(g))),([g[key][1] for key in sort(collect(keys(g)))]),c="grey",alpha=0.1);=#
	subplot(414)
	Tₚ=sort(rand(Ξₚ,0,Tobs))
	gₚ=FFBS2(g,Tₚ,łg,ρ²g)
	f₀ₚ=FFBS2(f₀,Tₚ,ł₀,ρ²₀)
	f₁ₚ=FFBS2(f₁,Tₚ,ł₁,ρ²₁)
	plot(Tₚ,[Φ(μg+gₚ[key][1])*Λ*Φ(μ₀+f₀ₚ[key][1])+(1-Φ(μg+gₚ[key][1]))*Λ*Φ(μ₀+f₁ₚ[key][1]) for key in Tₚ],c="grey",alpha=0.1)
	#=plot(sort(collect(keys(ξ₀ₐₐ))),-2000*ones(length(ξ₀ₐₐ)),c="white",marker="|",linestyle="None",markersize=10)=#
	#=plot(sort(collect(keys(ξ₁ₐₐ))),-2000*ones(length(ξ₁ₐₐ)),c="white",marker="|",linestyle="None",markersize=10)=#
	#=plot(sort(collect(keys(ξ₀ₐₐ))),-2000*ones(length(ξ₀ₐₐ)),c="white",marker="|",linestyle="None",markersize=10)=#
	#=plot(sort(collect(keys(ξ₁ₐₐ))),-2000*ones(length(ξ₁ₐₐ)),c="white",marker="|",linestyle="None",markersize=10)=#
	#=plot(sort(collect(keys(ξ₀ₐₐ))),-2000*ones(length(ξ₀ₐₐ)),c="blue",marker="|",linestyle="None",markersize=10)=#
	#=plot(sort(collect(keys(ξ₁ₐₐ))),-2000*ones(length(ξ₁ₐₐ)),c="red",marker="|",linestyle="None",markersize=10)=#
end
figure()
subplot(511)
plot(ρ²Trace)
subplot(512)
plot(μTrace)
subplot(513)
plot(łTrace)
subplot(514)
plot(gᵧTrace)
ylim(-0.5,1.5)
subplot(515)
plot(logoddsTrace)
mean(ρ²Trace)
mean(μTrace)
mean(łTrace)

#=tdataStream=open("tjuliaHigh.txt","w")=#
#=gdataStream=open("gjuliaHigh.txt","w")=#
#=gammadataStream=open("gammajuliaHigh.txt","w")=#
#=zgdataStream=open("zgjuliaHigh.txt","w")=#
#=writedlm(tdataStream,t)=#
#=writedlm(gdataStream,gTrace)=#
#=writedlm(zgdataStream,zgTrace)=#
#=writedlm(gammadataStream,γTrace)=#
#=close(tdataStream)=#
#=close(gdataStream)=#
#=close(zgdataStream)=#
#=close(gammadataStream)=#
