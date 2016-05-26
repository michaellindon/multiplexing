importall Base.Random
importall Distributions
using DataStructures
importall PyPlot
using ProgressMeter

srand(1)
#=cd("/home/grad/msl33/Dropbox/pprocess/")=#
include("classes.jl")
include("poissonpointprocess.jl")
include("statespace.jl")
#=include("eigen.jl")=#
include("singledatagen.jl")
#=include("hierdatagen.jl")=#
#=include("fakerealdata.jl")=#



niter=500
iter=1
σ²ₘ=300
accepted=0
r=2
σ²=1.0

y₀=SortedDict(Dict{Float64}{Float64}())

#initial values
μ₀=0.0
empty!(y₀)
[merge!(y₀,trial.y₀) for trial in Atrials];
f₀=SortedDict(Dict{Float64,Array{Float64,2}}())
for key in union(0,collect(keys(y₀)))
	f₀[key]=zeros(Float64,3,1)
end
ρ²₀=9.0
ł₀=0.1
Λ=1000
srand(2)
trace=Dict()
trace["ρ²₀"]=zeros(Float64,niter)
trace["ł₀"]=zeros(Float64,niter)
trace["μ₀"]=zeros(Float64,niter)
trace["Λ"]=zeros(Float64,niter)
trace["f₀"]=Dict()
trace["λ₀mean"]=Dict{Float64,Float64}()
trace["λ₀var"]=Dict{Float64,Float64}()
for t in vcat(map(x->collect(keys(x.ξ₀ₐ)),Atrials)...)
	trace["λ₀mean"][t]=0.0
	trace["λ₀var"][t]=0.0
end
srand(2)
@showprogress 5 "Computing..." for iter=1:niter
	#=println(iter)=#
	Ξₚ=PPProcess(Λ)
	#Sample Realizations
	for trial in Atrials
		Ξ!(trial,μ₀,f₀,zhutransition,zhuinnovation(ł₀,ρ²₀),Ξₚ)
	end
	empty!(y₀)
	[merge!(y₀,trial.y₀) for trial in Atrials];
	if(!all(isfinite(collect(values(y₀)))))
		break
	end
	łₚ=ł₀+rand(Normal(0,1000));
	if(łₚ>0)
		if(log(rand(Uniform(0,1)))<sslogdensity(y₀,1.0,μ₀,σ²,1000*eye(3),zhutransition,zhuinnovation(łₚ,ρ²₀))-sslogdensity(y₀,1.0,μ₀,σ²,1000*eye(3),zhutransition,zhuinnovation(ł₀,ρ²₀))+logpdf(Gamma(0.001,1/0.0000001),łₚ)-logpdf(Gamma(0.001,1/0.0000001),ł₀))
			ł₀=łₚ
		end
	end
	ρ²ₚ=ρ²₀+rand(Normal(0,1));
	if(ρ²ₚ>0)
		if(log(rand(Uniform(0,1)))<sslogdensity(y₀,1.0,μ₀,σ²,1000*eye(3),zhutransition,zhuinnovation(ł₀,ρ²ₚ))-sslogdensity(y₀,1.0,μ₀,σ²,1000*eye(3),zhutransition,zhuinnovation(ł₀,ρ²₀))+logpdf(Gamma(1,1/0.001),ρ²ₚ)-logpdf(Gamma(1,1/0.001),ρ²₀))
			ρ²₀=ρ²ₚ
		end
	end
	#=μ₀=rand(mu(y₀,1,σ²,ł₀,ρ²₀,σ²ₘ))=#

	f₀=FFBS(y₀,μ₀,σ²,0,1000*eye(3),zhutransition,zhuinnovation(ł₀,ρ²₀))
	Λshape=sum(map(x->length(x.ξ₀ₐ)+length(x.ξ₀ᵣ),Atrials))
	Λrate=sum(map(x->x.Tobs,Atrials))
	Λ=rand(Gamma(Λshape+0.00001,1/(Λrate+0.00001)))
	trace["ρ²₀"][iter]=ρ²₀
	trace["ł₀"][iter]=ł₀
	trace["μ₀"][iter]=μ₀
	trace["Λ"][iter]=Λ
	trace["f₀"][iter]=f₀
	if(iter>10)
		for t in unique(vcat(map(x->collect(keys(x.ξ₀ₐ)),Atrials)...))
			trace["λ₀mean"][t]=trace["λ₀mean"][t]+Λ*Φ(μ₀+f₀[t][1])
		end
	end
end
figure(3)
subplot(311)
plot(trace["μ₀"][1:iter],c="grey")
subplot(312)
plot(trace["ł₀"][1:iter],c="grey")
subplot(313)
plot(trace["ρ²₀"][1:iter],c="grey")

figure(2)
plot(trace["Λ"][1:iter],c="grey")


figure(43)
subplot(211)
for i=1:iter
	plot(collect(keys(trace["f₀"][i])),[trace["Λ"][i]*Φ(trace["μ₀"][i]+trace["f₀"][i][key][1]) for key in collect(keys(trace["f₀"][i]))],c="blue",alpha=0.01)
end

λ₀mean=copy(trace["λ₀mean"])
λ₀var=copy(trace["λ₀var"])
λ₀u=copy(λ₀mean)
λ₀b=copy(λ₀mean)
for key in keys(λ₀mean)
	λ₀mean[key]=λ₀mean[key]/iter
	λ₀var[key]=λ₀var[key]/iter
	λ₀u[key]=λ₀mean[key]+2*sqrt(λ₀var[key])
	λ₀b[key]=λ₀mean[key]-2*sqrt(λ₀var[key])
end

figure(43)
subplot(411)
plot(sort(collect(keys(λ₀mean))),[λ₀mean[key] for key in sort(collect(keys(λ₀mean)))],c="blue")
#=plot(sort(collect(keys(λ₀u))),[λ₀u[key] for key in sort(collect(keys(λ₀mean)))],linestyle="--",c="blue")=#
	#=plot(sort(collect(keys(λ₀b))),[λ₀b[key] for key in sort(collect(keys(λ₀mean)))],c="blue",linestyle="--")=#

		subplot(412)
		plt[:hist](vcat(map(x->collect(keys(x.ξ₀ₐ)),Atrials)...),100)
		subplot(413)
		plot(sort(collect(keys(λ₁mean))),[λ₁mean[key] for key in sort(collect(keys(λ₁mean)))],c="red")
		#=plot(sort(collect(keys(λ₁u))),[λ₁u[key] for key in sort(collect(keys(λ₁mean)))],linestyle="--",c="red")=#
			#=plot(sort(collect(keys(λ₁b))),[λ₁b[key] for key in sort(collect(keys(λ₁mean)))],c="red",linestyle="--")=#
				subplot(414)
				plt[:hist](vcat(map(x->collect(keys(x.ξ₁ₐ)),Btrials)...),100)



				figure(44)
				subplot(211)
				plt[:hist](vcat(map(x->collect(keys(x.ξ₀ₐ)),Atrials)...),100)
				ylim(0,50)
				subplot(212)
				plt[:hist](vcat(map(x->collect(keys(x.ξ₁ₐ)),Btrials)...),100)
				ylim(0,50)
