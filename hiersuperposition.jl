importall Base.Random
importall Distributions
using DataStructures
importall PyPlot
using ProgressMeter

srand(1)
cd("/home/grad/msl33/Dropbox/pprocess/")
include("classes.jl")
include("poissonpointprocess.jl")
include("statespace.jl")
include("eigen.jl")
include("hierdatagen.jl")
#=include("fakerealdata.jl")=#



niter=3000
iter=1
gᵧ=1
σ²ₘ=300
accepted=0
r=2
gᵧTrace=zeros(niter)
ρ²Trace=zeros(niter)
łTrace=zeros(niter)
μTrace=zeros(niter)
logodds=1
logoddsTrace=zeros(niter)
σ²=1.0
#=for key in keys(f₀)=#
	#=f₀[key]=zeros(d,1)=#
#=end=#
#=for key in keys(f₁)=#
	#=f₁[key]=zeros(d,1)=#
#=end=#
#=for key in keys(g)=#
	#=g[key]=zeros(d,1)=#
#=end=#
#=for trial in Atrials=#
	#=empty!(trial.y₀)=#
	#=empty!(trial.ξ₀ₐ)=#
#=end=#
#=for trial in Btrials=#
	#=empty!(trial.y₁)=#
	#=empty!(trial.ξ₁ₐ)=#
#=end=#
#=for trial in ABtrials=#
	#=empty!(trial.y₀)=#
	#=empty!(trial.y₁)=#
	#=empty!(trial.yg)=#
	#=empty!(trial.ξ₀ₐᵣ)=#
	#=empty!(trial.ξ₀ᵣᵣ)=#
	#=empty!(trial.ξ₁ₐᵣ)=#
	#=empty!(trial.ξ₁ᵣᵣ)=#
	#=empty!(trial.ξ₀ₐₐ)=#
	#=empty!(trial.ξ₁ₐₐ)=#
	#=for key in keys(trial.g)=#
		#=g[key]=zeros(d,1)=#
	#=end=#
#=end=#

y₀=SortedDict(Dict{Float64}{Float64}())
y₁=SortedDict(Dict{Float64}{Float64}())


#initial values
μ₀=0.0
μ₁=0.0
empty!(y₀)
[merge!(y₀,trial.y₀) for trial in Atrials];
[merge!(y₀,trial.y₀) for trial in ABtrials];
f₀=SortedDict(Dict{Float64,Array{Float64,1}}())
for key in collect(keys(y₀))
	f₀[key]=zeros(Float64,3)
end
empty!(y₁)
[merge!(y₁,trial.y₁) for trial in Btrials];
[merge!(y₁,trial.y₁) for trial in ABtrials];
f₁=SortedDict(Dict{Float64,Array{Float64,1}}())
for key in collect(keys(y₁))
	f₁[key]=zeros(Float64,3)
end
ł₀=84.0
ł₁=88.0
łg=80.0
ρ²₀=0.5
ρ²₁=0.5
ρ²g=1.0
Λ=1.0
srand(2)
ł₀=0.1
ł₁=0.1
łg=30.0
ρ²₀=1.0
ρ²₁=1.0
ρ²g=3.0
Λ=400.0

trace=Dict()
trace["ABtrials"]=Dict()
for i=1:length(ABtrials)
	trace["ABtrials"][i]=zeros(Int64,niter)
end
trace["ρ²g"]=zeros(Float64,niter)
trace["ρ²₀"]=zeros(Float64,niter)
trace["ρ²₁"]=zeros(Float64,niter)
trace["łg"]=zeros(Float64,niter)
trace["ł₀"]=zeros(Float64,niter)
trace["ł₁"]=zeros(Float64,niter)
trace["μ₀"]=zeros(Float64,niter)
trace["μ₁"]=zeros(Float64,niter)
trace["Λ"]=zeros(Float64,niter)
trace["μ₀+f₀"]=zeros(Float64,niter)
trace["f₀"]=Dict()
trace["f₁"]=Dict()
trace["g"]=Dict()
trace["μg"]=Dict()
trace["λ₀mean"]=Dict{Float64,Float64}()
trace["λ₀var"]=Dict{Float64,Float64}()
trace["λ₁mean"]=Dict{Float64,Float64}()
trace["λ₁var"]=Dict{Float64,Float64}()
for t in vcat(map(x->collect(keys(x.ξ₀ₐ)),Atrials)...)
	trace["λ₀mean"][t]=0.0
	trace["λ₀var"][t]=0.0
end
for t in vcat(map(x->collect(keys(x.ξ₁ₐ)),Btrials)...)
	trace["λ₁mean"][t]=0.0
	trace["λ₁var"][t]=0.0
end
srand(2)
@showprogress 5 "Computing..." for iter=1:niter
	#=println(iter)=#
	Ξₚ=PPProcess(Λ)

	#Sample Realizations
	for trial in Atrials
		Ξ!(trial,μ₀,μ₀ₜ,f₀,ł₀,ρ²₀,Ξₚ)
	end

	for trial in Btrials
		Ξ!(trial,μ₁,μ₁ₜ,f₁,ł₁,ρ²₁,Ξₚ) 
	end

	for trial in ABtrials
		Ξ!(trial,μ₀,μ₀ₜ,f₀,ł₀,ρ²₀,μ₁,μ₁ₜ,f₁,ł₁,ρ²₁,łg,ρ²g,Ξₚ)
	end
	#=@time foo=pmap(x->Ξ(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12]),[[trial,μ₀,f₀,ł₀,ρ²₀,μ₁,f₁,ł₁,ρ²₁,łg,ρ²g,Ξₚ] for trial in ABtrials]);=#
	#=@time foo=pmap(x->Ξ(x,μ₀,f₀,ł₀,ρ²₀,μ₁,f₁,ł₁,ρ²₁,łg,ρ²g,Ξₚ),ABtrials);=#

	empty!(y₀)
	[merge!(y₀,trial.y₀) for trial in Atrials];
	if(!all(isfinite(collect(values(y₀)))))
		break
	end
	[merge!(y₀,trial.y₀) for trial in ABtrials];
	łₚ=ł₀+rand(Normal(0,0.1));
	if(log(rand(Uniform(0,1)))<sslogdensity(ms(y₀,μ₀ₜ),1.0,μ₀,σ²,łₚ,ρ²₀)-sslogdensity(ms(y₀,μ₀ₜ),1.0,μ₀,σ²,ł₀,ρ²₀)+logpdf(Gamma(2,0.5),łₚ)-logpdf(Gamma(2,0.5),ł₀))
		ł₀=łₚ
	end
	ρ²ₚ=ρ²₀+rand(Normal(0,0.1));
	if(log(rand(Uniform(0,1)))<sslogdensity(ms(y₀,μ₀ₜ),1.0,μ₀,σ²,ł₀,ρ²ₚ)-sslogdensity(ms(y₀,μ₀ₜ),1.0,μ₀,σ²,ł₀,ρ²₀)+logpdf(Gamma(2,1),ρ²ₚ)-logpdf(Gamma(2,1),ρ²₀))
		ρ²₀=ρ²ₚ
	end
	μ₀=rand(mu(ms(y₀,μ₀ₜ),1,σ²,ł₀,ρ²₀,σ²ₘ))
	f₀=FFBS(ms(y₀,μ₀ₜ),μ₀,σ²,ł₀,ρ²₀);

	empty!(y₁)
	[merge!(y₁,trial.y₁) for trial in Btrials];
	if(!all(isfinite(collect(values(y₁)))))
		break
	end
	[merge!(y₁,trial.y₁) for trial in ABtrials];
	łₚ=ł₁+rand(Normal(0,0.1));
	if(log(rand(Uniform(0,1)))<sslogdensity(ms(y₁,μ₁ₜ),1.0,μ₁,σ²,łₚ,ρ²₁)-sslogdensity(ms(y₁,μ₁ₜ),1.0,μ₁,σ²,ł₁,ρ²₁)+logpdf(Gamma(2,0.5),łₚ)-logpdf(Gamma(2,0.5),ł₁))
		ł₁=łₚ
	end
	ρ²ₚ=ρ²₁+rand(Normal(0,0.1));
	if(log(rand(Uniform(0,1)))<sslogdensity(ms(y₁,μ₁ₜ),1.0,μ₁,σ²,ł₁,ρ²ₚ)-sslogdensity(ms(y₁,μ₁ₜ),1.0,μ₁,σ²,ł₁,ρ²₁)+logpdf(Gamma(2,1),ρ²ₚ)-logpdf(Gamma(2,1),ρ²₁))
		ρ²₁=ρ²ₚ
	end
	μ₁=rand(mu(ms(y₁,μ₁ₜ),1,σ²,ł₁,ρ²₁,σ²ₘ))
	f₁=FFBS(ms(y₁,μ₁ₜ),μ₁,σ²,ł₁,ρ²₁);

	for trial in ABtrials
		𝐺!(trial,σ²,σ²ₘ,łg,ρ²g)
	end
	łgₚ=łg+rand(Normal(0,0.5));
	if(log(rand(Uniform(0,1)))<sum([sslogdensity(trial,σ²,łgₚ,ρ²g) for trial in ABtrials])-sum([sslogdensity(trial,σ²,łg,ρ²g) for trial in ABtrials])+logpdf(Gamma(2,1),łgₚ)-logpdf(Gamma(2,1),łg))
		łg=łgₚ
	end
	ρ²gₚ=ρ²g+rand(Normal(0,0.5));
	if(log(rand(Uniform(0,1)))<sum([sslogdensity(trial,σ²,łg,ρ²gₚ) for trial in ABtrials])-sum([sslogdensity(trial,σ²,łg,ρ²g) for trial in ABtrials])+logpdf(Gamma(2,1),ρ²gₚ)-logpdf(Gamma(2,1),ρ²g))
		ρ²g=ρ²gₚ
	end

	Λshape=sum(map(x->length(x.ξ₀ₐ)+length(x.ξ₀ᵣ),Atrials))+sum(map(x->length(x.ξ₁ₐ)+length(x.ξ₁ᵣ),Btrials))+sum(map(x->length(x.ξ₀ₐᵣ)+length(x.ξ₀ᵣᵣ)+length(x.ξ₁ₐᵣ)+length(x.ξ₁ᵣᵣ)+length(x.ξ₀ₐₐ)+length(x.ξ₁ₐₐ),ABtrials))
	Λrate=sum(map(x->x.Tobs,Atrials))+sum(map(x->x.Tobs,Btrials))+sum(map(x->2*x.Tobs,ABtrials))
	Λ=rand(Gamma(Λshape+0.001,1/(Λrate+0.001)))
	for i=1:length(ABtrials)
		trace["ABtrials"][i][iter]=ABtrials[i].gᵧ
	end
	trace["ρ²g"][iter]=ρ²g
	trace["ρ²₀"][iter]=ρ²₀
	trace["ρ²₁"][iter]=ρ²₁
	trace["łg"][iter]=łg
	trace["ł₀"][iter]=ł₀
	trace["ł₁"][iter]=ł₁
	trace["μ₀"][iter]=μ₀
	trace["μ₁"][iter]=μ₁
	trace["Λ"][iter]=Λ
	trace["f₀"][iter]=f₀
	trace["f₁"][iter]=f₁
	trace["g"][iter]=ABtrials[1].g
	trace["μg"][iter]=ABtrials[1].μg
	#=if(iter>10)=#
	#=for t in unique(vcat(map(x->collect(keys(x.ξ₀ₐ)),Atrials)...))=#
		#=trace["λ₀mean"][t]=trace["λ₀mean"][t]+Λ*Φ(μ₀+f₀[t][1])=#
	#=end=#
	#=for t in unique(vcat(map(x->collect(keys(x.ξ₁ₐ)),Btrials)...))=#
		#=trace["λ₁mean"][t]=trace["λ₁mean"][t]+Λ*Φ(μ₁+f₁[t][1])=#
	#=end=#
	#=end=#
end
	figure(3)
	subplot(611)
	plot(trace["μ₀"][1:iter],c="grey")
	subplot(612)
	plot(trace["ł₀"][1:iter],c="grey")
	subplot(613)
	plot(trace["ρ²₀"][1:iter],c="grey")
	subplot(614)
	plot(trace["μ₁"][1:iter],c="grey")
	subplot(615)
	plot(trace["ł₁"][1:iter],c="grey")
	subplot(616)
	plot(trace["ρ²₁"][1:iter],c="grey")

	figure(1)
	subplot(211)
	plot(trace["łg"][1:iter],c="grey")
	subplot(212)
	plot(trace["ρ²g"][1:iter],c="grey")

	figure(2)
	plot(trace["Λ"][1:iter],c="grey")

	figure(4)
	subplot(611)
	plot(trace["ABtrials"][1][1:iter],c="grey")
	subplot(612)
	plot(trace["ABtrials"][2][1:iter],c="grey")
	subplot(613)
	plot(trace["ABtrials"][3][1:iter],c="grey")
	subplot(614)
	plot(trace["ABtrials"][4][1:iter],c="grey")
	subplot(615)
	plot(trace["ABtrials"][5][1:iter],c="grey")
	subplot(616)
	plot(trace["ABtrials"][6][1:iter],c="grey")

	figure()
	foo=zeros(0)
	for i=1:length(ABtrials)
		println(mean(trace["ABtrials"][i]))
		push!(foo,mean(trace["ABtrials"][i]))
	end
	plt[:hist](foo)

	figure(1)
	subplot(411)
	for i=1000:niter
		plot(collect(keys(trace["f₀"][i])),[trace["Λ"][i]*Φ(trace["μ₀"][i]+trace["f₀"][i][key][1]) for key in collect(keys(trace["f₀"][i]))],c="grey",alpha=0.01)
	end
	subplot(413)
	for i=1000:niter
		plot(collect(keys(trace["f₁"][i])),[trace["Λ"][i]*Φ(trace["μ₁"][i]+trace["f₁"][i][key][1]) for key in collect(keys(trace["f₁"][i]))],c="grey",alpha=0.01)
	end
	figure(4345)
	subplot(211)
	for i=1000:niter
		plot(collect(keys(trace["g"][i])),[(trace["g"][i][key][1]) for key in collect(keys(trace["g"][i]))],c="blue",alpha=0.01)
		#=plot(collect(keys(trace["g"][i])),[Φ(trace["μg"][i]+trace["g"][i][key][1]) for key in collect(keys(trace["g"][i]))],c="blue",alpha=0.01)=#
	end
	subplot(212)
	plot(1:niter,[trace["μg"][key] for key in keys(trace["μg"])])
	λ₀mean=copy(trace["λ₀mean"])
	λ₁mean=copy(trace["λ₁mean"])
	λ₀var=copy(trace["λ₀var"])
	λ₁var=copy(trace["λ₁var"])
	λ₀u=copy(λ₀mean)
	λ₀b=copy(λ₀mean)
	λ₁u=copy(λ₁mean)
	λ₁b=copy(λ₁mean)
	for key in keys(λ₀mean)
		λ₀mean[key]=λ₀mean[key]/iter
		λ₀var[key]=λ₀var[key]/iter
		λ₀u[key]=λ₀mean[key]+2*sqrt(λ₀var[key])
		λ₀b[key]=λ₀mean[key]-2*sqrt(λ₀var[key])
	end
	for key in keys(λ₁mean)
		λ₁mean[key]=λ₁mean[key]/iter
		λ₁var[key]=λ₁var[key]/iter
		λ₁u[key]=λ₁mean[key]+2*sqrt(λ₁var[key])
		λ₁b[key]=λ₁mean[key]-2*sqrt(λ₁var[key])
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
