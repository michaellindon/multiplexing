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
#=include("hierdatagen.jl")=#
#=include("fakerealdata.jl")=#



cell=cells["YKIC140130Loc_DoubleSound"]

fooA=SortedDict(Dict(zip(collect(0:0.005:1),cell["A"]["mean"])))
fooB=SortedDict(Dict(zip(collect(0:0.005:1),(cell["B"]["mean"]))))
μ₀ₜ(t)=invlogcdf(Normal(0,1),log((1/Λ)*deref((fooA,searchsortedlast(fooA,t)))[2]))
μ₁ₜ(t)=invlogcdf(Normal(0,1),log((1/Λ)*(deref((fooB,searchsortedlast(fooB,t)))[2])))


niter=1000
iter=1
σ²ₘ=9
k=0.1
a=1
σ²=1.0

y₀=SortedDict(Dict{Float64}{Float64}())
y₁=SortedDict(Dict{Float64}{Float64}())


#initial values
μ₀=0.0
μ₁=0.0
empty!(y₀)
[merge!(y₀,trial.y₀) for trial in ABtrials];
f₀=SortedDict(Dict{Float64,Array{Float64,1}}())
for key in collect(keys(y₀))
	f₀[key]=zeros(Float64,3)
end
empty!(y₁)
[merge!(y₁,trial.y₁) for trial in ABtrials];
f₁=SortedDict(Dict{Float64,Array{Float64,1}}())
for key in collect(keys(y₁))
	f₁[key]=zeros(Float64,3)
end
ł₀=0.1
ł₁=0.1
p=0.5
łg=30.0
ρ²₀=cell["A"]["variance"]
ρ²₁=cell["B"]["variance"]
ρ²g=3.0

trace=Dict()
trace["ABtrials"]=Dict()
for i=1:length(ABtrials)
	trace["ABtrials"][i]=zeros(Int64,niter)
end
trace["ρ²g"]=zeros(Float64,niter)
trace["p"]=zeros(Float64,niter)
trace["ρ²₀"]=zeros(Float64,niter)
trace["ρ²₁"]=zeros(Float64,niter)
trace["łg"]=zeros(Float64,niter)
trace["ł₀"]=zeros(Float64,niter)
trace["ł₁"]=zeros(Float64,niter)
trace["μ₀"]=zeros(Float64,niter)
trace["μ₁"]=zeros(Float64,niter)
trace["Λ"]=zeros(Float64,niter)
trace["μ₀+f₀"]=zeros(Float64,niter)
trace["λ₀"]=Dict()
trace["λ₁"]=Dict()
trace["g"]=Dict()
trace["μg"]=Dict()
trace["μprior"]=Dict()
srand(2)
ftmp₀=lazyGP(f₀,ł₀)
m₀ = x -> μ₀ₜ(x)+ftmp₀(x)
λ₀= x-> Λ*Φ(m₀(x))
ftmp₁=lazyGP(f₁,ł₁)
m₁ = x-> μ₁ₜ(x)+ftmp₁(x)
λ₁ = x-> Λ*Φ(m₁(x))
trueABtrials=deepcopy(ABtrials)
@showprogress 5 "Computing..." for iter=1:niter

	#Sample Realizations
	for trial in ABtrials
		Ξ!(trial,m₀,λ₀,m₁,λ₁,łg,ρ²g,Λ)
	end

	empty!(y₀)
	[merge!(y₀,trial.y₀) for trial in ABtrials];
	łₚ=ł₀+rand(Normal(0,0.1));
	if(log(rand(Uniform(0,1)))<sslogdensity(ms(y₀,μ₀ₜ),1.0,μ₀,σ²,łₚ,ρ²₀)-sslogdensity(ms(y₀,μ₀ₜ),1.0,μ₀,σ²,ł₀,ρ²₀)+logpdf(Gamma(2,0.5),łₚ)-logpdf(Gamma(2,0.5),ł₀))
		ł₀=łₚ
	end
	ρ²ₚ=ρ²₀+rand(Normal(0,0.1));
	if(log(rand(Uniform(0,1)))<sslogdensity(ms(y₀,μ₀ₜ),1.0,μ₀,σ²,ł₀,ρ²ₚ)-sslogdensity(ms(y₀,μ₀ₜ),1.0,μ₀,σ²,ł₀,ρ²₀)+logpdf(Gamma(1,1/100),ρ²ₚ)-logpdf(Gamma(1,1/100),ρ²₀))
		ρ²₀=ρ²ₚ
	end
	f₀=lazyGP(FFBS(ms(y₀,μ₀ₜ),μ₀,σ²,ł₀,ρ²₀),ł₀)
	m₀= x->μ₀ₜ(x)+f₀(x)
	λ₀= x-> Λ*Φ(m₀(x))

	empty!(y₁)
	[merge!(y₁,trial.y₁) for trial in ABtrials];
	łₚ=ł₁+rand(Normal(0,0.1));
	if(log(rand(Uniform(0,1)))<sslogdensity(ms(y₁,μ₁ₜ),1.0,μ₁,σ²,łₚ,ρ²₁)-sslogdensity(ms(y₁,μ₁ₜ),1.0,μ₁,σ²,ł₁,ρ²₁)+logpdf(Gamma(2,0.5),łₚ)-logpdf(Gamma(2,0.5),ł₁))
		ł₁=łₚ
	end
	ρ²ₚ=ρ²₁+rand(Normal(0,0.1));
	if(log(rand(Uniform(0,1)))<sslogdensity(ms(y₁,μ₁ₜ),1.0,μ₁,σ²,ł₁,ρ²ₚ)-sslogdensity(ms(y₁,μ₁ₜ),1.0,μ₁,σ²,ł₁,ρ²₁)+logpdf(Gamma(1,1/100),ρ²ₚ)-logpdf(Gamma(1,1/100),ρ²₁))
		ρ²₁=ρ²ₚ
	end
	f₁=lazyGP(FFBS(ms(y₁,μ₁ₜ),μ₁,σ²,ł₁,ρ²₁),ł₁)
	m₁= x->μ₁ₜ(x)+f₁(x)
	λ₁ = x-> Λ*Φ(m₁(x))

	for trial in ABtrials
		mulist=map(x->x.μprior,ABtrials)
		weights=map(x->logpdf(Normal(x,sqrt(k*σ²ₘ)),trial.μg),mulist)
		indx=find(map(x->x.id==trial.id,ABtrials))[1]
		weights[indx]=logpdf(Normal(0,sqrt((1+k)*σ²ₘ)),trial.μg)
		weights=exp(weights-maximum(weights))
		weights[indx]=weights[indx]*a
		weights=weights/sum(weights)
		sindx=rand(Categorical(weights))
		if(indx==sindx)
			trial.μprior=rand(Normal((1/(1/(k*σ²ₘ) + 1/σ²ₘ))*(trial.μg/(k*σ²ₘ)),sqrt(1/(1/(k*σ²ₘ) + (1/σ²ₘ))) ))
		else
			trial.μprior=mulist[sindx]
		end
	end

	for trial in ABtrials
		𝐺!(trial,σ²,k*σ²ₘ,łg,ρ²g,p)
	end
	p=rand(Beta(1+sum(map(x->x.gᵧ,ABtrials)),1+length(ABtrials)-sum(map(x->x.gᵧ,ABtrials))))
	łgₚ=łg+rand(Normal(0,0.5));
	if(log(rand(Uniform(0,1)))<sum([sslogdensity(trial,σ²,łgₚ,ρ²g) for trial in ABtrials])-sum([sslogdensity(trial,σ²,łg,ρ²g) for trial in ABtrials])+logpdf(Gamma(1,1/5),łgₚ)-logpdf(Gamma(1,1/5),łg))
		łg=łgₚ
	end
	ρ²gₚ=ρ²g+rand(Normal(0,0.5));
	if(log(rand(Uniform(0,1)))<sum([sslogdensity(trial,σ²,łg,ρ²gₚ) for trial in ABtrials])-sum([sslogdensity(trial,σ²,łg,ρ²g) for trial in ABtrials])+logpdf(Gamma(3,1),ρ²gₚ)-logpdf(Gamma(3,1),ρ²g))
		ρ²g=ρ²gₚ
	end

	for i=1:length(ABtrials)
		trace["ABtrials"][i][iter]=ABtrials[i].gᵧ
	end
	trace["ρ²g"][iter]=ρ²g
	trace["ρ²₀"][iter]=ρ²₀
	trace["ρ²₁"][iter]=ρ²₁
	trace["p"][iter]=p
	trace["łg"][iter]=łg
	trace["ł₀"][iter]=ł₀
	trace["ł₁"][iter]=ł₁
	trace["μ₀"][iter]=μ₀
	trace["μ₁"][iter]=μ₁
	trace["Λ"][iter]=Λ
	trace["λ₀"][iter]=map(x->λ₀(x),collect(0:0.01:1))
	trace["λ₁"][iter]=map(x->λ₁(x),collect(0:0.01:1))
	trace["g"][iter]=map(x->(x.gᵧ==1 ? FFBS2(x.g,collect(0:0.01:1),łg,ρ²g) : SortedDict(Dict(map(x->(x,zeros(Float64,3)),collect(0:0.01:1))))),ABtrials)
	#=trace["g"][iter]=map(x->FFBS2(x.g,collect(0:0.01:1),łg,ρ²g),ABtrials)=#
	trace["μg"][iter]=map(x->x.μg,ABtrials)
	trace["μprior"][iter]=map(x->x.μprior,ABtrials)
end

n=length(ABtrials)
#=plt[:hist](map(x-> (@as _ rand(Categorical([ones(n)/(a+n);a/(a+n)])) (_ == (n+1) ? rand(Normal(0,sqrt((1+k)*σ²ₘ))) : rand(Normal(trace["μprior"][x][_],sqrt(k*σ²ₘ))))) , collect(1:niter)),100)=#
for iter=1:3030
	plot(collect(-3:0.01:3),map(y->(a/(a+n))*pdf(Normal(0,sqrt(1+k*σ²ₘ)),y)+sum(map(x-> (1/(a+n))*pdf(Normal(x,sqrt(k*σ²ₘ)),y),trace["μprior"][iter])),collect(-3:0.01:3)),alpha=0.1,c="blue")
end

	vec(map(x->trace["μg"][x],1:iter))
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

	figure(11)
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

	figure()
	foo=zeros(0)
	for i=100:iter-1
		append!(foo,trace["μg"][i])
	end
	subplot(211)
	plt[:hist](foo,100)
	subplot(212)
	plt[:hist](Φ(foo),100)

	figure(1010)
	subplot(211)
	for i=collect(1:1:niter)
		#=plot(collect(keys(trace["f₀"][i])),[trace["f₀"][i][key][1] for key in collect(keys(trace["f₀"][i]))],c="grey",alpha=0.01)=#
		plot(trace["λ₀"][i],c="blue",alpha=0.01)
	end
	subplot(212)
	for i=collect(1:1:niter)
		#=plot(collect(keys(trace["f₁"][i])),[trace["f₁"][i][key][1] for key in collect(keys(trace["f₁"][i]))],c="grey",alpha=0.01)=#
		plot(trace["λ₁"][i],c="blue",alpha=0.01)
	end

	figure(1)
	subplot(211)
	for i=100:iter
		plot(collect(keys(trace["f₀"][i])),[trace["Λ"][i]*Φ(trace["μ₀"][i]+μ₀ₜ(key)+trace["f₀"][i][key][1]) for key in collect(keys(trace["f₀"][i]))],c="grey",alpha=0.01)
	end
	subplot(212)
	for i=100:iter
		plot(collect(keys(trace["f₁"][i])),[trace["Λ"][i]*Φ(trace["μ₁"][i]+μ₁ₜ(key)+trace["f₁"][i][key][1]) for key in collect(keys(trace["f₁"][i]))],c="grey",alpha=0.01)
	end
	for intAB=1:length(ABtrials)
	figure(intAB)
	subplot(311)
	for i=collect(1:1:niter)
		plot(collect(keys(trace["g"][i][intAB])),[(trace["g"][i][intAB][key][1]) for key in collect(keys(trace["g"][i][intAB]))],c="blue",alpha=0.01)
	end
	plot(collect(keys(trueABtrials[intAB].g)),[trueABtrials[intAB].g[key][1]*trueABtrials[intAB].gᵧ for key in collect(keys(trueABtrials[intAB].g))],c="red")
	subplot(312)
	for i=collect(1:1:niter)
		plot(collect(1/1000:1/1000:1),trace["μg"][i][intAB]*ones(1000),c="blue",alpha=0.01)
	end
	plot(collect(1/1000:1/1000:1),trueABtrials[intAB].μg*ones(1000),c="red")
	subplot(313)
	for i=collect(1:1:niter)
		plot(collect(keys(trace["g"][i][intAB])),[Φ(trace["μg"][i][intAB]+trace["g"][i][intAB][key][1]) for key in collect(keys(trace["g"][i][intAB]))],c="blue",alpha=0.01)
	end
	plot(collect(keys(trueABtrials[intAB].g)),[Φ(trueABtrials[intAB].μg+trueABtrials[intAB].g[key][1]) for key in collect(keys(trueABtrials[intAB].g))],c="red")
	ylim(-0.1,1.1)
	savefig(string("plots/16/",intAB,".pdf"))
	close(intAB)
	end


	mean(trace["ABtrials"][intAB][1:iter])

	for i=100:iter
		plot(collect(keys(trace["g"][1][i])),[trace["Λ"][i]*Φ(trace["μ₀"][i]+μ₀ₜ(key)+trace["f₀"][i][key][1]) for key in collect(keys(trace["f₀"][i]))],c="grey",alpha=0.01)
	end

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
