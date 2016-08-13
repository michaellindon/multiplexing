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
μ₀(t)=invlogcdf(Normal(0,1),log((1/Λ)*deref((fooA,searchsortedlast(fooA,t)))[2]))
μ₁(t)=invlogcdf(Normal(0,1),log((1/Λ)*(deref((fooB,searchsortedlast(fooB,t)))[2])))
ρ²₀=cell["A"]["variance"]
ρ²₁=cell["B"]["variance"]



Atrials=map(trial->filter(z->z>0,trial)/1000,map(y->azeem[1][3]["TIMES"][find(x-> x==y,azeem[1][3]["TRIAL2"])],1:50))
Btrials=map(trial->filter(z->z>0,trial)/1000,map(y->azeem[1][3]["TIMES"][find(x-> x==y,azeem[1][3]["TRIAL2"])],51:100))
ABtrials=map(trial->filter(z->z>0,trial)/1000,map(y->azeem[1][3]["TIMES"][find(x-> x==y,azeem[1][3]["TRIAL2"])],101:150))
Atrials=[Atrial(1,trial,1) for trial in Atrials]
Btrials=[Btrial(1,trial,1) for trial in Btrials]
ABtrials=[ABtrial(1,trial,1) for trial in ABtrials]
μ₀,μ₁,Λ,ρ²₀,ρ²₁=preprocess(Atrials,Btrials)

niter=10000
iter=1
k=0.1
a=1
σ²=1.0

y₀=SortedDict(Dict{Float64}{Float64}())
y₁=SortedDict(Dict{Float64}{Float64}())


μ₀ = x->0
μ₁ = x->0
ρ²₀=1
ρ²₁=1

trace=Dict()
trace["ρ²g"]=zeros(Float64,niter)
trace["p"]=zeros(Float64,niter)
trace["ρ²₀"]=zeros(Float64,niter)
trace["ρ²₁"]=zeros(Float64,niter)
trace["łg"]=zeros(Float64,niter)
trace["ł₀"]=zeros(Float64,niter)
trace["ł₁"]=zeros(Float64,niter)
trace["Λ"]=zeros(Float64,niter)
trace["μ₀+f₀"]=zeros(Float64,niter)
trace["λ₀"]=Dict()
trace["μ₀bar"]=Dict()
trace["μ₁bar"]=Dict()
trace["λ₁"]=Dict()
trace["g"]=Dict()
trace["μg"]=Dict()
trace["ABtrials"]=Dict()
trace["μprior"]=Dict()
trace[:α]=Dict()
trace[:g]=Dict()
trace[:μg]=Dict()
trace[:gᵧ]=Dict()
λ₀ = x->Λ
λ₁ = x->Λ
m₀ = x->0
m₁ = x->0
ł₀=1
ł₁=1
łg=1
ρ²g=1
σ²ₘ=100
k=0.1
a=1
p=0.5
prior=Dict(:p => Beta(1,1), :ł₀ => Gamma(2,1), :ρ²₀ => Gamma(1,1) ,:ł₁ => Gamma(2,1),:ρ²₁ => Gamma(1,1) ,:łg => Gamma(5,5/100),:ρ²g => Gamma(2,1) ,:μg => Normal(0,sqrt(3)))
#=trueABtrials=deepcopy(ABtrials)=#
@showprogress 5 "Computing..." for iter=1:niter

	#Sample Realizations
	for trial in Atrials
		Ξ!(trial,m₀,λ₀,Λ)
	end
	for trial in Btrials
		Ξ!(trial,m₁,λ₁,Λ)
	end
	for trial in ABtrials
		Ξ!(trial,m₀,λ₀,m₁,λ₁,Λ)
	end

	y₀=reduce(merge,map(x->x.y₀,[Atrials;ABtrials]))
	łₚ=ł₀+rand(Normal(0,0.1));
	ł₀= (log(rand(Uniform(0,1))) < -(map(z->sslogdensity(ms(y₀,μ₀),1,0,σ²,z,ρ²₀)+logpdf(prior[:ł₀],z),[łₚ,ł₀])...)) ?  łₚ : ł₀
	ρ²ₚ=ρ²₀+rand(Normal(0,0.1));
	ρ²₀= (log(rand(Uniform(0,1))) < -(map(z->sslogdensity(ms(y₀,μ₀),1,0,σ²,ł₀,z)+logpdf(prior[:ρ²₀],z),[ρ²ₚ,ρ²₀])...)) ? ρ²ₚ : ρ²₀
	μ₀bar = rand(mu(y₀,1,σ²,ł₀,ρ²₀,σ²ₘ,0))
	μ₀ = x-> μ₀bar 
	f₀=lazyGP(FFBS(ms(y₀,μ₀),0,σ²,ł₀,ρ²₀),ł₀)
	m₀= x->μ₀(x)+f₀(x)
	λ₀= x-> Λ*Φ(m₀(x))

	y₁=reduce(merge,map(x->x.y₁,[Btrials;ABtrials]))
	łₚ=ł₁+rand(Normal(0,0.1));
	ł₁= (log(rand(Uniform(0,1))) < -(map(z->sslogdensity(ms(y₁,μ₁),1,0,σ²,z,ρ²₁)+logpdf(prior[:ł₁],z),[łₚ,ł₁])...)) ? łₚ : ł₁
	ρ²ₚ=ρ²₁+rand(Normal(0,0.1));
	ρ²₁= (log(rand(Uniform(0,1))) < -(map(z->sslogdensity(ms(y₁,μ₁),1,0,σ²,ł₁,z)+logpdf(prior[:ρ²₁],z),[ρ²ₚ,ρ²₁])...)) ? ρ²ₚ : ρ²₁
	μ₁bar = rand(mu(y₁,1,σ²,ł₁,ρ²₁,σ²ₘ,0))
	μ₁ = x -> μ₁bar
	f₁=lazyGP(FFBS(ms(y₁,μ₁),0,σ²,ł₁,ρ²₁),ł₁)
	m₁= x->μ₁(x)+f₁(x)
	λ₁ = x-> Λ*Φ(m₁(x))

	#=for trial in ABtrials=#
		#=mulist=map(x->x.μprior,ABtrials)=#
		#=weights=map(x->logpdf(Normal(x,sqrt(k*σ²ₘ)),trial.μg),mulist)=#
		#=indx=find(map(x->x.id==trial.id,ABtrials))[1]=#
		#=weights[indx]=logpdf(Normal(0,sqrt((1+k)*σ²ₘ)),trial.μg)=#
		#=weights=exp(weights-maximum(weights))=#
		#=weights[indx]=weights[indx]*a=#
		#=weights=weights/sum(weights)=#
		#=sindx=rand(Categorical(weights))=#
		#=if(indx==sindx)=#
			#=trial.μprior=rand(Normal((1/(1/(k*σ²ₘ) + 1/σ²ₘ))*(trial.μg/(k*σ²ₘ)),sqrt(1/(1/(k*σ²ₘ) + (1/σ²ₘ))) ))=#
		#=else=#
			#=trial.μprior=mulist[sindx]=#
		#=end=#
	#=end=#

	for trial in ABtrials
		𝐺!(trial,σ²,k*σ²ₘ,łg,ρ²g,p)
	end
	p=rand(Beta(1+sum(map(x->x.gᵧ,ABtrials)),1+length(ABtrials)-sum(map(x->x.gᵧ,ABtrials))))
	if(any(map(x->x.gᵧ==1,ABtrials)))
		łgₚ=łg+rand(Normal(0,0.5));
		łg=(log(rand(Uniform(0,1))) < -(map(z->logpdf(prior[:łg],z)+sum(map(x->sslogdensity(x.yg,1,x.μg,σ²,z,ρ²g),filter(x-> x.gᵧ==1,ABtrials))),[łgₚ,łg])...)) ? łgₚ : łg
		ρ²gₚ=ρ²g+rand(Normal(0,0.5));
		ρ²g=(log(rand(Uniform(0,1))) < -(map(z->logpdf(prior[:ρ²g],z)+sum(map(x->sslogdensity(x.yg,1,x.μg,σ²,łg,z),filter(x-> x.gᵧ==1,ABtrials))),[ρ²gₚ,ρ²g])...))?ρ²gₚ : ρ²g
	else
		łgₚ=łg+rand(Normal(0,0.5));
		łg=(log(rand(Uniform(0,1))) < -(map(z->logpdf(prior[:łg],z),[łgₚ,łg])...)) ? łgₚ : łg
		ρ²gₚ=ρ²g+rand(Normal(0,0.5));
		ρ²g=(log(rand(Uniform(0,1))) < -(map(z->logpdf(prior[:ρ²g],z),[ρ²gₚ,ρ²g])...))?ρ²gₚ : ρ²g
	end

	trace["ρ²g"][iter]=ρ²g
	trace["ρ²₀"][iter]=ρ²₀
	trace["ρ²₁"][iter]=ρ²₁
	trace["p"][iter]=p
	trace["łg"][iter]=łg
	trace["ł₀"][iter]=ł₀
	trace["μ₀bar"][iter]=μ₀bar
	trace["μ₁bar"][iter]=μ₁bar
	trace["ł₁"][iter]=ł₁
	#=trace["Λ"][iter]=Λ=#
	trace["λ₀"][iter]=map(x->λ₀(x),0:0.01:1)
	trace["λ₁"][iter]=map(x->λ₁(x),0:0.01:1)
	trace[:α][iter]=(map(x->map(y->x.α(y),0:0.01:1),ABtrials))
	trace[:g][iter]=(map(x->map(y->x.g(y),0:0.01:1),ABtrials))
	trace[:μg][iter]=map(x->x.μg,ABtrials)
	trace[:gᵧ][iter]=map(x->x.gᵧ,ABtrials)
end

plot(map(trial->trial.μg,trueABtrials),map(y->mean(map(x->trace[:μg][x][y],1:niter)),1:length(ABtrials)),linestyle="None",marker="o")

map(x->plot(trace["λ₀"][x],alpha=0.1,c="blue"),1:niter)
plot(map(x->Λ*Φ(μ₀(x)),0:0.01:1),c="red")

map(x->plot(trace["λ₁"][x],alpha=0.1,c="blue"),1:niter)
plot(map(x->Λ*Φ(μ₁(x)),0:0.01:1),c="red")

plot(trace["μ₀bar"])
plot(trace["ł₀"])
plot(trace["ρ²₀"])
plot(trace["ł₁"])
plot(trace["ρ²₁"])


plot(trace["p"])
plot(trace["łg"])
plot(trace["ρ²g"])
for intAB=1:length(ABtrials)
	figure(intAB)
	#=plot(map(x->ABtrials[intAB].α(x),0:0.01:1))=#
	for iter=1:niter
		plot(trace[:α][iter][intAB],alpha=0.01,c="red")
	end
	ylim(-0.1,1.1)
	savefig(string("plots/23/",intAB,".pdf"))
	close(intAB)
end
figure(44)
plot(map(x->λ₀(x),0:0.01:1))
plot(map(x->λ₁(x),0:0.01:1))
savefig("plots/22/44.pdf")
close(44)
figure(45)
plot(collect(1/2048:1/2048:1),kde(trace["p"]).density)
savefig("plots/22/45.pdf")
close(45)


plot(map(x->trueABtrials[1].α(x),0:0.01:1))
for iter=1:niter
	plot(map(x->trace[:α][iter][1](x),0:0.01:1),alpha=0.1,c="red")
end
