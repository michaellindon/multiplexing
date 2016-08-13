importall Base.Random
importall Distributions
importall PyPlot
using DataStructures

srand(1)
include("classes.jl")
include("poissonpointprocess.jl")
include("statespace.jl")
include("eigen.jl")
#=include("hierdatagen.jl")=#
include("gaussianprocess.jl")
ρ²=5.0
ł=0.1
p=2
ν=p+0.5
ρ²=1.0
d=p+1
σ²=1.0

###LogMarginalLikelihood###
inputs=collect(0:0.01:1)
figure()
K=Kernel(inputs,inputs,ł,ρ²)
f=realization(statcov(ł,ρ²),transition(ł),innovation(ł,ρ²),inputs)
#=f=realization(ł,ρ²,inputs)=#
μ=10
y=[f[key][1]+μ+rand(Normal(0,√σ²)) for key in sort(collect(keys(f)))]
y=convert(Array{Float64},y)
logpdf(MvNormal(μ*ones(length(y)),σ²*eye(length(y))+K),y)
sslogdensity(SortedDict(Dict(zip(sort(collect(keys(f))),y))),1.0,μ,σ²,ł,ρ²)
logpdf(MvNormal(μ*ones(length(y)),σ²*eye(length(y))),y)
sslogdensity(SortedDict(Dict(zip(sort(collect(keys(f))),y))),0.0,μ,σ²,ł,ρ²)

σ²ₘ=1
mulogdensity(SortedDict(Dict(zip(sort(collect(keys(f))),y))),1.0,σ²,ł,ρ²,σ²ₘ)
logpdf(MvNormal(zeros(length(y)),σ²*eye(length(y))+K+σ²ₘ*ones(length(y),length(y))),y)
mulogdensity(SortedDict(Dict(zip(sort(collect(keys(f))),y))),0.0,σ²,ł,ρ²,σ²ₘ)
logpdf(MvNormal(zeros(length(y)),σ²*eye(length(y))+σ²ₘ*ones(length(y),length(y))),y)

###marginal mu full conditional###
σ²ₘ=0.1
mean0=3
oo=ones(length(inputs))
sqrt((1/(oo'*inv(eye(length(inputs))+K)*oo+1/σ²ₘ)[1])[1])
((1/(oo'*inv(eye(length(inputs))+K)*oo+1/σ²ₘ)[1])*(oo'*inv(eye(length(inputs))+K)*y+mean0/σ²ₘ))[1]
res=mu(SortedDict(Dict(zip(inputs,y))),1,σ²,ł,ρ²,σ²ₘ,mean0)
sqrt((1/(oo'*inv(eye(length(inputs)))*oo+1/σ²ₘ)[1])[1])
((1/(oo'*inv(eye(length(inputs)))*oo+1/σ²ₘ)[1])*(oo'*inv(eye(length(inputs)))*y+mean0/σ²ₘ))[1]
res=mu(SortedDict(Dict(zip(inputs,y))),0,σ²,ł,ρ²,σ²ₘ,mean0)


###REALIZATION###
inputs=collect(0:0.01:1)
figure()
K=Kernel(inputs,inputs,ł,ρ²)
for i=1:1000
	subplot(211)
	plot(rand(MvNormal(K)),c="blue")
	subplot(212)
	f=realization(ł,ρ²,inputs)
	plot(sort(collect(keys(f))),[f[key][1] for key in sort(collect(keys(f)))],c="blue")
end

###FFBS###	
σ²=0.03
xo=[0.0,0.3,0.5,0.8,0.9]
μ=100
yo=μ+rand(5)*10
xoyo=SortedDict(Dict(zip(xo,yo)))
xp=collect(-0.5:0.01:1.5)

ł=0.1
ρ²=1.5
Kₒₒ=Kernel(xo,xo,ł,ρ²)
Mean=Kₒₒ*\(Kₒₒ+σ²*eye(length(xo)),yo-μ)
Cov=Kₒₒ-Kₒₒ*\(Kₒₒ+σ²*eye(length(xo)),Kₒₒ')
figure()
niter=1000
subplot(211)
plot(collect(keys(xoyo)),collect(values(xoyo)),linestyle="none",marker="o",c="red")
subplot(212)
plot(collect(keys(xoyo)),collect(values(xoyo)),linestyle="none",marker="o",c="red")
for iter=1:niter
	subplot(211)
	plot(xo,rand(MvNormal(Mean,Cov)),c="blue",alpha=0.1)
	subplot(212)
	f=FFBS(xoyo,μ,σ²,ł,ρ²)
	plot(collect(keys(f)),[f[key][1] for key in collect(keys(f))],c="blue",alpha=0.1)
end


###FFBS and Prediction Functions###	
σ²=1.0
xo=[0.0,0.3,0.5,0.8,0.9]
μ=10.0
yo=μ+rand(5)*10
xoyo=SortedDict(Dict(zip(xo,yo)))
xp=collect(-0.5:0.01:1.5)

ł=0.1
ρ²=1.5
Kₒₒ=Kernel(xo,xo,ł,ρ²)
Kₚₚ=Kernel(xp,xp,ł,ρ²)
Kₚₒ=Kernel(xp,xo,ł,ρ²)
Mean=Kₚₒ*\(Kₒₒ+σ²*eye(length(xo)),yo-μ)
Cov=Kₚₚ-Kₚₒ*\(Kₒₒ+σ²*eye(length(xo)),Kₚₒ')
figure()
niter=1000
#=subplot(131)=#
#=plot(collect(keys(xoyo)),collect(values(xoyo)),linestyle="none",marker="o",c="red")=#
#=subplot(132)=#
#=plot(collect(keys(xoyo)),collect(values(xoyo)),linestyle="none",marker="o",c="red")=#
for iter=1:niter
	subplot(131)
	plot(xp,rand(MvNormal(Mean,Cov)),c="blue",alpha=0.1)
	subplot(132)
	f=FFBS(xoyo,μ,σ²,ł,ρ²)
	f=FFBS2(f,setdiff(xp,collect(keys(xoyo))),ł,ρ²)
	plot(collect(keys(f)),[f[key][1] for key in collect(keys(f))],c="blue",alpha=0.1)
	subplot(133)
	plot(xp,rand(MvNormal(Mean,Cov)),c="blue",alpha=0.01)
	plot(collect(keys(f)),[f[key][1] for key in collect(keys(f))],c="red",alpha=0.01)
end

subplot(211)
plot(xo,yo,linestyle="None",marker="o",c="red")
ylim(-5,12)
subplot(212)
plot(xo,yo,linestyle="None",marker="o",c="red")
ylim(-5,12)




###FFBS2###  ###Cannot compare against the non state space representation###	
xp=collect(-0.5:0.01:1.5)
xo=sort(rand(20))
#=realization(ł,ρ²,xo)=#
xc=sort(union(xo,xp))
figure()
niter=1000
for iter=1:niter
	subplot(211)
	full=realization(statcov(ł,ρ²),transition(ł),innovation(ł,ρ²),xc)
	plot(sort(collect(keys(full))),[full[key][1] for key in sort(collect(keys(full)))],c="blue")
	subplot(212)
	half=realization(statcov(ł,ρ²),transition(ł),innovation(ł,ρ²),xo)
	full=FFBS2(half,xc,ł,ρ²)
	plot(sort(collect(keys(full))),[full[key][1] for key in sort(collect(keys(full)))],c="blue")
end


###MORE FFBS2###
σ²=0.00003
xo=[0.0,0.3,0.5,0.8,0.9]
yo=rand(5)*10
xoyo=Dict(zip(xo,yo))
xp=collect(-0.5:0.01:1.5)
xc=sort(union(xo,xp))
figure()
niter=1000
for iter=1:niter
	subplot(211)
	f=FFBS(xoyo,xc,0.0,σ²,ł,ρ²)
	plot(sort(collect(keys(f))),[f[key][1] for key in sort(collect(keys(f)))],c="blue",alpha=0.1)
	subplot(212)
	f=FFBS(xoyo,xo,0.0,σ²,ł,ρ²)
	g=FFBS2(f,xc,ł,ρ²)
	plot(sort(collect(keys(g))),[g[key][1] for key in sort(collect(keys(g)))],c="blue",alpha=0.1)
end
subplot(211)
plot(xo,yo,linestyle="None",marker="o",c="red")
ylim(-5,12)
subplot(212)
plot(xo,yo,linestyle="None",marker="o",c="red")
ylim(-5,12)



























srand(1)
foo=innovation(0.1,ł)
yfoo=rand(d)
\(foo,yfoo)
chol(foo)
foo1=copy(foo)
yfoo1=copy(yfoo)
foo2=copy(foo)
yfoo2=copy(yfoo)
LAPACK.potrf!('L',foo1)
LAPACK.potrs!('L',foo1,yfoo1)
LAPACK.pstrf!('L',foo2,sqrt(eps(real(float(one(eltype(foo2)))))))
LAPACK.pstrs!('L',foo2,yfoo2,sqrt(eps(real(float(one(eltype(foo2)))))))


###NormalGeneration###
srand(1);
foo=rand(3,3);
foo=foo'*foo;
L=chol(foo,:L);
L=full(L);
z=rand(d);
res=zeros(Float64,3);
foo3=copy(foo); #Cholesky factor stored in place
LAPACK.potrf!('L',foo3);
#Note that the lower triangular of L and foo3 agree
L
foo3
#Method 1
L*z
#Method 2
BLAS.gemv!('N', 1.0, L, z, 0.0, res)
#Method 3
BLAS.trmv('L', 'N', 'N', foo3, z)


srand(1)                                                                                                                                                                                                          
a=rand(3,10)                                                                                                                                                                                                      
A=a'*a                                                                                                                                                                                                            
                                                                                                                                                                                                                   
rank(A)                                                                                                                                                                                                            
                                                                                                                                                                                                                   
M=A+eps()*eye(10)
F=cholfact(M,:L,Val{true},tol=eps()) 
L=F[:L]                                                                                                                                                                                                            
P=F[:P]                                                                                                                                                                                                            
P'*A*P-L*L'
P*L*L'*P'-A

F=cholfact(A+eps()*eye(10),:L,Val{true})
L=F[:L]                                                                                                                                                                                                            
P=F[:P]                                                                                                                                                                                                            
P'*A*P-L*L'


F=cholfact!(M,:L,Val{true})
F=cholfact(A+eps()*eye(10),:L,Val{false})
L=F[:L]                                                                                                                                                                                                            
P=F[:P]                                                                                                                                                                                                            
P'*A*P-L*L'



y=[y₀[key] for key in sort(collect(keys(y₀)))]
t=sort(collect(keys(y₀)))
x=Array{Array{Float64,2},1}()
m=Array{Array{Float64,2},1}()
M=Array{Array{Float64,2},1}()
A=Array{Array{Float64,2},1}()
AMAQ=Array{Array{Float64,2},1}()

x=[(zeros(d,1)) for i=1:length(y)]
m=[(zeros(d,1)) for i=1:length(y)]
M=[(zeros(d,d)) for i=1:length(y)]
A=[(zeros(d,d)) for i=1:length(y)]
AMAQ=[(zeros(d,d)) for i=1:length(y)]

x=[Mat(zeros(d)) for i=1:length(y)]
m=[Mat(zeros(d)) for i=1:length(y)]
M=[Mat(zeros(d,d)) for i=1:length(y)]
A=[Mat(zeros(d,d)) for i=1:length(y)]
AMAQ=[Mat(zeros(d,d)) for i=1:length(y)]
FFBS(y,x,t,m,M,A,AMAQ,μ₀,σ²,10,ρ²₀);

y=[y₀[key] for key in sort(collect(keys(y₀)))];
y=convert(Array{Float64},y);
t=sort(collect(keys(y₀)));
n=length(t);
@time ccall((:FFBS3d, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Int32,Float64,Float64,Float64,Float64), y,t,n,ł₀,σ²,ρ²₀,μ₀)

y=[y₀[key] for key in sort(collect(keys(y₀)))];
y=convert(Array{Float64},y);
t=sort(collect(keys(y₀)));
n=length(t);
@time ccall((:LogDensity3d, "./eigen.so"), Float64, (Ref{Cdouble},Ref{Cdouble},Int32,Int32,Float64,Float64,Float64,Float64), y,t,1,n,μ₀,σ²,ł₀,ρ²₀)


tout=rand(10000)
xout=Array{Float64}(3,10000)
xin=Array{Float64}(3,length(xin))

tin=sort(collect(keys(f₀)))
[xin[:,i]=f₀[tin[i]] for i=1:length(tin)]

srand(4)
xin=rand(3*100)
xout=ones(3*100)
tin=sort(rand(100))
tout=sort(rand(100))
tout[100]=1.0
tout[99]=0.999
tout[98]=0.998
ł₀=1.0
ρ²₀=1.0

g=Dict(zip(tin,xin))
@time gₚ=FFBS2(g,tout,ł₀,ρ²₀)

for i=1:100
ν=2.5
d=3
ł₀=1.0
ρ²₀=1.0
g=realization(ł₀,ρ²₀,collect(0.01:1/100:1))
xin=Array{Float64}(0)
tin=sort(collect(keys(g)))
tout=sort(rand(100))
xout=ones(length(tout)*3)
[xin=vcat(xin,vec(g[key])) for key in sort(collect(keys(g)))];
cxin=copy(xin)
#=@time ccall((:Predict3d, "./eigen.so"), Void, (Ref{Cdouble},Ref{Cdouble},Int32,Ref{Cdouble},Ref{Cdouble},Int32, Int32,Float64,Float64), xin,tin,length(tin),xout,tout,length(tout),1,ł₀,ρ²₀)=#
xout=FFBS2(g,tout,ł₀,ρ²₀)
xin=reshape(xin,3,length(tin))
xout=reshape(xout,3,length(tout))
plot(tin,vec(xin[1,:]),linestyle="none",marker=".",c="red")
plot(tout,vec(xout[1,:]),c="blue")
end

ν=2.5
d=3
ł₀=1.0
ρ²₀=1.0
μ=100.0
f=realization(ł₀,ρ²₀,collect(0.01:1/100:1))
y=SortedDict(Dict(zip(collect(keys(f)),(Float64)[μ+f[key][1] for key in collect(keys(f))])))
n=length(y)
mus=zeros(10000)
for i=1:10000
f=FFBS(y,μ,σ²,ł₀,ρ²₀)
μ=rand(Normal((n/σ²)*mean([y[key]-f[key][1] for key in collect(keys(f))])*(1/((n/σ²)+(1/σ²ₘ))),sqrt(1/((n/σ²)+(1/σ²ₘ)))))
mus[i]=μ
end


tout=sort(rand(100))
xout=FFBS2(g,tout,ł₀,ρ²₀)
plot(collect(keys(g)),[g[key][1] for key in collect(keys(g))],linestyle="none",marker=".",c="red")
plot(tout,[xout[key][1] for key in collect(keys(xout))],c="blue")



@async a=rand(Gamma(3,3),10000000)
@async b=rand(Gamma(3,3),10000000)
@async c=rand(Gamma(3,3),10000000)
@async d=rand(Gamma(3,3),10000000)
@async h=rand(Gamma(3,3),10000000)
@async f=rand(Gamma(3,3),10000000)
@async g=rand(Gamma(3,3),10000000)










@everywhere f(s,count)=(println("process id = $(myid()) s = $s count = $count");repeat(s,count))
pmap((a1,a2)->f(a1,a2),{"a","b","c"},{2,1,3})
[Ξ!(trial,μ₀,f₀,ł₀,ρ²₀,Ξₚ) for trial in Atrials];
pmap((args)->Ξ!(args...),[[trial,μ₀,f₀,ł₀,ρ²₀,Ξₚ] for trial in Atrials]) ; 

using Distributions

@everywhere function makenormals!(a,μ,σ)
	for i=1:length(a)
		a[i]=rand(Normal(μ,σ))
	end
end

myarrays=[zeros(Float64,10) for i=1:10]
pmap((args)->makenormals!(args...),[[myarrays[i],i,i] for i=1:length(myarrays)]) ; 
pmap((a1,a2,a3)->makenormals!(a1,a2,a3),myarrays,[0 for i=1:10],[1 for i=1:10])

@everywhere function testmod!(a,μ)
	for i=1:length(a)
		a[i]=i*μ
	end
	b=copy(a)
	return b
end

pmap((a1,a2)->testmod!(a1,a2),myarrays,[i for i=1:10])


using DistributedArrays

@everywhere type Atrial
	id::Int64
	Tobs::Float64
	y₀::Dict
	ξ₀ₐ::Dict
	ξ₀ᵣ::Dict
end

@everywhere function Atrial(id,ξ₀ₐ,Tobs) 
	return Atrial(id,Tobs,Dict(),ξ₀ₐ,Dict())
end

Atrials=[Atrial(i,Dict(zip(rand(10),rand(10)),),1) for i=1:10];
DAtrials=@DArray [trial for trial in Atrials];
LAtrials=convert(Array,DAtrials);


function test1()
ud=Dict{Int64,Int64}()
sizehint!(ud,1000)
for i=1:1000
	ud[i]=i
end
ud=SortedDict(ud)
end

function test2()
ud=SortedDict(Dict{Int64,Int64}())
for i=1:1000
	ud[i]=i
end
end



function test1()
ud=Set{Float64}()
sizehint!(ud,100000)
for i=1:100000
	push!(ud,i)
end
end

function test2()
ud=Array{Float64,1}()
sizehint!(ud,100000)
for i=1:100000
	push!(ud,i)
end
end




for trial in Atrials
	println(trial.id)
	plot(collect(keys(trial.y₀)),collect(values(trial.y₀)))
	foo=readline(STDIN)
end





ρ²₁=3
ł₁=0.2
foo=realization(statcov(ł₁,ρ²₁),transition(ł₁),innovation(ł₁,ρ²₁),collect(0:1/1000:1))
empty!(y₁)
[y₁[key]=foo[key][1] for key in collect(keys(foo))]
foob=map(z->sslogdensity(y₁,1,0,σ²,z,ρ²₁)+logpdf(prior[:ł₁],z),0.01:0.01:1)
foob=foob-maximum(foob)
plot(collect(0.01:0.01:1),exp(foob))
