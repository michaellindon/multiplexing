using Distributions #Random Variable Package
using PyPlot #Plotting Package


srand(10); #Set Seed
d=3; #Dimensionality of State space

###Create Generator Matrix A###
A=rand(Beta(100,100),(d,d)); 
for r=1:d #Rows of Generator Matrix Must sum to 0
	A[r,r]=-(sum(A[r,1:d])-A[r,r])
end

Omega=2*maximum(abs(A))
B=eye(d)+A/Omega
tobs=8
n=rand(Poisson(tobs*Omega))
T=sort(rand(Uniform(0,tobs),n))

S=Array(Int64,1)
S[1]=1
T_thin=Array(Float64,0)
T_tran=Array(Float64,0)
S_thin=Array(Int64,0)
S_tran=Array(Int64,0)
for i=1:n
	push!(S,rand(Categorical(vec(B[S[end],:]))))
	if(S[end]==S[end-1])
		push!(S_thin,S[end])
		push!(T_thin,T[i])
	else 
		push!(S_tran,S[end])
		push!(T_tran,T[i])
	end
end
unshift!(T,0); #set t'_0=0
push!(T,tobs); #set t_n+1=tobs
push!(S,S[end])


###Simulate CTMC###
state(arg)=S[maximum(find( x->(x <= arg), T))]
plot([0:0.001:tobs],map(state,[0:0.001:tobs])); PyPlot.ylim([-0.5,3.5])
plot(T_thin,S_thin,marker="o",linestyle="None",c="white")
plot(T_tran,S_tran,marker="o",linestyle="None",c="black")

###Simulate Spike Train###
###Simulate heterogeneous pprocess from dominating homogeneous poisson process by thinning###
t=0
lam=[1.0,300.0,600]; #state rates
lam_dom=600
n_dom=rand(Poisson(lam_dom*tobs))
o_dom=sort(rand(Uniform(0,tobs),n_dom))
o=Array(Float64,0)
o_thin=Array(Float64,0)
for i=1:n_dom
	u=rand(Uniform(0,1))
	if(u<=lam[state(o_dom[i])]/lam_dom)
		#Accept
		push!(o,o_dom[i])
	else 
		#Thin
		push!(o_thin,o_dom[i])
	end
end

PyPlot.plot(o,0.5*ones(length(o)),c="red",marker="|",linestyle="None")
PyPlot.plot(o_thin,-0.5*ones(length(o_thin)),marker="|",linestyle="None",c="green")



prior_initial=ones(d)/d; #Prior on Zeroeth State S0


niter=100
for iter=1:niter

	U=Array(Float64,0);
	for i=1:length(T)-1
		nu=rand(Poisson(Omega-abs(A[S[i],S[i]])))
		ru=rand(Uniform(T[i],T[i+1]),nu)
		append!(U,ru)
	end
	W=sort(union(T,U))
	nw=length(W)
	V=Array(Int64,nw)
	prob=Array(Float64,d,nw);

	R=zeros(d,d,nw); #T1,...,Tn+1 Transition Matrices
	logL=zeros(d,nw-1); #L1,...,Ln+1
	logH=zeros(d,nw); #A1,...,An+1
	uflow=ones(nw)
	exponent=Array(Float64,d)

	#Backwards Filtering#
	logH[:,end]=zeros(d)
	for k=(nw-1):-1:1
		for i=1:d
			logL[i,k]=count(x->(x>=W[k]&&x<W[k+1]),o)*log(lam[i]) - lam[i]*(W[k+1]-W[k]);
		end
		R[:,:,k]=eye(d)+A/Omega

		for i=1:d
			for j=1:d
				exponent[j]=log(R[i,j,k])+logL[j,k]+logH[j,k+1];
			end
			M=maximum(exponent)
			logH[i,k]=M+log(sum(exp(exponent-M)))
		end
		#uflow[k]=maximum(logH[:,k])
		#logH[:,k]/=uflow[k]
	end

	S=Array(Int64,0)
	T=zeros(Float64,1)
	#Forward Sampling of S_i=X_t_i for i=0,1,...,n+1#

	prob[:,1]=zeros(d)
	for i=1:d
		for j=1:d
			prob[i,1]+=exp(logL[i,1]+log(R[j,i,1])+logH[i,2]-logH[j,1]+log(prior_initial[i]))
		end
	end
#	prob[:,1]/=uflow[1]

	V[1]=rand(Categorical(vec(prob[:,1])))
	push!(S,V[1])
	for k=2:(nw-1)
		prob[:,k]=exp(log(R[V[k-1],:,k])' + logL[:,k] + logH[:,k+1] - logH[[V[k-1]],k][1])
#		prob[:,k]/=uflow[k]
		V[k]=rand(Categorical(prob[:,k]))
		if(V[k]!=V[k-1])
			push!(S,V[k])
			push!(T,W[k])
		end
	end
	R[:,:,nw]=eye(d)+A/Omega
	prob[:,nw]=vec(R[V[nw-1],:,nw])
	V[nw]=rand(Categorical(prob[:,nw]))
	push!(T,tobs)
	push!(S,V[nw])

	state(arg)=S[maximum(find( x->(x <= arg), T))]
	plot([0:0.001:tobs],map(state,[0:0.001:tobs]),c="red",alpha=0.1); PyPlot.ylim([-0.5,3.5])

end

plot([0:0.001:tobs],map(state,[0:0.001:tobs])); PyPlot.ylim([-1,7])
PyPlot.plot(O,zeros(length(O)),c="blue",marker="|",linestyle="None")
PyPlot.plot(all_times, US, c="red", marker=".",linestyle="None")
PyPlot.plot(U_times, 1.5*ones(length(U_times)), c="red", marker=".", linestyle="None")
PyPlot.plot(O, S[1:n], c="red", marker=".", linestyle="None")


