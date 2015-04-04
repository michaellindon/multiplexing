using Distributions #Random Variable Package
using PyPlot #Plotting Package

O=1 #Zero Indexing

srand(4); #Set Seed
d=2; #Dimensionality of State space

###Create Generator Matrix Q###
Q=rand(Beta(100,100),(d,d)); 
for r=1:d #Rows of Generator Matrix Must sum to 0
	Q[r,r]=-(sum(Q[r,1:d])-Q[r,r])
end
rho=-diag(Q) #Poisson Process Leaving Rates for CTMC states

###Create Transition Matrix P###
P=copy(Q) 
for r=1:d
	P[r,1:d]=P[r,1:d]/rho[r]
	P[r,r]=0
end

###Simulate CTMC###
Tobs=8
t=0
transition_times=Array(Float64,1)
s=Array(Int64,1)
transition_times[1]=0
s[1]=1
while(t<=Tobs)
	t=t+rand(Exponential(rho[1]))
	if(t<Tobs)
		push!(transition_times,t)
		push!(s,rand(Categorical(vec(P[s[end],1:d]))))
	end
end
state(arg)=s[maximum(find( x->(x <= arg), transition_times))]
plot([0:0.01:Tobs],map(state,[0:0.01:Tobs])); PyPlot.ylim([0,3])

###Simulate Spike Train###
###Simulate heterogeneous pprocess from dominating homogeneous poisson process by thinning###
t=0
lam=[20.0,100.0]; #state rates
y_times=Array(Float64,0)
while(t<=Tobs)
	u=rand(Uniform(0,1))
	t=t-log(u)/maximum(lam)
	if(t<=Tobs)
		if(rand(Uniform(0,1))<=lam[state(t)]/maximum(lam)) 
			push!(y_times,t)
		end
	end
end
plot([0:0.01:Tobs],map(state,[0:0.01:Tobs])); PyPlot.ylim([0,3])
PyPlot.plot(y_times, map(state,y_times), c="red", marker=".", linestyle="None")


n=length(y_times); #Number of Y events
tp=copy(y_times); #t'1,...,t'n
unshift!(tp,0); #set t'_0=0
push!(tp,Tobs); #set t_n+1=Tobs
t=tp[2:end]-tp[1:(end-1)]; #t1,t2,...,#tn+1

T=zeros(d+1,d+1,n+2); #T1,...,Tn+1 Transition Matrices
L=zeros(d+1,d+1,n+2); #L1,...,Ln+1
A=zeros(d+1,d+1,n+2); #A1,...,An+1
uflow=ones(n+2)
Gw=vcat(hcat(Q-diagm(lam),lam),zeros(d+1)') #Generator for Meta Process W
NU=zeros(Int64,length(t)) #For each Interval have number of 

prob=Array(Float64,d+1,n+2); #P1,...,Pn+1
S=Array(Float64,n+2); #S1,...,Sn+1
onevec=ones(d+1);
onevec[end]=0;
priorS0=ones(d+1)/d; #Prior on Zeroeth State S0
priorS0[end]=0;

niter=1000
for iter=1:niter

	####Sample States at Observed Times###

	#Backwards Filtering#
	A[:,:,end]=expm(Gw*(t[end]));
	for k=n:-1:0
		L[:,:,O+k]=diagm(vcat(lam,0));
		T[:,:,O+k]=expm(Gw*(t[O+k]))
		A[:,:,O+k]=T[:,:,O+k]*L[:,:,O+k]*A[:,:,O+k+1]
		uflow[O+k]=maximum(A[:,:,O+k])
		A[:,:,O+k]/=uflow[O+k]
	end

	#Forward Sampling of S_i=X_t_i for i=0,1,...,n+1#
	prob[:,O]=priorS0.*(A[:,:,O+1]*onevec)/(priorS0'A[:,:,O+1]*onevec)
	S[O]=rand(Categorical(vec(prob[:,O])))
	for k=1:n
		prob[:,O+k]=T[S[O+k-1],:,O+k]'.*diag(L[:,:,O+k]).*(A[:,:,O+k+1]*onevec)/((A[:,:,O+k]*onevec)[S[O+k-1]])
		prob[:,O+k]/=uflow[O+k]
		S[O+k]=rand(Categorical(prob[:,O+k]))
	end
	prob[:,O+n+1]=vec(A[S[O+n],:,O+n+1])
	S[O+n+1]=rand(Categorical(prob[:,O+n+1]))

	###Sample CTMC Conditional on Observed States###
	#Number of Dominating U-Events in intervals t_i=tp_i-tp_i-1 for i=1,2,...n+1, tp_0=0, tp_n+1=Tobs#
	rho=maximum(-diag(Gw))
	M=(1/rho)*Gw+eye(d+1)
	U_times=Array(Float64,0)
	for i=1:length(t)
		probvec=Array(Float64,0)
		r=0
		while sum(probvec)<0.9999
			push!(probvec,exp(-rho*t[i])*((rho*t[i])^r)*((M^r)[S[O+i-1],S[O+i]])/((factorial(BigInt(r)))*expm(Gw*t[i])[S[O+i-1],S[O+i]]))
			r=r+1
		end
		NU[i]=rand(Categorical(probvec/sum(probvec)))-1
		append!(U_times,sort(rand(Uniform(tp[O+i-1],tp[O+i]),NU[i])))
	end

	US=Array(Int64,0);
	push!(US,S[O]);
	for interval=1:(n+1)
		St=S[O+interval]
		r=NU[interval]
		for j=1:r
			probvec=M[US[end],:]'.*(M^(r-j))[:,St]/(M^(r-j+1))[US[end],St]
			push!(US,rand(Categorical(vec(probvec))));
		end
		push!(US,St);
	end
	all_times=sort(vcat(tp,U_times))
	intervals=all_times[2:end]-all_times[1:end-1]

	count(x->x==2,S)

	for i=1:d
	lam[i]=rand(Gamma(0.001+count(x->x==i,S),1/(0.001+sum(intervals[find(x->x==i,US[1:end-1])]))))
	end

	

	PyPlot.plot(all_times, US+3, c="red", alpha=0.01)

end

plot([0:0.001:Tobs],map(state,[0:0.001:Tobs])); PyPlot.ylim([-1,7])
PyPlot.plot(all_times, US, c="red", alpha=0.1)
PyPlot.plot(y_times,zeros(length(y_times)),c="blue",marker=".",linestyle="None")
PyPlot.plot(U_times, 1.5*ones(length(U_times)), c="red", marker=".", linestyle="None")
PyPlot.plot(y_times, S[1:n], c="red", marker=".", linestyle="None")


