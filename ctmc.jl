using Distributions #Random Variable Package
using PyPlot #Plotting Package

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
lam=[20,100]; #state rates
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

T=zeros(d+1,d+1,n+1); #T1,...,Tn+1 Transition Matrices
L=zeros(d+1,d+1,n+1); #L1,...,Ln+1
A=zeros(d+1,d+1,n+1); #A1,...,An+1
uflow=ones(n+1)
Gw=vcat(hcat(Q-diagm(lam),lam),zeros(d+1)') #Generator for Meta Process W
NU=zeros(Int64,length(t)) #For each Interval have number of 

prob=Array(Float64,d+1,n+1); #P1,...,Pn+1
S=Array(Float64,n+1); #S1,...,Sn+1
onevec=ones(d+1);
onevec[end]=0;
priorS0=ones(d+1)/d; #Prior on Zeroeth State S0
priorS0[end]=0;
probS0=copy(priorS0)

#Backwards Filtering#
A[:,:,n+1]=expm(Gw*(t[n+1]));
for k=n:-1:1
	L[:,:,k]=diagm(vcat(lam,0));
	T[:,:,k]=expm(Gw*(t[k]))
	A[:,:,k]=T[:,:,k]*L[:,:,k]*A[:,:,k+1]
	uflow[k]=maximum(A[:,:,k])
	A[:,:,k]/=uflow[k]
end

#Forward Sampling of S_i=X_t_i for i=0,1,...,n+1#
probS0=priorS0.*(A[:,:,1]*onevec)/(priorS0'A[:,:,1]*onevec)
S0=rand(Categorical(vec(probS0)))
prob[:,1]=T[S0,:,1]'.*diag(L[:,:,1]).*(A[:,:,2]*onevec)/((A[:,:,1]*onevec)[S0])
prob[:,1]/=uflow[1]
S[1]=rand(Categorical(prob[:,1]))
for k=2:n
	prob[:,k]=T[S[k-1],:,k]'.*diag(L[:,:,k]).*(A[:,:,k+1]*onevec)/((A[:,:,k]*onevec)[S[k-1]])
	prob[:,k]/=uflow[k]
	S[k]=rand(Categorical(prob[:,k]))
end
prob[:,n+1]=vec(A[S[n],:,n+1])
S[n+1]=rand(Categorical(prob[:,n+1]))

#Number of Dominating U-Events in intervals t_i=tp_i-tp_i-1 for i=1,2,...n+1, tp_0=0, tp_n+1=Tobs#
rho=maximum(-diag(Gw))
M=(1/rho)*Gw+eye(d+1)
U_times=Array(Float64,0)
for i=2:length(t)
	probvec=Array(Float64,0)
	r=0
	while sum(probvec)<0.9999
		push!(probvec,exp(-rho*t[i])*((rho*t[i])^r)*((M^r)[S[i-1],S[i]])/((factorial(BigInt(r)))*expm(Gw*t[i])[S[i-1],S[i]]))
		r=r+1
	end
	NU[i]=rand(Categorical(probvec/sum(probvec)))-1
	append!(U_times,rand(Uniform(tp[i-1],tp[i]),NU[i]))
end




plot([0:0.01:Tobs],map(state,[0:0.01:Tobs])); PyPlot.ylim([0,3])
PyPlot.plot(y_times, S[1:n], c="red", marker=".", linestyle="None")
