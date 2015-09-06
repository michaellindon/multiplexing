using Distributions
using PyPlot

srand(4)

x= [0:0.1:8]
nx=length(x)
rho2=1
psi2=0.1
C=Array(Float64,nx,nx)
for i=1:nx
	for j=1:nx
		C[i,j]=rho2*exp(-psi2*(x[i]-x[j])^2)
	end
end

L=chol(C+0.0001*eye(nx))'
s2=1
z=rand(Normal(0,1),nx)
g=sqrt(s2)*L*z
lam_star=100
intensity1=lam_star*cdf(Normal(0,1),g)
plot(x,intensity1)
srand(8)
for i=1:nx
	for j=1:nx
		C[i,j]=rho2*exp(-psi2*(x[i]-x[j])^2)
	end
end
L=chol(C+0.0001*eye(nx))'
s2=1
z=rand(Normal(0,1),nx)
g=sqrt(s2)*L*z
lam_star=100
intensity2=lam_star*cdf(Normal(0,1),g)
plot(x,intensity2)


function intensity(time,lamda)
	if(lamda==1)
		return intensity1[maximum(find( y->(y <= time), x))]
	else
		return intensity2[maximum(find( y->(y <= time), x))]
	end
end

set.seed(1)
###Create Generator Matrix A###
d=2; #Dimensionality of State space
A=5*rand(Beta(100,100),(d,d)); 
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
plot([0:0.001:tobs],map(state,[0:0.001:tobs]));# PyPlot.ylim([-0.5,3.5])
#plot(T_thin,S_thin,marker="o",linestyle="None",c="white")
plot(T_tran,S_tran,marker="o",linestyle="None",c="black")



tobs=8
lam_dom=lam_star
no=rand(Poisson(lam_dom*tobs))
o=sort(rand(Uniform(0,tobs),no))
oacc=Array(Float64,0)
othin=Array(Float64,0)
binary=Array(Int64,0)
for i=1:no
	u=rand(Uniform(0,1))
	if(u<=intensity(o[i],state(o[i]))/lam_dom)
		#Accept
		push!(oacc,o[i])
		push!(binary,1)
	else 
		#Thin
		push!(othin,o[i])
		push!(binary,0)
	end
end
na=length(oacc)

PyPlot.plot(oacc,0.5*ones(length(oacc)),c="red",marker="|",linestyle="None")
#PyPlot.plot(othin,-0.5*ones(length(othin)),marker="|",linestyle="None",c="green")

ix1=find(x->x==1,map(state,o))
binary1=binary[ix1]
o1=o[ix1]
no1=length(o1)
C1=Array(Float64,no1,no1)
for i=1:no1
	for j=1:no1
		C1[i,j]=rho2*exp(-psi2*(o1[i]-o1[j])^2)
	end
end

C1=C1+0.0001*eye(no1)
L1=chol(C1)'
z=rand(Normal(0,1),no1)
g1=sqrt(s2)*L1*z

ix2=find(x->x==2,map(state,o))
binary2=binary[ix2]
o2=o[ix2]
no2=length(o2)
C2=Array(Float64,no2,no2)
for i=1:no2
	for j=1:no2
		C2[i,j]=rho2*exp(-psi2*(o2[i]-o2[j])^2)
	end
end

C2=C2+0.0001*eye(no2)
L2=chol(C2)'
z=rand(Normal(0,1),no2)
g2=sqrt(s2)*L2*z

g=Array(Float64,0)

i1=1
i2=1
while i1<=no1 && i2<=no2
	if(o1[i1]<o2[i2])
		push!(g,g1[i1])
		i1=i1+1
	else
		push!(g,g2[i2])
		i2=i2+1
	end
end
while i1<=no1
	push!(g,g1[i1])
	i1=i1+1
end
while i2<=no2
	push!(g,g2[i2])
	i2=i2+1
end


prior_initial=ones(d)/d
for iter=1:10
	Z1=Array(Float64,no1);
	for i=1:no1
		if(binary1[i]==1)
			Z1[i]=rand(Truncated(Normal(g1[i],1), 0, Inf),1)[1]
		else
			Z1[i]=rand(Truncated(Normal(g1[i],1), -Inf, 0),1)[1]
		end
	end
	C1=Array(Float64,no1,no1);
	for i=1:no1
		for j=1:no1
			C1[i,j]=rho2*exp(-psi2*(o1[i]-o1[j])^2)
		end
	end
	C1=C1+0.0001*eye(no1);
	M1=C1*inv(eye(no1)+C1);
	L1=chol(M1)';
	g1=M1*Z1+L1*rand(Normal(0,1),no1);
	plot(o1,lam_star*cdf(Normal(0,1),g1),c="grey",alpha=0.1);
	Z2=Array(Float64,no2);
	for i=1:no2
		if(binary2[i]==1)
			Z2[i]=rand(Truncated(Normal(g2[i],1), 0, Inf),1)[1]
		else
			Z2[i]=rand(Truncated(Normal(g2[i],1), -Inf, 0),1)[1]
		end
	end
	C2=Array(Float64,no2,no2);
	for i=1:no2
		for j=1:no2
			C2[i,j]=rho2*exp(-psi2*(o2[i]-o2[j])^2)
		end
	end
	C2=C2+0.0001*eye(no2);
	M2=C2*inv(eye(no2)+C2);
	L2=chol(M2)';
	g2=M2*Z2+L2*rand(Normal(0,1),no2);
	plot(o2,lam_star*cdf(Normal(0,2),g2),c="grey",alpha=0.1);
	n_prop=rand(Poisson(lam_dom*tobs));
	oprop=sort(rand(Uniform(0,tobs),n_prop));
	ix1p=find(x->x==1,map(state,oprop))
	oprop1=oprop[ix1p]
	np1=length(oprop1)
	C1pp=Array(Float64,np1,np1);
	for i=1:np1
		for j=1:np1
			C1pp[i,j]=rho2*exp(-psi2*(oprop1[i]-oprop1[j])^2)
		end
	end
	C1po=Array(Float64,np1,no1);
	for i=1:np1
		for j=1:no1
			C1po[i,j]=rho2*exp(-psi2*(oprop1[i]-o1[j])^2)
		end
	end
	g1prop=Array(Float64,np1);
	gprop1=C1po*\(C1,g1)+chol(C1pp-C1po*\(C1,C1po')+0.0001*eye(np1))'*rand(Normal(0,1),np1);
	g1thin=Array(Float64,0);
	o1thin=Array(Float64,0);
	for i=1:np1
		u=rand(Uniform(0,1))
		if(u<=lam_star*cdf(Normal(0,1),g1prop[i])/lam_dom)
			#Accept
		else 
			#Thin
			push!(o1thin,oprop1[i])
			push!(g1thin,gprop1[i])
		end
	end
	ix2p=find(x->x==2,map(state,oprop))
	oprop2=oprop[ix2p]
	np2=length(oprop2)
	C2pp=Array(Float64,np2,np2);
	for i=1:np2
		for j=1:np2
			C2pp[i,j]=rho2*exp(-psi2*(oprop2[i]-oprop2[j])^2)
		end
	end
	C2po=Array(Float64,np2,no2);
	for i=1:np2
		for j=1:no2
			C2po[i,j]=rho2*exp(-psi2*(oprop2[i]-o2[j])^2)
		end
	end
	g2prop=Array(Float64,np2);
	gprop2=C2po*\(C2,g2)+chol(C2pp-C2po*\(C2,C2po')+0.0002*eye(np2))'*rand(Normal(0,1),np2);
	g2thin=Array(Float64,0);
	o2thin=Array(Float64,0);
	for i=1:np2
		u=rand(Uniform(0,1))
		if(u<=lam_star*cdf(Normal(0,1),g2prop[i])/lam_dom)
			#Accept
		else 
			#Thin
			push!(o2thin,oprop2[i])
			push!(g2thin,gprop2[i])
		end
	end
	#Collect thinned o's
	nt=length(o1thin)+length(o2thin)
	gthin=Array(Float64,0)
	othin=Array(Float64,0)
	i1=1
	i2=1
	while i1<=length(o1thin) && i2<=length(o2thin)
		if(o1thin[i1]<o2thin[i2])
			push!(othin,o1thin[i1])
			push!(gthin,g1thin[i1])
			i1=i1+1
		else
			push!(othin,o2thin[i2])
			push!(gthin,g2thin[i2])
			i2=i2+1
		end
	end
	while i1<=length(o1thin)
		push!(othin,o1thin[i1])
		push!(gthin,g1thin[i1])
		i1=i1+1
	end
	while i2<=length(o2thin)
		push!(othin,o2thin[i2])
		push!(gthin,g2thin[i2])
		i2=i2+1
	end
	no=nt+na;
	o=Array(Float64,0);
	gacc=g[find(x->x==1,binary)];
	g=Array(Float64,0);
	binary=Array(Int64,0);
	ai=1;
	ti=1;
	while ti<=nt && ai<=na
		if(othin[ti]<oacc[ai])
			push!(o,othin[ti])
			push!(binary,0)
			push!(g,gthin[ti])
			ti=ti+1
		else
			push!(o,oacc[ai])
			push!(binary,1)
			push!(g,gacc[ai])
			ai=ai+1
		end
	end
	while ti<=nt
		push!(o,othin[ti])
		push!(g,gthin[ti])
		push!(binary,0)
		ti=ti+1
	end
	while ai<=na
		push!(o,oacc[ai])
		push!(binary,1)
		push!(g,gacc[ai])
		ai=ai+1
	end

	g1=g[find(x->x==1,map(state,o))]
	g2=g[find(x->x==2,map(state,o))]
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

	ix1=find(x->x==1,map(state,o))
	binary1=binary[ix1]
	o1=o[ix1]
	no1=length(o1)
	C1=Array(Float64,no1,no1)
	for i=1:no1
		for j=1:no1
			C1[i,j]=rho2*exp(-psi2*(o1[i]-o1[j])^2)
		end
	end
	C1=C1+0.0001*eye(no1)
	ix2=find(x->x==2,map(state,o))
	binary2=binary[ix2]
	o2=o[ix2]
	no2=length(o2)
	C2=Array(Float64,no2,no2)
	for i=1:no2
		for j=1:no2
			C2[i,j]=rho2*exp(-psi2*(o2[i]-o2[j])^2)
		end
	end
	C2=C2+0.0001*eye(no2)
	#Backwards Filtering#
	logH[:,end]=zeros(d)
	for k=(nw-1):-1:1
		accepted=oacc[find(x->(x>=W[k]&&x<W[k+1]),oacc)]
		ga=gacc[find(x->(x>=W[k]&&x<W[k+1]),gacc)]
		oa=oacc[find(x->(x>=W[k]&&x<W[k+1]),oacc)]
		na=length(accepted)
		thinned=othin[find(x->(x>=W[k]&&x<W[k+1]),othin)]
		gt=gthin[find(x->(x>=W[k]&&x<W[k+1]),gthin)]
		ot=othin[find(x->(x>=W[k]&&x<W[k+1]),othin)]
		nt=length(thinned)
		C1a=Array(Float64,na,na);
		for i=1:na
			for j=1:na
				C1a[i,j]=rho2*exp(-psi2*(oa[i]-oa[j])^2)
			end
		end
		C2a=Array(Float64,na,na);
		for i=1:na
			for j=1:na
				C2a[i,j]=rho2*exp(-psi2*(oa[i]-oa[j])^2)
			end
		end
		C1t=Array(Float64,nt,nt);
		for i=1:nt
			for j=1:nt
				C1t[i,j]=rho2*exp(-psi2*(ot[i]-ot[j])^2)
			end
		end
		C2t=Array(Float64,nt,nt);
		for i=1:nt
			for j=1:nt
				C2t[i,j]=rho2*exp(-psi2*(ot[i]-ot[j])^2)
			end
		end
		C1ao=Array(Float64,na,no1);
		for i=1:na
			for j=1:no1
				C1ao[i,j]=rho2*exp(-psi2*(oa[i]-o1[j])^2)
			end
		end
		C2ao=Array(Float64,na,no2);
		for i=1:na
			for j=1:no2
				C2ao[i,j]=rho2*exp(-psi2*(oa[i]-o2[j])^2)
			end
		end
		C1to=Array(Float64,nt,no1);
		for i=1:nt
			for j=1:no1
				C1to[i,j]=rho2*exp(-psi2*(ot[i]-o1[j])^2)
			end
		end
		C2to=Array(Float64,nt,no2);
		for i=1:nt
			for j=1:no2
				C2to[i,j]=rho2*exp(-psi2*(ot[i]-o2[j])^2)
			end
		end
		for i=1:d
			if(i==1)
				if(na!=0)
					ga=C1ao*\(C1,g1)+chol(C1a-C1ao*\(C1,C1ao')+0.002*eye(na))'*rand(Normal(0,1),na);
				end
				if(nt!=0)
					gt=C1to*\(C1,g1)+chol(C1t-C1to*\(C1,C1to')+0.002*eye(nt))'*rand(Normal(0,1),nt);
				end
			end
			if(i==2)
				if(na!=0)
					ga=C2ao*\(C2,g2)+chol(C2a-C2ao*\(C2,C2ao')+0.002*eye(na))'*rand(Normal(0,1),na);
				end
				if(nt!=0)
					gt=C2to*\(C2,g2)+chol(C2t-C2to*\(C2,C2to')+0.002*eye(nt))'*rand(Normal(0,1),nt);
				end
			end
			likeobs=0
			if(na!=0)
				likeobs=sum(log(cdf(Normal(0,1),ga)))
			end
			likethin=0
			if(nt!=0)
				likethin=sum(log(1-log(cdf(Normal(0,1),gt))))
			end
			logL[i,k]=(na+nt)*log(lam_dom)-lam_dom*(W[k+1]-W[k])+likeobs+likethin;
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
end
