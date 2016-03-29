#include <random>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>
#include <map>
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Map;
using Eigen::JacobiSVD;



inline void Transition3d(Matrix3d & A, double d, double l){
	double d2=d*d;
	double d4=d*d*d*d;
	double l2=l*l;
	double l3=l*l*l;
	double l5=l*l*l*l*l;
	double emdl=exp(-d*l);
	//Transition Matrix
	A(0,0)=0.5*emdl*(2.0+2.0*d*l+d2*l2);
	A(0,1)=emdl*d*(1.0+d*l);
	A(0,2)=0.5*emdl*d2;
	A(1,0)=-0.5*emdl*d2*l3;
	A(1,1)=-emdl*(-1.0-d*l+d2*l2);
	A(1,2)=-0.5*emdl*d*(-2.0+d*l);
	A(2,0)=0.5*emdl*d*l3*(-2.0+d*l);
	A(2,1)=emdl*d*l2*(-3+d*l);
	A(2,2)=0.5*emdl*(2.0-4.0*d*l+d2*l2);
}

inline void RevTransition3d(Matrix3d & A, double d, double l){
	double d2=d*d;
	double d4=d*d*d*d;
	double l2=l*l;
	double l3=l*l*l;
	double l5=l*l*l*l*l;
	double emdl=exp(-d*l);
	//Transition Matrix
	A(0,0)=0.5*emdl*(2.0+2.0*d*l+d2*l2);
	A(0,1)=-1.0*emdl*d*(1.0+d*l);
	A(0,2)=0.5*emdl*d2;
	A(1,0)=0.5*emdl*d2*l3;
	A(1,1)=-emdl*(-1.0-d*l+d2*l2);
	A(1,2)=0.5*emdl*d*(-2.0+d*l);
	A(2,0)=0.5*emdl*d*l3*(-2.0+d*l);
	A(2,1)=-1.0*emdl*d*l2*(-3+d*l);
	A(2,2)=0.5*emdl*(2.0-4.0*d*l+d2*l2);
}

inline void Innovation3d(Matrix3d & Q, double d, double l){
	double d2=d*d;
	double d4=d*d*d*d;
	double l2=l*l;
	double l3=l*l*l;
	double l5=l*l*l*l*l;
	double em2dl=exp(-2.0*d*l);
	//double e2dl=exp(2.0*d*l);
	double q=(16.0*l5)/3.0;
	//Innovation Covariance Matrix
	Q(0,0)=(q*(3.0+em2dl*(-3.0-2.0*d*l*(3.0+d*l*(3.0+d*l*(2.0+d*l))))))/(16.0*l5);
	Q(0,1)=(1.0/8.0)*em2dl*q*d4;
	Q(0,2)=-((q*(-em2dl+1.0+2.0*d*l*em2dl*(-1.0+d*l*(-1.0+d*l*(-2.0+d*l)))))/(16.0*l3));
	Q(1,0)=Q(0,1);
	Q(1,1)=(q*(1.0+em2dl*(-1.0-2.0*d*l*(1.0+d*l*(-1.0+d*l)*(-1.0+d*l)))))/(16.0*l3);
	Q(1,2)=(1.0/8.0)*em2dl*q*d2*(-2.0+d*l)*(-2.0+d*l);
	Q(2,0)=Q(0,2);
	Q(2,1)=Q(1,2);
	Q(2,2)=(q*(3.0+em2dl*(-3.0-2.0*d*l*(-5.0+d*l*(11.0+d*l*(-6.0+d*l))))))/(16.0*l);
}


inline void RevInnovation3d(Matrix3d & Q, double d, double l){
	double d2=d*d;
	double d4=d*d*d*d;
	double l2=l*l;
	double l3=l*l*l;
	double l5=l*l*l*l*l;
	double em2dl=exp(-2.0*d*l);
	//double e2dl=exp(2.0*d*l);
	double q=(16.0*l5)/3.0;
	//Innovation Covariance Matrix
	Q(0,0)=(q*(3.0+em2dl*(-3.0-2.0*d*l*(3.0+d*l*(3.0+d*l*(2.0+d*l))))))/(16.0*l5);
	Q(0,1)=-1.0*(1.0/8.0)*em2dl*q*d4;
	Q(0,2)=-((q*(-em2dl+1.0+2.0*d*l*em2dl*(-1.0+d*l*(-1.0+d*l*(-2.0+d*l)))))/(16.0*l3));
	Q(1,0)=Q(0,1);
	Q(1,1)=(q*(1.0+em2dl*(-1.0-2.0*d*l*(1.0+d*l*(-1.0+d*l)*(-1.0+d*l)))))/(16.0*l3);
	Q(1,2)=-1.0*(1.0/8.0)*em2dl*q*d2*(-2.0+d*l)*(-2.0+d*l);
	Q(2,0)=Q(0,2);
	Q(2,1)=Q(1,2);
	Q(2,2)=(q*(3.0+em2dl*(-3.0-2.0*d*l*(-5.0+d*l*(11.0+d*l*(-6.0+d*l))))))/(16.0*l);
}
extern "C" void FFBS3d(double * y, double * xout, double * t, int n, double ls, double s2, double p2, double mu)
{
	Map<MatrixXd> x(xout,3,n);
	std::vector< Vector3d> m(n);
	std::vector< Matrix3d> M(n);
	std::vector< Matrix3d> AMAQ(n);
	std::vector< Matrix3d> A(n);
	Matrix3d C;
	Matrix3d Q;
	double l=sqrt(5.0)/ls;
	double l2=l*l;
	double l3=l*l*l;
	double l5=l*l*l*l*l;
	double q=(16.0*l5)/3.0;

	//Stationary Correlation Matrix
	C(0,0)=(3.0*q)/(16.0*l5);
	C(0,1)=0.0;
	C(0,2)=-(q/(16.0*l3));
	C(1,0)=0.0;
	C(1,1)=q/(16.0*l3);
	C(1,2)=0.0;
	C(2,0)=-(q/(16.0*l3));
	C(2,1)=0.0;
	C(2,2)=(3.0*q)/(16.0*l);
	m[0]=p2*C.col(0)*(y[0]-mu)/(s2+p2*C(0,0));
	M[0]=p2*C-p2*C.col(0)*C.row(0)*p2/(s2+p2*C(0,0));
	//Forward Filtering
	for(int i=1;i<n;++i){
		Innovation3d(Q,t[i]-t[i-1],l);
		Transition3d(A[i-1],t[i]-t[i-1],l);
		AMAQ[i-1]=A[i-1]*M[i-1]*A[i-1].transpose()+p2*Q;
		m[i]=A[i-1]*m[i-1]+AMAQ[i-1].col(0)*(y[i]-mu-A[i-1].row(0)*m[i-1])/(s2+AMAQ[i-1](0,0));
		M[i]=AMAQ[i-1]-(AMAQ[i-1].col(0)*AMAQ[i-1].row(0))/(s2+AMAQ[i-1](0,0));

	}

	//Backward Sampling
	Matrix3d S=0.5*(M[n-1]+M[n-1].transpose());
	Vector3d Z;
	std::random_device rd;std::mt19937 engine(rd());
	std::normal_distribution<double> N(0,1);
	for(int i=0;i<3;++i) Z(i)=N(engine);
	Eigen::LDLT<Matrix3d> ldlt;
	ldlt.compute(S);
	Matrix3d L(ldlt.matrixL());
	Vector3d D(ldlt.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);

	x.col(n-1)=m[n-1]+ldlt.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*Z;
	for(int i=(n-2); i>=0;--i){
		Eigen::FullPivLU<Matrix3d> lu(AMAQ[i]);
		x.col(i)=m[i]+M[i]*A[i].transpose()*lu.solve(x.col(i+1)-A[i]*m[i]);
		S=M[i]-M[i]*A[i].transpose()*lu.solve(A[i]*M[i]);
		S=0.5*(S+S.transpose());
		ldlt.compute(S);
		Matrix3d L(ldlt.matrixL());
		Vector3d D(ldlt.vectorD());
		for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
		for(int i=0;i<3;++i) Z(i)=N(engine);
		x.col(i)+=ldlt.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*Z;
	}
	Map<MatrixXd>( xout, x.rows(), x.cols() ) =  x;
}

extern "C" void ThreePoint(double * jx, double * jz, double t, double tl, double * jvl, double tr, double * jvr, double l){
	Map<Vector3d> x(jx);
	Map<Vector3d> z(jz);
	Map<Vector3d> vr(jvr);
	Map<Vector3d> vl(jvl);
	Eigen::LDLT<Matrix3d> Q;
	Matrix3d Qrl,Ql,Al,Ar;
	Innovation3d(Qrl,tr-tl,l);
	Innovation3d(Ql,t-tl,l);
	Transition3d(Al,t-tl,l);
	Transition3d(Ar,tr-t,l);
	Eigen::FullPivLU<Matrix3d> lu(Qrl);
	x=Al*vl+Ql*Ar.transpose()*lu.solve(vr-Ar*Al*vl);
	Q.compute(Ql-Ql*Ar.transpose()*lu.solve(Ar*Ql));
	Matrix3d L(Q.matrixL());
	Vector3d D(Q.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	x+=Q.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*z;
	Map<Vector3d>( jx,3,1) =   x;
}

inline Vector3d ThreePoint(std::mt19937 & engine, double t, double tl, Vector3d vl, double tr, Vector3d vr, double l){
	std::normal_distribution<double> N(0,1);
	Eigen::LDLT<Matrix3d> Q;
	Matrix3d Qrl,Ql,Al,Ar;
	Innovation3d(Qrl,tr-tl,l);
	Innovation3d(Ql,t-tl,l);
	Transition3d(Al,t-tl,l);
	Transition3d(Ar,tr-t,l);
	Eigen::FullPivLU<Matrix3d> lu(Qrl);
	Vector3d x=Al*vl+Ql*Ar.transpose()*lu.solve(vr-Ar*Al*vl);
	Q.compute(Ql-Ql*Ar.transpose()*lu.solve(Ar*Ql));
	Matrix3d L(Q.matrixL());
	Vector3d D(Q.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	Vector3d Z;
	for(int j=0;j<3;++j) Z(j)=N(engine);
	x+=Q.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*Z;
	return x;
}
inline Vector3d Backwards(std::mt19937 & engine,double t, double tr, Vector3d vr, double l){
	std::normal_distribution<double> N(0,1);
	Matrix3d Qr,Ar;
	RevTransition3d(Ar,tr-t,l);
	RevInnovation3d(Qr,tr-t,l);
	Vector3d x=Ar*vr;
	Eigen::LDLT<Matrix3d> S;
	S.compute(Qr);
	Matrix3d L(S.matrixL());
	Vector3d D(S.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	Vector3d Z;
	for(int j=0;j<3;++j) Z(j)=N(engine);
	x+=S.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*Z;
	return x;
}
extern "C" void Backwards(double * jx, double * jz,double t, double tr, double * jvr, double l){
	Map<Vector3d> x(jx);
	Map<Vector3d> z(jz);
	Map<Vector3d> vr(jvr);
	Matrix3d Qr,Ar;
	RevTransition3d(Ar,tr-t,l);
	RevInnovation3d(Qr,tr-t,l);
	x=Ar*vr;
	Eigen::LDLT<Matrix3d> S;
	S.compute(Qr);
	Matrix3d L(S.matrixL());
	Vector3d D(S.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	x+=S.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*z;
	Map<Vector3d>(jx,3,1)=x;
}
inline Vector3d Forwards(std::mt19937 & engine, double t, double tl, Vector3d vl, double l){
	std::normal_distribution<double> N(0,1);
	Matrix3d Al,Ql;
	Innovation3d(Ql,t-tl,l);
	Transition3d(Al,t-tl,l);
	Vector3d x=Al*vl;
	Eigen::LDLT<Matrix3d> S;
	S.compute(Ql);
	Matrix3d L(S.matrixL());
	Vector3d D(S.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	Vector3d Z;
	for(int j=0;j<3;++j) Z(j)=N(engine);
	x+=S.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*Z;
	return x;
}
extern "C" void Forwards(double * jx, double * jz, double t, double tl, double * jvl, double l){
	Map<Vector3d> x(jx);
	Map<Vector3d> z(jz);
	Map<Vector3d> vl(jvl);
	Matrix3d Al,Ql;
	Innovation3d(Ql,t-tl,l);
	Transition3d(Al,t-tl,l);
	x=Al*vl;
	Eigen::LDLT<Matrix3d> S;
	S.compute(Ql);
	Matrix3d L(S.matrixL());
	Vector3d D(S.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	x+=S.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*z;
	Map<Vector3d>(jx,3,1)=x;
}


extern "C" double LogDensity3d(double * y, double * t, int gamma, int n, double mu, double s2, double ls, double p2)
{
	if(ls<0) return -1.0/0.0;
	if(p2<0) return -1.0/0.0;

	if(gamma==1){
		std::vector< Vector3d> m(n);
		std::vector< Matrix3d> M(n);
		std::vector< Matrix3d> AMAQ(n);
		std::vector< Matrix3d> A(n);
		Matrix3d C;
		Matrix3d Q;
		double l=sqrt(5.0)/ls;
		double l2=l*l;
		double l3=l*l*l;
		double l5=l*l*l*l*l;
		double q=(16.0*l5)/3.0;

		//Stationary Correlation Matrix
		C(0,0)=(3.0*q)/(16.0*l5);
		C(0,1)=0.0;
		C(0,2)=-(q/(16.0*l3));
		C(1,0)=0.0;
		C(1,1)=q/(16.0*l3);
		C(1,2)=0.0;
		C(2,0)=-(q/(16.0*l3));
		C(2,1)=0.0;
		C(2,2)=(3.0*q)/(16.0*l);
		m[0]=p2*C.col(0)*(y[0]-mu)/(s2+p2*C(0,0));
		M[0]=p2*C-p2*C.col(0)*C.row(0)*p2/(s2+p2*C(0,0));
		//Forward Filtering
		for(int i=1;i<n;++i){
			Innovation3d(Q,t[i]-t[i-1],l);
			Transition3d(A[i-1],t[i]-t[i-1],l);
			AMAQ[i-1]=A[i-1]*M[i-1]*A[i-1].transpose()+p2*Q;
			m[i]=A[i-1]*m[i-1]+AMAQ[i-1].col(0)*(y[i]-mu-A[i-1].row(0)*m[i-1])/(s2+AMAQ[i-1](0,0));
			M[i]=AMAQ[i-1]-(AMAQ[i-1].col(0)*AMAQ[i-1].row(0))/(s2+AMAQ[i-1](0,0));

		}
		double logdensity=-0.5*log(2*M_PI*(s2+p2*C(0,0)))-0.5*(y[0]-mu)*(y[0]-mu)/(s2+p2*C(0,0));
		for(int i=1; i<n; ++i){
			logdensity+=-0.5*log(2*M_PI*(s2+AMAQ[i-1](0,0))) -0.5*(y[i]-(mu+A[i-1].row(0)*m[i-1]))*(y[i]-(mu+A[i-1].row(0)*m[i-1]))/(s2+AMAQ[i-1](0,0));
		}
		return(logdensity);
	}else{
		double logdensity=0;
		for(int i=0; i<n;++i){
			logdensity+=-0.5*log(2*M_PI*s2)-0.5*(y[i]-mu)*(y[i]-mu)/s2;

		}
		return(logdensity);
	}
}



extern "C" void Predict3d(double * xin, double * tc, int nc, double * xout, double * tp, int np, int gamma, double ls, double p2)
{
	double l=sqrt(5.0)/ls;
	Map<MatrixXd> xc(xin,3,nc);
	Map<MatrixXd> xp(xout,3,np);
	Matrix3d Al;
	Matrix3d Ar;
	Matrix3d Ql;
	Matrix3d QBE;
	Eigen::LDLT<Matrix3d> Q;
	Eigen::LDLT<Matrix3d> S;
	Vector3d Z;
	std::random_device rd;std::mt19937 engine(rd());
	std::normal_distribution<double> N(0,1);
	std::map<char,int> bound;
	std::map<char,int> sub;

	if((tc[0]<tp[0] && tp[0]< tc[nc-1] ) || (tp[0]<tc[0] && tc[0]<tp[np-1] )){ // [tc[0],tc[nc-1]] intersect [tp[0],tp[np-1]] not empty set
		//Find Sections
		int p=0;
		int c=0;
		while(c<nc && p<np){
			if(tp[p]<tc[c]){
				++p;
				if(p<np) if(tp[p]>tc[c]){ //tp[p] was before tc[c], now after tc[c], Just jumped over tc[c], subsequence bounded on the right by tc[c]
					bound['r']=c;
					sub['e']=p-1;
					if(bound.count('l')){ //If bounded on the left, then execute sandwich algorithm
						p=sub['b'];
						xp.col(p)=ThreePoint(engine,tp[p], tc[bound['l']], xc.col(bound['l']), tc[bound['r']], xc.col(bound['r']),l);
						++p;
						while(p<=sub['e']){
							xp.col(p)=ThreePoint(engine,tp[p], tp[p-1], xp.col(p-1), tc[bound['r']], xc.col(bound['r']),l);
							++p;
						}
						sub.empty();
						bound.empty();
					}else{ //Bounded only on the right
						--p;
						xp.col(p)=Backwards(engine,tp[p],tc[bound['r']],xc.col(bound['r']),l);
						--p;
						while(p>=0){
							xp.col(p)=Backwards(engine,tp[p],tp[p+1],xp.col(p+1),l);
							--p;
						}
						p=sub['e']+1;
						sub.empty();
						bound.empty();
					}
				}
			}else{
				++c;
				if(c<nc) if(tp[p]<tc[c]){ //tc[c] was before tp[p], now after tp[p], Just jumped over tp[p], subsequence bounded by tc[c-1],tc[c]
					bound['l']=c-1;
					bound['r']=c;
					sub['b']=p;

				}
			}
		}
		if(c==nc){ //xo...o
			bound['l']=nc-1;
			sub['b']=p;
			xp.col(p)=Forwards(engine,tp[p],tc[bound['l']],xc.col(bound['l']),l);
			++p;
			while(p<np){
				xp.col(p)=Forwards(engine,tp[p],tp[p-1],xp.col(p-1),l);
				++p;
			}
		}else{ // xo...ox...x
			p=sub['b'];
			bound['r']=c;
			xp.col(p)=ThreePoint(engine,tp[p], tc[bound['l']], xc.col(bound['l']), tc[bound['r']], xc.col(bound['r']),l);
			++p;
			while(p<np){
				xp.col(p)=ThreePoint(engine,tp[p], tp[p-1], xp.col(p-1), tc[bound['r']], xc.col(bound['r']),l);
				++p;
			}
		}
	}else if(tp[np-1]<tc[0]){ //All Predictive Points are on the left of the condition points
		int p=np-1;
		xp.col(p)=Backwards(engine,tp[p],tc[0],xc.col(0),l);
		--p;
		while(p>=0){
			xp.col(p)=Backwards(engine,tp[p],tp[p+1],xp.col(p+1),l);
			--p;
		}
	}else{ //All Condition points are on the left of the prediction points
		int p=0;
		xp.col(p)=Forwards(engine,tp[p],tc[nc-1],xc.col(nc-1),l);
		++p;
		while(p<np){
			xp.col(p)=Forwards(engine,tp[p],tp[p-1],xp.col(p-1),l);
			++p;
		}
	}

	Map<MatrixXd>( xout, xp.rows(), xp.cols() ) =   xp;
}


