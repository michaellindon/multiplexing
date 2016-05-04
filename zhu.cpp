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

inline void priorVariance(Matrix3d & C){
	//Transition Matrix
	C(0,0)=1000.0;
	C(0,1)=0.0;
	C(0,2)=0.0;
	C(1,0)=0.0;
	C(1,1)=1000.0;
	C(1,2)=0.0;
	C(2,0)=0.0;
	C(2,1)=0.0;
	C(2,2)=1000.0;
}


inline void Transition3d(Matrix3d & A, double d, bool approx){
	//Transition Matrix
	if(!approx){
	A(0,0)=1;
	A(0,1)=d;
	A(0,2)=0.5*d*d;
	A(1,0)=0;
	A(1,1)=1;
	A(1,2)=d;
	A(2,0)=0;
	A(2,1)=0;
	A(2,2)=1;
	}else{
	A(0,0)=1;
	A(0,1)=d;
	A(0,2)=0;
	A(1,0)=0;
	A(1,1)=1;
	A(1,2)=d;
	A(2,0)=0;
	A(2,1)=0;
	A(2,2)=1;
	}
}

inline void RevTransition3d(Matrix3d & A, double d){
	//Transition Matrix
	A(0,0)=1;
	A(0,1)=-d;
	A(0,2)=0.5*d*d;
	A(1,0)=-0;
	A(1,1)=1;
	A(1,2)=-d;
	A(2,0)=0;
	A(2,1)=-0;
	A(2,2)=1;
}


inline void Innovation3d(Matrix3d & Q, double d, double s2u, double s2a,bool approx){
	if(!approx){
	double d2=d*d;
	double d3=d2*d;
	double d4=d3*d;
	double d5=d4*d;
	Q(0,0)=(1/20)*d5*s2a+(1/3)*d3*s2u;
	Q(0,1)=(1/8)*d4*s2a+(1/2)*d2*s2u;
	Q(0,2)=(1/6)*d3*s2a; 
	Q(1,0)=(1/8)*d4*s2a+(1/2)*d2*s2u; 
	Q(1,1)=(1/3)*d3*s2a+d*s2u;
	Q(1,2)=(1/2)*d2*s2a;
	Q(2,0)=(1/6)*d3*s2a;
	Q(2,1)=(1/2)*d2*s2a;
	Q(2,2)=d*s2a;
	}else{
	Q(0,0)=0;
	Q(0,1)=0;
	Q(0,2)=0;
	Q(1,0)=0;
	Q(1,1)=d*s2u;
	Q(1,2)=0;
	Q(2,0)=0;
	Q(2,1)=0;
	Q(2,2)=d*s2a;
	}	
}


inline void RevInnovation3d(Matrix3d & Q, double d, double s2u, double s2a){
	double d2=d*d;
	double d3=d2*d;
	double d4=d3*d;
	double d5=d4*d;
	Q(0,0)=(1/20)*d5*s2a+(1/3)*d3*s2u;
	Q(0,1)=-(1/8)*d4*s2a+(1/2)*d2*s2u;
	Q(0,2)=(1/6)*d3*s2a; 
	Q(1,0)=Q(0,1);
	Q(1,1)=(1/3)*d3*s2a+d*s2u;
	Q(1,2)=-(1/2)*d2*s2a;
	Q(2,0)=Q(0,2);
	Q(2,1)=Q(1,2);
	Q(2,2)=d*s2a; 
}

extern "C" void FFBS3d(double * y, double * xout, double * t, int n, double s2, double s2u, double s2a, double p2, double mu, double * jz, bool approx)
{
	if(approx){
		std::cout << "true" << std::endl;
	}else{
		std::cout << "false" << std::endl;
	}
	Map<MatrixXd> x(xout,3,n);
	Map<MatrixXd> Z(jz,3,n);
	std::vector< Vector3d> m(n);
	std::vector< Matrix3d> M(n);
	std::vector< Matrix3d> AMAQ(n);
	std::vector< Matrix3d> A(n);
	Matrix3d C;
	Matrix3d Q;

	//Stationary Correlation Matrix
	priorVariance(C);
	m[0]=p2*C.col(0)*(y[0]-mu)/(s2+p2*C(0,0));
	M[0]=p2*C-p2*C.col(0)*C.row(0)*p2/(s2+p2*C(0,0));
	//Forward Filtering
	for(int i=1;i<n;++i){
		Innovation3d(Q,t[i]-t[i-1],s2u,s2a, approx);
		Transition3d(A[i-1],t[i]-t[i-1],approx);
		AMAQ[i-1]=A[i-1]*M[i-1]*A[i-1].transpose()+p2*Q;
		m[i]=A[i-1]*m[i-1]+AMAQ[i-1].col(0)*(y[i]-mu-A[i-1].row(0)*m[i-1])/(s2+AMAQ[i-1](0,0));
		M[i]=AMAQ[i-1]-(AMAQ[i-1].col(0)*AMAQ[i-1].row(0))/(s2+AMAQ[i-1](0,0));

	}

	//Backward Sampling
	Matrix3d S=0.5*(M[n-1]+M[n-1].transpose());
	Eigen::LDLT<Matrix3d> ldlt;
	ldlt.compute(S);
	Matrix3d L(ldlt.matrixL());
	Vector3d D(ldlt.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);

	x.col(n-1)=m[n-1]+ldlt.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*Z.col(n-1);
	for(int i=(n-2); i>=0;--i){
		Eigen::JacobiSVD<Matrix3d> mysvd(AMAQ[i], Eigen::ComputeFullU | Eigen::ComputeFullV);
		Vector3d singularValues=mysvd.singularValues();
		double pinvtol=1.e-7;
		for(int v=0; v<3; ++v){
			if(singularValues(v)<pinvtol){
				singularValues(v)=0.0;
			}else{
				singularValues(v)=1.0/singularValues(v);
			}
		}
		Matrix3d Qinv=mysvd.matrixV()*singularValues.asDiagonal()*mysvd.matrixU().transpose();
		x.col(i)=m[i]+M[i]*A[i].transpose()*Qinv*(x.col(i+1)-A[i]*m[i]);
		S=M[i]-M[i]*A[i].transpose()*Qinv*(A[i]*M[i]);
		S=0.5*(S+S.transpose());
		ldlt.compute(S);
		Matrix3d L(ldlt.matrixL());
		Vector3d D(ldlt.vectorD());
		for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
		x.col(i)+=ldlt.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*Z.col(i);
	}
	Map<MatrixXd>( xout, x.rows(), x.cols() ) =  x;
}

extern "C" void ThreePoint(double * jx, double * jz, double t, double tl, double * jvl, double tr, double * jvr, double s2u, double s2a){
	Map<Vector3d> x(jx);
	Map<Vector3d> z(jz);
	Map<Vector3d> vr(jvr);
	Map<Vector3d> vl(jvl);
	Eigen::LDLT<Matrix3d> Q;
	Matrix3d Qrl,Ql,Al,Ar;
	Innovation3d(Qrl,tr-tl,s2u,s2a,false);
	Innovation3d(Ql,t-tl,s2u,s2a,false);
	Transition3d(Al,t-tl,false);
	Transition3d(Ar,tr-t,false);
	Eigen::JacobiSVD<Matrix3d> mysvd(Qrl, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Vector3d singularValues=mysvd.singularValues();
	double pinvtol=1.e-7;
	for(int v=0; v<3; ++v){
		if(singularValues(v)<pinvtol){
			singularValues(v)=0.0;
		}else{
			singularValues(v)=1.0/singularValues(v);
		}
	}
	Matrix3d Qinv=mysvd.matrixV()*singularValues.asDiagonal()*mysvd.matrixU().transpose();
	x=Al*vl+Ql*Ar.transpose()*Qinv*(vr-Ar*Al*vl);
	Q.compute(Ql-Ql*Ar.transpose()*Qinv*(Ar*Ql));
	Matrix3d L(Q.matrixL());
	Vector3d D(Q.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	x+=Q.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*z;
	Map<Vector3d>( jx,3,1) =   x;
}

extern "C" void Backwards(double * jx, double * jz,double t, double tr, double * jvr, double s2u, double s2a){
	Map<Vector3d> x(jx);
	Map<Vector3d> z(jz);
	Map<Vector3d> vr(jvr);
	Matrix3d Qr,Ar;
	RevTransition3d(Ar,tr-t);
	RevInnovation3d(Qr,tr-t,s2u,s2a);
	x=Ar*vr;
	Eigen::LDLT<Matrix3d> S;
	S.compute(Qr);
	Matrix3d L(S.matrixL());
	Vector3d D(S.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	x+=S.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*z;
	Map<Vector3d>(jx,3,1)=x;
}
extern "C" void Forwards(double * jx, double * jz, double t, double tl, double * jvl, double s2u, double s2a){
	Map<Vector3d> x(jx);
	Map<Vector3d> z(jz);
	Map<Vector3d> vl(jvl);
	Matrix3d Al,Ql;
	Innovation3d(Ql,t-tl,s2u,s2a,false);
	Transition3d(Al,t-tl,false);
	x=Al*vl;
	Eigen::LDLT<Matrix3d> S;
	S.compute(Ql);
	Matrix3d L(S.matrixL());
	Vector3d D(S.vectorD());
	for(int j=0;j<3;++j) D(j)=std::max(D(j),0.0);
	x+=S.transpositionsP().transpose()*L*D.cwiseSqrt().asDiagonal()*z;
	Map<Vector3d>(jx,3,1)=x;
}


extern "C" double LogDensity3d(double * y, double * t, int gamma, int n, double mu, double s2, double s2u, double s2a, double p2)
{
	if(s2u<0) return -1.0/0.0;
	if(s2a<0) return -1.0/0.0;
	if(p2<0) return -1.0/0.0;

	if(gamma==1){
		std::vector< Vector3d> m(n);
		std::vector< Matrix3d> M(n);
		std::vector< Matrix3d> AMAQ(n);
		std::vector< Matrix3d> A(n);
		Matrix3d C;
		priorVariance(C);
		Matrix3d Q;
		//Stationary Correlation Matrix
		m[0]=p2*C.col(0)*(y[0]-mu)/(s2+p2*C(0,0));
		M[0]=p2*C-p2*C.col(0)*C.row(0)*p2/(s2+p2*C(0,0));
		//Forward Filtering
		for(int i=1;i<n;++i){
			Innovation3d(Q,t[i]-t[i-1],s2u,s2a,false);
			Transition3d(A[i-1],t[i]-t[i-1],false);
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

extern "C" void mu3d(double * y, double * t, int n, double s2, double s2u, double s2a, double p2, double * mean, double * prec)
{

	std::vector< Matrix3d> M(n);
	std::vector< Matrix3d> AMAQ(n);
	std::vector< Matrix3d> A(n);
	Matrix3d Cor;
	Matrix3d Q;

	//Stationary Correlation Matrix
	priorVariance(Cor);
	Vector3d C=p2*Cor.col(0)*(y[0])/(s2+p2*Cor(0,0));
	Vector3d D=-p2*Cor.col(0)/(s2+p2*Cor(0,0));
	M[0]=p2*Cor-p2*Cor.col(0)*Cor.row(0)*p2/(s2+p2*Cor(0,0));
	*prec=1/(s2+p2*Cor(0,0));
	*mean=y[0]/(s2+p2*Cor(0,0));
	for(int i=2;i<n;++i){
		Innovation3d(Q,t[i]-t[i-1],s2u,s2a,false);
		Transition3d(A[i-1],t[i]-t[i-1],false);
		AMAQ[i-1]=A[i-1]*M[i-1]*A[i-1].transpose()+p2*Q;
		M[i]=AMAQ[i-1]-(AMAQ[i-1].col(0)*AMAQ[i-1].row(0))/(s2+AMAQ[i-1](0,0));
		Vector3d F=AMAQ[i-1].col(0)/(s2+AMAQ[i-1](0,0));
		C=C-F*A[i-1].row(0)*C+F*y[i-1];
		D=D-F*A[i-1].row(0)*D-F;
		*prec+=(1+A[i-1].row(0)*D)*(1+A[i-1].row(0)*D)/(s2+AMAQ[i-1](0,0));
		*mean+=(y[i]-A[i-1].row(0)*C)*(1+A[i-1].row(0)*D)/(s2+AMAQ[i-1](0,0));
	}
}

