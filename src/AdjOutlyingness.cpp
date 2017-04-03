#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <limits>
#include <inttypes.h>

#include <Eigen/Dense>
#include <Eigen/QR>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

template <typename Iterator>
extern inline bool next_combination(const Iterator first,Iterator k,const Iterator last);
extern double GetUniform(unsigned int *seed);
extern double GetNormal(unsigned int *seed);
extern VectorXi SampleD(const int& m,const int& p,VectorXi& y,unsigned int *seed);
extern VectorXi SampleR(const int& m,const int& p,VectorXi& ind,unsigned int *seed);
extern double quantiles(const Ref<const VectorXd>& x,const double quant);
extern void xad(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void xrd(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void aed(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void red(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void rsd(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void pCalc(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coeff,VectorXi& QIndexnin,unsigned int *seed,void (*pCalcMethod)(const MatrixXd&,const int&,const int&,int&,const double&,VectorXd&,VectorXi&,unsigned int *));
extern double sign(double x);
extern double h_kern(double a,double b,int ai,int bi,int ab,double eps);
extern double whimed_i(double a[],int w[],int n,double a_cand[],double a_srt[],int w_cand[]);
extern double mlmccN2(double x[],const int n);
extern double mlmccN(const double z[],const int n,const int dr);

int CAPD(VectorXd& m_resd,const double cupper,const double clower,const double coef,const int n1,const double tol){
	const int64_t n=m_resd.size();
	double tupD=0,tloD=0,param=0,tup=0,tlo=0,q1=0,q3=0,m_medn=0,tmc=0.0,IQR;
	const double Q1=0.25,Q2=0.5,Q3=0.75;
	int EF=0;
	VectorXd Y_resd=VectorXd::Zero(n);
	VectorXd W_resd=VectorXd::Zero(n);

	m_medn=quantiles(m_resd.head(n1),Q2);
	m_resd.array()-=m_medn;
	q1=quantiles(m_resd.head(n1),Q1);
	q3=quantiles(m_resd.head(n1),Q3);
	IQR=std::abs(q3-q1);
	if(IQR>tol){
		tmc=mlmccN(m_resd.head(n1).data(),n1,0);
		param=(tmc>=0)?(cupper*tmc):(tmc*clower);
		tup=q3+coef*IQR*exp(param);
		param=(tmc>=0)?(-1*clower*tmc):(-1*tmc*cupper);
		tlo=q1-coef*IQR*exp(param);

		tupD=m_resd.head(n1).minCoeff();
		Y_resd=(m_resd.array()>tup).select(tupD,m_resd);
		tupD=Y_resd.head(n1).maxCoeff();
		tloD=m_resd.head(n1).maxCoeff();
		Y_resd=(m_resd.array()<tlo).select(tloD,m_resd);
		tloD=-1.0*Y_resd.head(n1).minCoeff();
		W_resd.setConstant(tloD);
		Y_resd=(m_resd.array()>=0).select(tupD,W_resd);
		if(Y_resd.minCoeff()<tol){
		  VectorXd w_resd=VectorXd::Zero(n);
		  VectorXd d_resd=m_resd.array().abs();
		  m_resd=(d_resd.array()>tol).select(1.0,w_resd);
		  EF=1;
		} else {
		  m_resd=m_resd.array().abs();
		  m_resd.array()/=Y_resd.array();
		}
	} else {
		VectorXd w_resd=VectorXd::Zero(n);
		VectorXd d_resd=m_resd.array().abs();
		m_resd=(d_resd.array()>tol).select(1.0,w_resd);
		EF=1;
	}
	return(EF);
}
void Mainadjprojout(const MatrixXd& x,int& ndir,double& mcaap,int& sdr,int& type,const int& n1,VectorXd& Aoutlyingness,const int& exs,unsigned int *seed,  VectorXd&  projvector){
	void 	(*pFoo[])(const MatrixXd&,const int&,const int&,int&,const double&,VectorXd&,VectorXi&,unsigned int *)={&aed,&red,&rsd,&xad,&xrd};
	const int64_t p=x.cols(),n=x.rows();
	int m_r=0,j=0,EF=0;
	const double tol=std::numeric_limits<float>::min(),cupper=3.0,clower=4.0,coef=1.5;
	double tmc=0;
	VectorXi QIndexnin(n1);
	type+=(type<2 && exs)?3:0;
	VectorXd m_coef=VectorXd::Ones(p);
	VectorXd m_resd=VectorXd::Zero(n);
	QIndexnin.setLinSpaced(n1,0,n1-1);
	if(p>1){
		while(j<ndir){
			j++;
			pCalc(x,p,n1,m_r,tol,m_coef,QIndexnin,seed,pFoo[type]);
      if(m_r==p){
        m_resd=x*m_coef;
				EF=CAPD(m_resd,cupper,clower,coef,n1,tol);
				if(EF==0){
					Aoutlyingness=Aoutlyingness.cwiseMax(m_resd);
				} else {
        j=ndir;
  			EF=0;
				Aoutlyingness=m_resd;
        projvector = m_coef;
				}
			} else {
				sdr++;
			}
		}
	} else {
	  m_resd=x*m_coef;
		EF=CAPD(m_resd,cupper,clower,coef,n1,tol);
		if(EF==0){
		  Aoutlyingness=Aoutlyingness.cwiseMax(m_resd);
		} else {
		  j=ndir;
		  EF=0;
		  Aoutlyingness=m_resd;
		  projvector = m_coef;
		}
	}
	tmc=mlmccN(Aoutlyingness.head(n1).data(),n1,0);
	mcaap=tmc;
}
extern "C"{
	void adjprojout(int* n,int* p,int* ndir,double* xi,double* Adjprojout,double* mcaap,int* sdr,int* type,int* n1,int* exs, double* projectionvector, unsigned int* seed){
		int isdr=0,itype=*type-1,indir=*ndir;
		const int in1=*n1,xs=*exs;
		double inmcaap=0.0;
    unsigned int* iseed = seed;
		MatrixXd x=Map<MatrixXd>(xi,*n,*p);
		VectorXd outlyingness=VectorXd::Zero(*n);
    VectorXd projvector=VectorXd::Zero(*p);
		Mainadjprojout(x,indir,inmcaap,isdr,itype,in1,outlyingness,xs,iseed, projvector);
 		Map<VectorXd>(Adjprojout,*n)=outlyingness.array();
    Map<VectorXd>(projectionvector,*p)=projvector.array();
		*mcaap=inmcaap;
		*sdr=isdr;
	}
}
