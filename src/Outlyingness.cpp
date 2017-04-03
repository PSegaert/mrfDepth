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

#include <Eigen/Dense>
#include <Eigen/QR>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

template <typename Iterator>
extern inline bool next_combination(const Iterator first, Iterator k, const Iterator last);
extern double GetUniform(unsigned int *seed);
extern double GetNormal(unsigned int *seed);
extern VectorXi SampleR(const int& m,const int& p,VectorXi& ind,unsigned int *seed);
extern VectorXi SampleD(const int& p,VectorXi& y,unsigned int *seed);
extern double quantiles(const Ref<const VectorXd>& x,const double quant);
extern void xad(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void xrd(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void aed(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void red(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed);
extern void rsd(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexin,unsigned int *seed);
extern void pCalc(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed,void (*pCalcMethod)(const MatrixXd&,const int&,const int&,int&,const double&,VectorXd&,VectorXi&,unsigned int *));
extern double MedMad(VectorXd& m_resd,const int& n1,const int& h,const double& factor);
extern double cMed(VectorXd& m_resd,const int& n1,const int& h,const double& factor);
extern double cMcd(VectorXd& m_resd,const int& n1,const int& h,const double& factor);
extern double unimcd(VectorXd& m_resd,const int& n1,const int& h,const double& factor);
extern double pScal(VectorXd& m_resd,const int& n1,const int& h,const double& factor,double (*qCalcMethod)(VectorXd&,const int&,const int&,const double&));

void Mainprojoutlyingness(const MatrixXd& x,int& ndir,int& sdr,int& type,const int& n1,VectorXd& outlyingness,int& scaleT,const int& h,const int& center,const double& factor,const int& exs,unsigned int *seed, VectorXd& projvector){
//carries the actual computation of the projection outlyingness.
	void 	(*pFoo[])(const MatrixXd&,const int&,const int&,int&,const double&,VectorXd&,VectorXi&,unsigned int *)={&aed,&red,&rsd,&xad,&xrd};
	double 	(*qFoo[])(VectorXd&,const int&,const int&,const double&)={&MedMad,&unimcd,&cMed,&cMcd};
	const int p=x.cols(),n=x.rows();
	int m_r=0,j=0;
	const double tol=std::numeric_limits<float>::min();
	double m_medn=.0;
	VectorXi QIndexnin(n1);
	type+=(type<2 && exs)?3:0;
	scaleT+=(center==1)?0:2;
	VectorXd m_coef=VectorXd::Ones(p);
	VectorXd m_resd=VectorXd::Zero(n);
	QIndexnin.setLinSpaced(n1,0,n1-1);
	if(p>1){
	  while(j<ndir){
		j++;
		pCalc(x,p,n1,m_r,tol,m_coef,QIndexnin,seed,pFoo[type]);
		if(m_r==p){
	    m_resd=x*m_coef;
			m_medn=pScal(m_resd,n1,h,factor,qFoo[scaleT]);
			if(m_medn>tol){
				m_resd.array()/=m_medn;
				outlyingness=outlyingness.cwiseMax(m_resd);
			} else {
				j=ndir;
  			VectorXd w_resd=VectorXd::Zero(n);
				outlyingness=(m_resd.array()>tol).select(1.0,w_resd);
        projvector = m_coef;
			}
		} else {
			sdr++;
		}
	  }
	} else {
	  m_resd=x*m_coef;
	  m_medn=pScal(m_resd,n1,h,factor,qFoo[scaleT]);
	  if(m_medn>tol){
	    m_resd.array()/=m_medn;
	    outlyingness=outlyingness.cwiseMax(m_resd);
	  } else {
	    j=ndir;
	    VectorXd w_resd=VectorXd::Zero(n);
	    outlyingness=(m_resd.array()>tol).select(1.0,w_resd);
	    projvector = m_coef;
	  }
	}
}
extern "C"{
	void projoutlyingness(int* n,int* p,int* ndir,double* xi,double* projoutlyingness,int* sdr,int* type,int* n1,int* itypes,int* ih,int* icenter,double* rfactor,int* ex, double* projectionvector, unsigned int* seed){
//C wrapper for link to R.
		const int in1=*n1,h=*ih,center=*icenter,exs=*ex;
		const double factor=*rfactor;
		int itype=*type-1,indir=*ndir,types=*itypes-1,isdr=0;
		unsigned int* iseed = seed;
		MatrixXd x=Map<MatrixXd>(xi,*n,*p);
		VectorXd outlyingness=VectorXd::Zero(*n);
    VectorXd projvector=VectorXd::Zero(*p);
		Mainprojoutlyingness(x,indir,isdr,itype,in1,outlyingness,types,h,center,factor,exs,iseed, projvector);
 		Map<VectorXd>(projoutlyingness,*n)=outlyingness.array();
    Map<VectorXd>(projectionvector,*p)=projvector.array();
    *sdr=isdr;
	}
}
