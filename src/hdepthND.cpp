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
#include <time.h>

#include <Eigen/Dense>
#include <Eigen/QR>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

struct IdLess {					//internal function.
    template <typename T>
    IdLess(T iter) : values(&*iter) {}
    bool operator()(int left,int right){
		const double tol = std::numeric_limits<float>::min();
		// We make sure that the points from X come before the points from Z:
		// if points are equal we sort based on position in original vector.
		// We use the fact that by construction in the original vector contains
		// first the data points from the data group and secondly fron the Z group.
		if (abs(values[left]-values[right])<tol){
			return left < right;
		}
		else{
			return values[left]<values[right];
		}
    }
    double const* values;
};
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

void hddepth(const VectorXd& x,const int& n1,VectorXi& SIndx1){
	const int n=x.size();
	const double tol=std::numeric_limits<float>::min();
	const int nX = n1;
	int smallerCount;
	int biggerCount;
	VectorXi SIndx2(n);
	VectorXi BiggerDepth(n);

	SIndx2.setLinSpaced(n,0,n-1);
	std::sort(SIndx2.data(),SIndx2.data()+SIndx2.size(),IdLess(x.data()));
	smallerCount = 0;
	biggerCount = 0;
	// loop trough the vector and count the x encountered up to and the current position.
	for (int i = 0; i < n; i++){
		if ((SIndx2(i) + 1) <= nX){ smallerCount = smallerCount + 1; }
		SIndx1(SIndx2(i)) = smallerCount;
		//Watch out for zero indexing
		if ((SIndx2( (n-1) - i ) + 1) <= nX){ biggerCount = biggerCount + 1; }
		BiggerDepth(SIndx2((n - 1) - i)) = biggerCount;
	}

	//correct for ties: since x's are now before z's.
	//we can simple loop from the back and when there are ties set smallerCount the latest smallerCount.
	//for the biggerCount we loop from the front and when there are ties we set biggerCount to the first biggerCount.
	for (int i = (n - 1); i > 0; i--){
		if (abs( x(SIndx2(i)) - x(SIndx2(i - 1)) ) < tol ) { SIndx1(SIndx2(i - 1)) = SIndx1(SIndx2(i)); }
		//Check Indices !
		if (abs(x(SIndx2((n - 1) - i + 1)) - x(SIndx2((n - 1) - i))) <tol){ BiggerDepth(SIndx2( (n - 1) - i + 1 )) = BiggerDepth(SIndx2( (n - 1) - i )); }

	}

	SIndx1 = SIndx1.cwiseMin(BiggerDepth);

}
void MainprojHSD(const MatrixXd& x,int& ndir,int& sdr,int& type,const int& n1,VectorXi& outlyingness,const int& exs,unsigned int *seed){
//carries the actual computation of the projection outlyingness.
	void 	(*pFoo[])(const MatrixXd&,const int&,const int&,int&,const double&,VectorXd&,VectorXi&,unsigned int *)={&aed,&red,&rsd,&xad,&xrd};
	const int p=x.cols(),n=x.rows();
	int m_r=0,j=0;
	const double tol=std::numeric_limits<float>::min();
	VectorXi QIndexnin(n1);
	VectorXi SIndx1(n);
	VectorXd m_coef=VectorXd::Ones(p);
	VectorXd m_resd=VectorXd::Zero(n);

	type+=(type<2 && exs)?3:0;
	QIndexnin.setLinSpaced(n1,0,n1-1);
	if(p>1){
	  while(j<ndir){
  		j++;
  		pCalc(x,p,n1,m_r,tol,m_coef,QIndexnin,seed,pFoo[type]);
  		if(m_r==p){
      	m_resd=x*m_coef;
  			hddepth(m_resd,n1,SIndx1);
  			outlyingness=outlyingness.cwiseMin(SIndx1);
  		} else {
  			sdr++;
  		}
  	}
	} else {
	  m_resd=x*m_coef;
	  hddepth(m_resd,n1,SIndx1);
	  outlyingness=outlyingness.cwiseMin(SIndx1);
	}
}

extern "C"{
	void HSDND(int* n,int* p,int* ndir,double* xi,int* projoutlyingness,int* sdr,int* type,int* n1,int* ex, unsigned int* seed){
//C wrapper for link to R.
		const int in1=*n1,exs=*ex,in=*n;
		int itype=*type-1,indir=*ndir,isdr=0;
		unsigned int* iseed = seed;
		MatrixXd x=Map<MatrixXd>(xi,in,*p);
		VectorXi HSdepth(in);
		HSdepth.setLinSpaced(in,in,in);
		MainprojHSD(x,indir,isdr,itype,in1,HSdepth,exs,iseed);
		Map<VectorXi>(projoutlyingness,in)=HSdepth.array();
		*sdr=isdr;
 	}
}

void MainprojHSDFAST(const MatrixXd& x,int& ndir,int& sdr,int& type,const int& n1,VectorXi& outlyingness,const int& exs,unsigned int *seed, MatrixXd& Directions, int& DirFlag){
  //carries the actual computation of the projection outlyingness.
  void 	(*pFoo[])(const MatrixXd&,const int&,const int&,int&,const double&,VectorXd&,VectorXi&,unsigned int *)={&aed,&red,&rsd,&xad,&xrd};
  const int p=x.cols(),n=x.rows();
  int m_r=0,j=0;
  const double tol=std::numeric_limits<float>::min();
  VectorXi QIndexnin(n1);
  VectorXi SIndx1(n);
  VectorXd m_coef(p);
  VectorXd m_resd(n);

  type+=(type<2 && exs)?3:0;
  QIndexnin.setLinSpaced(n1,0,n1-1);
  while(j<ndir){
    j++;
    if(DirFlag==0) {
      pCalc(x,p,n1,m_r,tol,m_coef,QIndexnin,seed,pFoo[type]);
      Directions.row(j-1) = m_coef;
    } else {
      m_coef = Directions.row(j-1);
      m_r = p;
    }
    if(m_r==p){
      m_resd=x*m_coef;
      hddepth(m_resd,n1,SIndx1);
      outlyingness=outlyingness.cwiseMin(SIndx1);
    } else {
      sdr++;
    }
  }
}

extern "C"{
  void HSDNDFast(int* n,int* p,int* ndir,double* xi,int* projoutlyingness,int* sdr,int* type,int* n1,int* ex, unsigned int* seed, double* Directions, int* DirFlag){
    //C wrapper for link to R.
    const int in1=*n1,exs=*ex,in=*n;
    int itype=*type-1,indir=*ndir,isdr=0;
    unsigned int* iseed = seed;
    MatrixXd x = Map<MatrixXd>(xi,in,*p);
    MatrixXd Dirs = Map<MatrixXd>(Directions,indir,*p);
    VectorXi HSdepth(in);
    HSdepth.setLinSpaced(in,in,in);
    MainprojHSDFAST(x,indir,isdr,itype,in1,HSdepth,exs,iseed, Dirs, *DirFlag);
    Map<VectorXi>(projoutlyingness,in)=HSdepth.array();
    Map<MatrixXd>(Directions,indir,*p)=Dirs.array();
    *sdr=isdr;
  }
}


