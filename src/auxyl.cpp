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

#include <Eigen/Dense>
#include <Eigen/QR>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

extern double GetUniform(unsigned int *seed);
extern double GetNormal(unsigned int *seed);

template <typename Iterator>
inline bool next_combination(const Iterator first, Iterator k, const Iterator last){
	/*
	Computes the next lexicographic permutation of 1:n.
	Credits: Thomas Draper Via http://stackoverflow.com/a/5097100/189035
	*/
	if ((first==last) || (first==k) || (last==k))	return false;
	Iterator itr1=first;
	Iterator itr2=last;
	++itr1;
	if (last==itr1)				return false;
	itr1=last;
	--itr1;
	itr1=k;
	--itr2;
	while(first!=itr1){
		if (*--itr1<*itr2){
			Iterator j=k;
			while(!(*itr1<*j)) ++j;
			std::iter_swap(itr1,j);
			++itr1;
			++j;
			itr2=k;
			std::rotate(itr1,j,last);
			while(last!=j){
				++j;
				++itr2;
			}
			std::rotate(k,itr2,last);
			return true;
		}
	}
	std::rotate(first,k,last);
	return false;
}
VectorXi SampleR(const int& m,const int& p,VectorXi& ind,unsigned int *seed){
//samples p elements in 1:n without replacement.
	int i,j,nn=m;
	VectorXi y(p);
	ind.setLinSpaced(m,0,m-1);
	for(i=0;i<p;i++){
		j=GetUniform(seed)*nn;
		y(i)=ind(j);
		ind(j)=ind(--nn);
	}
	return(y);
}
VectorXi SampleD(const int& p,VectorXi& y,unsigned int *seed){
//computes the next lexicographic combination of p out of 1:n elements.
	next_combination(y.data(),y.data()+p,y.data()+y.size());
	return(y.head(p));
}
double quantiles(const Ref<const VectorXd>& x,const double quant){
//computes the quantile 'quant' of x.
	const int n=x.size();
	double lq,uq,fq;
	const double q1=n*(double)quant+0.5;
	const int index1=floor(q1);
	const int index2=ceil(q1);
	const double index3=(double)index2-q1;
	VectorXd x1=x;
	std::nth_element(x1.data(),x1.data()+index1-1,x1.data()+x1.size());
	lq=x1(index1-1);
	if(index1==index2){
		fq=lq;
	} else {
		uq=x1.segment(index1,x1.size()-index1-1).minCoeff();
		fq=lq*index3+uq*(1.0-index3);
	}
	return(fq);
}
void xad(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed){
//QR-computes normed the slope of the hyperplane through p points (this is the exhaustive enumeration version).
	MatrixXd A(p,p);
	VectorXd y=VectorXd::Ones(p);
	VectorXi QIndexpin=QIndexnin.head(p);
	for(int i=0;i<p;i++)	A.row(i)=x.row(QIndexpin(i));
	ColPivHouseholderQR<MatrixXd>PQR(A);
	m_r=PQR.rank();
  if(m_r==p)		m_coef=PQR.solve(y);
	SampleD(p,QIndexnin,seed);
}
void xrd(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed){
//computes the normed slope of the hyperplane through 2 points (this is the exhaustive emuneration version).
	VectorXi QIndexpin=QIndexnin.head(2);
	m_coef=x.row(QIndexpin(0))-x.row(QIndexpin(1));
	double m_medn=m_coef.norm();
	if(m_medn>tol){
		m_coef.array()/=m_medn;
		m_r=p;
	}
  else{
    m_r=0;
  }
	SampleD(2,QIndexnin,seed);
}
void aed(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed){
//computes the normed slope of the hyperplane through p points sampled randomly.
	MatrixXd A(p,p);
	VectorXd y=VectorXd::Ones(p);
	VectorXi QIndexpin(p);
	QIndexpin=SampleR(n1,p,QIndexnin,seed);
	for(int i=0;i<p;i++)	A.row(i)=x.row(QIndexpin(i));
	ColPivHouseholderQR<MatrixXd>PQR(A);
	m_r=PQR.rank();
  if(m_r==p)		m_coef=PQR.solve(y);
}
void red(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed){
//computes the normed slope of the hyperplane through 2 points sampled randomly.
	VectorXi QIndexpin(2);
	QIndexnin.setLinSpaced(n1,0,n1-1);
	QIndexpin=SampleR(n1,2,QIndexnin,seed);
	m_coef=x.row(QIndexpin(0))-x.row(QIndexpin(1));
	double m_medn=m_coef.norm();
	if(m_medn>tol){
		m_coef.array()/=m_medn;
		m_r=p;
	}
  else{
    m_r=0;
  }
}
void rsd(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexin,unsigned int *seed){
//p points on the unit sphere.
	for(int i=0;i<p;i++)	m_coef(i)=GetNormal(seed);
	double m_medn=m_coef.norm();
  if(m_medn>tol){
		m_coef.array()/=m_medn;
		m_r=p;
	}
  else{
    m_r=0;
  }
}
void pCalc(const MatrixXd& x,const int& p,const int& n1,int& m_r,const double& tol,VectorXd& m_coef,VectorXi& QIndexnin,unsigned int *seed,void (*pCalcMethod)(const MatrixXd&,const int&,const int&,int&,const double&,VectorXd&,VectorXi&,unsigned int *)){
//pointer to the 5 directions generating functions above.
        pCalcMethod(x,p,n1,m_r,tol,m_coef,QIndexnin,seed);
}
double MedMad(VectorXd& m_resd,const int& n1,const int& h,const double& factor){
//robust Z-scores of a vector using the median and the mad.
	const double Q2=0.5;
	m_resd.array()-=quantiles(m_resd.head(n1),Q2);
	m_resd=m_resd.array().abs();
	return(1.4826*quantiles(m_resd.head(n1),Q2));
}
double cMed(VectorXd& m_resd,const int& n1,const int& h,const double& factor){
//robust Z-scores of a centered vector using the median of absolute values.
	const double Q2=0.5;
	m_resd=m_resd.array().abs();
	return(quantiles(m_resd.head(n1),Q2));
}
double cMcd(VectorXd& m_resd,const int& n1,const int& h,const double& factor){
//robust Z-scores of a centered vector using the sd of the h observations closest to 0.
	double initcov=0.0;
	if(h==n1){
		initcov=(m_resd.head(n1).array()).abs2().sum()/(double)(h-1);
		return(sqrt(initcov));
	}
	VectorXd y3(n1);
	VectorXd y2(n1);
	y3=m_resd.head(n1).array().abs();
	std::nth_element(y3.data(),y3.data()+h,y3.data()+y3.size());
	initcov=y3.head(h).array().square().sum()/(double)(h);
//RW
	y3=(m_resd.head(n1).array()).array().abs2()/initcov;
	std::nth_element(y3.data(),y3.data()+h-1,y3.data()+y2.size());
	initcov*=y3(h-1)/factor;
	y3=(m_resd.head(n1).array()).array().abs2()/initcov;
	y2.setOnes();
	y2=(y3.array()>5.023886).select(0.0,y2);
	initcov=(y2.array()*y3.array().abs2()).sum()/(y2.array().sum()-1.0);
	initcov=sqrt(initcov);
	m_resd=m_resd.array().abs();
	return(initcov);
}
double unimcd(VectorXd& m_resd,const int& n1,const int& h,const double& factor){
//robust Z-scores of a vector using the reweighted unimcd estimates of location and scatter.
	const int len=n1-h+1;
	double initmean=0.0,initcov=0.0,sumw=0.0;
	int minone;
	if(h==n1){
		initmean=m_resd.head(n1).sum()/(double)h;
		initcov=(m_resd.head(n1).array()-initmean).abs2().sum()/(double)(h-1);
		return(sqrt(initcov));
	}
	VectorXd y=m_resd.head(n1);
	VectorXd ay(len);
	VectorXd ay2(len);
	VectorXd sq(len);
	VectorXd y2(n1);

	std::sort(y.data(),y.data()+y.size());
	ay(0)=y.head(h).sum();
	for(int samp=1;samp<len;samp++) ay(samp)=ay(samp-1)-y(samp-1)+y(samp+h-1);
	ay2=ay.array().square()/(double)h;
	y2=y.array().square();
	sq(0)=y2.head(h).sum()-ay2(0);
	for(int samp=1;samp<len;samp++) sq(samp)=sq(samp-1)-y2(samp-1)+y2(samp+h-1)-ay2(samp)+ay2(samp-1);
	initcov=sq.minCoeff(&minone);
	initcov/=(double)(h-1);
	initmean=ay(minone)/(double)h;
//RW
	y2=(m_resd.head(n1).array()-initmean).array().abs2()/initcov;
	std::nth_element(y2.data(),y2.data()+h-1,y2.data()+y2.size());
	initcov*=y2(h-1)/factor;
	y=(m_resd.head(n1).array()-initmean).array().abs2()/initcov;
	y2.setOnes();
	y2=(y.array()>5.023886).select(0.0,y2);
	sumw=y2.array().sum();
	initmean=(y2.array()*m_resd.head(n1).array()).sum()/sumw;
	m_resd.array()-=initmean;
	initcov=(y2.array()*m_resd.head(n1).array().abs2()).sum()/(sumw-1.0);
	initcov=sqrt(initcov);
	m_resd=m_resd.array().abs();
	return(initcov);
}
double pScal(VectorXd& m_resd,const int& n1,const int& h,const double& factor,double (*qCalcMethod)(VectorXd&,const int&,const int&,const double&)){
//pointer to the Z score methods above.
        return(qCalcMethod(m_resd,n1,h,factor));
}
