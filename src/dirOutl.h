#ifndef dirOutl_H
#define dirOutl_H


#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
#include "omp.h"

RcppExport SEXP dirOutl_cpp(SEXP DATAX, SEXP DATAZ, SEXP TYPE,SEXP NDIR, SEXP RMZEROES, SEXP MAXRATIO,SEXP PRECSCALE);


#ifndef  ARMA_USE_CXX11
#define ARMA_USE_CXX11
#endif 

using namespace arma;

namespace dirout
{

struct Splitsampleresult
{
  vec samplea;
  vec sampleb;
  double med;
};

struct UnivDOresult 
{
  vec numerator;
  vec denominator;
  vec diroutl;
};

struct SaSbresult 
{
  double sa;
  double sb;
  double med;
};

struct generDirresult
{
  mat directions;
  int singularsubsets;
};

uvec SampleIndex(const int n, const int d);

void rhoHuber(vec &univariateSample, const double tuningConstant = 2.1);

double scale1StepM(vec univariateSample,
                   const double precScale =1e-10,
                   const double tuningConstant = 2.1);

Splitsampleresult FastSplitSample(vec univariateSample);

SaSbresult compScales(vec univariateSample,
                      const bool rmZeroes = false,
                      const double maxRatio = 0,
                      const double precScale =1e-10);

UnivDOresult DO_univ(vec univariateSamplex,
                     vec univariateSamplez,
                     const bool rmZeroes = false,
                     const double maxRatio = 0,
                     const double precScale = 1e-10);


generDirresult GenerDir(mat data, const int type,
                        const int ndir,
                        const double precScale = 1e-10);
}

#endif