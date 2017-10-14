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

uvec SampleIndex(const int n, const int d)
{
  //samples d elements in 1:n without replacement.
  int i, j, nn = n;
  uvec result(d);
  uvec ind = linspace<uvec>(n - 1, 0, n);
  vec Unifvals = randu<vec>(d);
  for (i = 0; i < d; i++)
  {
    j = Unifvals(i)*nn;
    result(i) = ind(j);
    ind(j) = ind(--nn);
  }
  return result;
}

void rhoHuber(vec &univariateSample, const double tuningConstant = 2.1)
{
  // changes univariatesample into rho(univariatesample)
  univariateSample.for_each([tuningConstant](mat::elem_type &value){
    double a = pow((value / tuningConstant), 2);
    value = a > 1 ? 1 : a * pow(1.54, 2);});
}

double scale1StepM(vec univariateSample, const double precScale =1e-10, const double tuningConstant = 2.1) 
{
  // Computes the first step of an algorithm for
  // a scale M - estimator using the given rho function.
  // The scatter is computed relative to zero.
  
  const double sigma0 = 1.4826 * median(abs(univariateSample));
  
  univariateSample = univariateSample / sigma0;
  rhoHuber(univariateSample, tuningConstant);
  return sigma0 * sqrt(sum(univariateSample) * 2 / univariateSample.size());
}

Splitsampleresult FastSplitSample(vec univariateSample) 
{
  // Centers sample by median, and divides in 2 equal halves.
  // Assumes that NAs have already been removed.
  // This function has time complexity O(n) if median has this complexity
  
  const int h = (int)(univariateSample.size() / 2);
  Splitsampleresult result = { zeros<vec>(h),zeros<vec>(h),0 };
  result.med = median(univariateSample);
  univariateSample = univariateSample - result.med;
  
  int sampleasize = 0;
  int samplebsize = 0;
  univariateSample.for_each([&](mat::elem_type &value)
  {
    if (value > 0)
    {
      result.samplea(sampleasize) = value;
      sampleasize = sampleasize + 1;
    }
    if (value < 0) 
    {
      result.sampleb(samplebsize) = value;
      samplebsize = samplebsize + 1;
    }
  });
  result.sampleb = abs(result.sampleb);
  
  return result;
}

SaSbresult compScales(vec univariateSample, const bool rmZeroes = false, const double maxRatio = 0, const double precScale =1e-10)
{
  // Computes the scales sa and sb(above and below the median).
  
  SaSbresult result = { 0,0,0 };
  Splitsampleresult splitresult = FastSplitSample(univariateSample);
  univariateSample.for_each([splitresult](mat::elem_type& val) { val = val - splitresult.med; });
  double scaleall = scale1StepM(univariateSample, precScale);
  if (rmZeroes) 
  {
    splitresult.samplea = splitresult.samplea(find(splitresult.samplea > precScale));
    splitresult.sampleb = splitresult.sampleb(find(splitresult.sampleb > precScale));
  }
  result.sa = scale1StepM(splitresult.samplea, precScale);
  result.sb = scale1StepM(splitresult.sampleb, precScale);
  result.med = splitresult.med;
  if (maxRatio >= 2)
  {
    result.sa = fmin(fmax(result.sa, (scaleall / maxRatio)), (scaleall * maxRatio));
    result.sb = fmin(fmax(result.sb, (scaleall / maxRatio)), (scaleall * maxRatio));
  }
  return result;
}

UnivDOresult DO_univ(vec univariateSamplex, vec univariateSamplez, const bool rmZeroes = false, const double maxRatio = 0, const double precScale = 1e-10)
{
  // Assumes that x1, x2 are arrays of numbers
  // Computes the directional outlyingness of each element in x2
  // with respect to x1
  // "rmZeroes" reduces the breakdown value but yields fewer implosions
  //"maxRatio" increases the breakdown value but reduces adaptation
  // to skewness.
  
  UnivDOresult result = { zeros<vec>(univariateSamplez.size()) ,zeros<vec>(univariateSamplez.size()) ,zeros<vec>(univariateSamplez.size()) };
  SaSbresult scales = compScales(univariateSamplex, rmZeroes, maxRatio, precScale);
  uvec inda = find(univariateSamplez > scales.med);
  uvec indb = find(univariateSamplez < scales.med);
  vec samplea = univariateSamplez(inda);
  vec sampleb = univariateSamplez(indb);
  result.numerator = abs(univariateSamplez - scales.med);
  result.denominator(inda).fill(scales.sa);
  result.diroutl(inda) = result.numerator(inda) / result.denominator(inda);
  result.denominator(indb).fill(scales.sb);
  result.diroutl(indb) = result.numerator(indb) / result.denominator(indb);
  return result;
}



generDirresult GenerDir(mat data, const int type,const int ndir, const double precScale = 1e-10) 
{
  // Generates `ndir' directions(with unit norm),
  // orthogonal to d - subsets of X
  // Assumes dim(data) = n x d
  
  mat A(ndir, data.n_cols, fill::zeros);
  int nbdirections = 0;
  int ssubsets = 0;
  
  if(type == 1){ // Affine invatiant
    vec B = ones<vec>(data.n_cols);
  while (nbdirections < ndir)
  {
    uvec indices = SampleIndex(data.n_rows, data.n_cols);
    mat Atemp = data.rows(indices);
    if (rank(Atemp, precScale) == data.n_cols) 
    {
      A.row(nbdirections) = solve(Atemp, B, solve_opts::fast).t();
      nbdirections = nbdirections + 1;
    }else{
      ssubsets = ssubsets + 1;
    }
    
  }
  }
  if(type == 2){ // Rotation invariant
    while (nbdirections < ndir)
    {
      uvec indices = SampleIndex(data.n_rows, 2);
      rowvec newdir = data.row(indices(1)) - data.row(indices(0));
      if(norm(newdir, 2) > precScale)
      {
        A.row(nbdirections) = newdir;
        nbdirections = nbdirections + 1;
      }else{
        ssubsets = ssubsets + 1;
      }
    }
  }
  if(type == 3){ // Shift invatiant
    while (nbdirections < ndir)
    {
      vec newdir(data.n_cols);
      newdir.randn();
      if(norm(newdir,2)> precScale)
      {
        A.row(nbdirections) = newdir.t();
        nbdirections = nbdirections + 1;
      }else{
        ssubsets = ssubsets + 1;
      }
    }
  }
  A = normalise(A, 2, 1); // normed directions
  generDirresult result = {A,ssubsets};
  
  return result;
}
}