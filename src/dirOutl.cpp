#include "dirOutl.h"
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace dirout;

/******************************************************/
/*       Main Directional Outlyingness function       */
/******************************************************/

SEXP dirOutl_cpp(SEXP DATAX, SEXP DATAZ, SEXP TYPE,SEXP NDIR, SEXP RMZEROES, SEXP MAXRATIO,SEXP PRECSCALE)
{
  
  try
  {
    // Input arguments
    mat datax = as<mat>(DATAX);
    mat dataz = as<mat>(DATAZ);
    const int type = as<int>(TYPE);
    const unsigned int ndir = as<unsigned int>(NDIR);
    const bool rmZeroes = as<bool>(RMZEROES); 
    const double maxRatio = as<double>(MAXRATIO);
    const double precScale = as<double>(PRECSCALE);
    
    // Create output vectors
    vec outlyingnessX = zeros<vec>(datax.n_rows);
    vec outlyingnessZ = zeros<vec>(dataz.n_rows);
    int error_code = 0; // 0 = no error, 1 = direction with zero scale
    int singularsubsets = 0;
    vec hyperplane =  zeros<vec>(datax.n_cols);
    
    if(type != 0) // Projection Pursuit
    {
      // Generate directions
      generDirresult dirtemp = GenerDir(datax, type, ndir, precScale);
      mat directions = dirtemp.directions;
      singularsubsets = dirtemp.singularsubsets;
      // Project data on the generated directions
      mat Yx = datax * directions.t(); 
      mat Yz = dataz * directions.t(); 
      // Check for zero scales
      uword i = 0;
      do 
      {
        if ((median(abs(Yx.col(i) - median(Yx.col(i))))) < precScale)
        {
          error_code = 1;
          hyperplane = directions.row(i).t();
        }
        i++;
      }
      while ((i < Yx.n_cols) && (error_code==0));
      
      // Continue if no zero scales
      if(error_code == 0)
      {
        for (uword i = 0; i < Yx.n_cols; i++)
        {
          vec tempresult(Yx.n_rows + Yz.n_rows, fill::zeros);
          tempresult.subvec(0, (Yx.n_rows - 1)) = Yx.col(i);
          tempresult.subvec(Yx.n_rows, (Yx.n_rows + Yz.n_rows - 1)) = Yz.col(i);
          tempresult = DO_univ(Yx.col(i), tempresult, rmZeroes, maxRatio, precScale).diroutl;
          Yx.col(i) = tempresult.subvec(0, (Yx.n_rows - 1));
          Yz.col(i) = tempresult.subvec(Yx.n_rows, (Yx.n_rows + Yz.n_rows - 1));
        }
        for (uword i = 0; i < Yx.n_rows; i++) 
        {
          outlyingnessX(i) = max(Yx.row(i));
        }
        for (uword i = 0; i < Yz.n_rows; i++)
        {
          outlyingnessZ(i) = max(Yz.row(i));
        }
      }
    }
    if(type == 0) // componentwise
    {
      for (uword i = 0; i < datax.n_cols; i++)
      {
        vec tempresult(datax.n_rows + dataz.n_rows, fill::zeros);
        tempresult.subvec(0, (datax.n_rows - 1)) = datax.col(i);
        tempresult.subvec(datax.n_rows, (datax.n_rows + dataz.n_rows - 1)) = dataz.col(i);
        tempresult = DO_univ(datax.col(i), tempresult, rmZeroes, maxRatio, precScale).diroutl;
        outlyingnessX = outlyingnessX + arma::pow(tempresult.subvec(0, (datax.n_rows - 1)),2);
        outlyingnessZ = outlyingnessZ + arma::pow(tempresult.subvec(datax.n_rows, (datax.n_rows + dataz.n_rows - 1)),2);
      }
      outlyingnessX = arma::sqrt(outlyingnessX);
      outlyingnessZ = arma::sqrt(outlyingnessZ);
    }
    
    return List::create( Named("outlyingnessX") = outlyingnessX,
                         Named("outlyingnessZ") = outlyingnessZ, 
                         Named("error_code") = error_code,
                         Named("singularsubsets") = singularsubsets,
                         Named("hyperplane") = hyperplane);
    
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return wrap(NA_REAL);
}



