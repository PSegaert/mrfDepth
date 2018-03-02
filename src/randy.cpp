#include <stdlib.h>
#include <iostream>
#include <limits>
#include <cmath>
/* This algorithm is mentioned in the ISO C standard, here extended
   for 32 bits.  */
int randy_r(unsigned int *seed){
  unsigned int next=*seed;
  int result;

  next*=1103515245;
  next+=12345;
  result=(unsigned int) (next/65536)%2048;

  next*=1103515245;
  next+=12345;
  result<<=10;
  result^=(unsigned int) (next/65536)%1024;

  next*=1103515245;
  next+=12345;
  result<<=10;
  result^=(unsigned int) (next/65536)%1024;

  *seed=next;

  return result;
}
double GetUniform(unsigned int *seed){
//generates a U[0,1( random variable.
    return randy_r(seed)/(double)std::numeric_limits<int>::max();
}

double GetNormal(unsigned int *seed){

  double mean = 0.0, stddev = 1.0;
  static double n2;
  static int n2_cached = 0;
  if(n2_cached){		        
    n2_cached = 0;
    return n2*stddev + mean;
  } else {
    
    double x, y, r;
    do
    {
      x = 2.0 * GetUniform(seed) - 1;
      y = 2.0 * GetUniform(seed) - 1;
      r = x*x + y*y;
    }
    while (r == 0.0 || r > 1.0);
    {
      double d = sqrt(-2.0*log(r)/r);
      double n1 = x*d;
      n2 = y*d;
      double result = n1*stddev + mean;
      n2_cached = 1;
      return result;
    }
  }
}
