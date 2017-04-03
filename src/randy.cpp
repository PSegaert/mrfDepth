/* Reentrant random function from POSIX.1c.
   Copyright (C) 1996, 1999, 2009 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@cygnus.com>, 1996.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */
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
double GetNormal(unsigned int *seed){	/* normal random variate generator, mean m, standard deviation s */
/* boxmuller.c           Implements the Polar form of the Box-Muller
                         Transformation
                      (c) Copyright 1994, Everett F. Carter Jr.
                          Permission is granted by the author to use
			  this software for any application provided this
			  copyright notice is preserved.
*/
	
	double x1,x2,w,y1,m=0.0,s=1.0;
	static double y2;
	static int use_last = 0;
	if(use_last){		        /* use value from previous call */
		y1=y2;
		use_last=0;
	} else {
		do{
			x1=2.0*GetUniform(seed)-1.0;
			x2=2.0*GetUniform(seed)-1.0;
			w=x1*x1+x2*x2;
		} while (w>=1.0);
		w=sqrt((-2.0*log(w))/w);
		y1=x1*w;
		y2=x2*w;
		use_last=1;
	}
	return(m+y1*s);
}
