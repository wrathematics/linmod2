#ifndef __MMULT_H__
#define __MMULT_H__


#include <stdbool.h>
#include "lapack.h"


// dgemm wrapper
static inline void matmult(const bool transx, const bool transy, const double alpha, const int mx, const int nx, const double *const restrict x, const int my, const int ny, const double *const restrict y, double *restrict ret)
{
  // m = # rows of op(x)
  // n = # cols of op(y)
  // k = # cols of op(x)
  int im, in, ik;
  char ctransx, ctransy;
  static const double zero = 0.;
  
  ctransx = transx ? 'T' : 'N';
  ctransy = transy ? 'T' : 'N';
  
  if (transx)
  {
    im = nx;
    ik = mx;
  }
  else
  {
    im = mx;
    ik = nx;
  }
  
  in = transy ? my : ny;
  
  dgemm_(&ctransx, &ctransy, &im, &in, &ik, &alpha, x, &mx, y, &my, &zero, ret, &im);
}


#endif
