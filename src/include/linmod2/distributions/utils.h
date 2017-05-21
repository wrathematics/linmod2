#ifndef __LINMOD2_DISTRIBUTIONS_UTILS__
#define __LINMOD2_DISTRIBUTIONS_UTILS__


#include <math.h>
#include "constants.h"
#include "normal.h"


static inline double erfinv(const double x)
{
  return qnorm((1. + x)/2., 0, 1, true, false) / ROOT2;
}



static inline double erfinvc(const double x)
{
  return erfinv(1. - x);
}



static inline double probit(const double x)
{
  return -ROOT2 * erfinvc(2. * x);
}


#endif
