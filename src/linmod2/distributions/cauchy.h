#ifndef __LINMOD2_DISTRIBUTIONS_CAUCHY__
#define __LINMOD2_DISTRIBUTIONS_CAUCHY__


#include <stdbool.h>
#include <math.h>

#include "constants.h"


// pdf
static inline double dcauchy(const double x, const double location, const double scale, const bool log_)
{
  double ret, tmp;
  
  tmp = (x - location)/scale;
  
  ret = 1.0 / (LINMOD2_PI*scale * (1.0 + tmp*tmp));
  
  if (log_ == true)
    ret = log(ret);
  
  return ret;
}



// cdf
static inline double pcauchy(const double q, const double location, const double scale, const bool lower_tail, const bool log_p)
{
  double ret;
  
  ret = LINMOD2_PIINV * atan((q-location)/scale) + 0.5;
  
  if (log_p == true)
    ret = log(ret);
  
  return ret;
}



// quantile
static inline double qcauchy(const double p, const double location, const double scale, const bool lower_tail, const bool log_p)
{
  double ret;
  
  ret = location + scale * tan(LINMOD2_PI*(p-0.5));
  
  if (log_p == true)
    ret = log(ret);
  
  return ret;
}


#endif
