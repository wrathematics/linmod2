#ifndef __LINMOD2_GLM_FAMILY_UTILS__
#define __LINMOD2_GLM_FAMILY_UTILS__


#include "glm_constants.h"

static inline int glm_check_fitted(const int family, const int n, double *const restrict mu)
{
  int i;
  
  if (family == GLM_FAMILY_BINOMIAL)
  {
    for (i=0; i<n; i++)
    {
      if (mu[i] >= 1.0 || mu[i] <= 0.0)
        return GLM_FAMILY_BADMU;
    }
  }
  else if (family == GLM_FAMILY_POISSON || family == GLM_FAMILY_GAMMA)
  {
    for (i=0; i<n; i++)
    {
      if (mu[i] <= 0.0)
        return GLM_FAMILY_BADMU;
  }
  
  return 0;
}



// initialize mu based on the error distribution
static inline void glm_initial_mu(const int family, const int n, const double *const restrict y, const double *const restrict wt, double *const restrict mu)
{
  int i;
  
  if (family == GLM_FAMILY_BINOMIAL)
  {
    for (i=0; i<n; i++)
      mu[i] = (wt[i] * y[i] + 0.5) / (wt[i] + 1.0);
  }
  else if (family == GLM_FAMILY_GAMMA || family == GLM_FAMILY_GAUSSIAN || family == GLM_FAMILY_INVERSEGAUSSIAN)
  {
    for (i=0; i<n; i++)
      mu[i] = y[i];
  }
  else if (family == GLM_FAMILY_POISSON)
  {
    for (i=0; i<n; i++)
      mu[i] = y[i] + 0.1;
  }
}



static inline void glm_variance(const int family, const int n, const double *const restrict mu, double *const restrict var)
{
  int i;
  
  if (family == GLM_FAMILY_BINOMIAL)
  {
    for (i=0; i<n; i++)
      var[i] = mu[i] * (1.0d0 - mu[i]);
  }
  else if (family == GLM_FAMILY_GAMMA)
  {
    for (i=0; i<n; i++)
      var[i] = mu[i] * mu[i];
  }
  else if (family == GLM_FAMILY_GAUSSIAN)
  {
    for (i=0; i<n; i++)
      var[i] = 1.0;
  }
  else if (family == GLM_FAMILY_POISSON)
  {
    for (i=0; i<n; i++)
      var[i] = mu[i];
  }
  else if (family == GLM_FAMILY_INVERSEGAUSSIAN)
  {
    for (i=0; i<n; i++)
      var[i] = mu[i] * mu[i] * mu[i];
  }
}


#endif
