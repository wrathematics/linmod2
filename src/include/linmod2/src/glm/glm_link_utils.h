#ifndef __LINMOD2_GLM_LINK_UTILS__
#define __LINMOD2_GLM_LINK_UTILS__


#include <math.h>
#include "glm_constants.h"

static inline void glm_link(const int link, const int n, const double *const restrict x, double *const restrict y)
{
  int i;
  
  if (link == GLM_LINK_CLOGLOG)
  {
    for (i=0; i<n; i++)
      y[i] = log(-log(1.0 - x[i]));
  }
  else if (link == GLM_LINK_IDENTITY)
  {
    for (i=0; i<n; i++)
      y[i] = x[i];
  }
  else if (link == GLM_LINK_INVERSE)
  {
    for (i=0; i<n; i++)
      y[i] = 1.0 / x[i];
  }
  else if (link == GLM_LINK_LOG)
  {
    for (i=0; i<n; i++)
      y[i] = log(x[i]);
  }
  else if (link == GLM_LINK_LOGIT)
  {
    for (i=0; i<n; i++)
      y[i] = log(x[i] / (1.0 - x[i]));
  }
  else if (link == GLM_LINK_SQRT)
  {
    for (i=0; i<n; i++)
      y[i] = sqrt(x[i]);
  }
  else if (link == GLM_LINK_PROBIT)
  {
    for (i=0; i<n; i++)
      y[i] = qnorm(x[i], 0.0, 1.0, true, false);
  }
  else if (link == GLM_LINK_CAUCHIT)
  {
    for (i=0; i<n; i++)
      y[i] = qcauchy(x[i], 0.0, 1.0, true, false);
  }
  else if (link == GLM_LINK_INVERSESQUARE)
  {
    for (i=0; i<n; i++)
      y[i] = 1.0 / (x[i] * x[i]);
  }
}



static inline void glm_linkinv(const int link, const int n, const double *const restrict x, double *const restrict y)
{
  int i;
  
  if (link == GLM_LINK_CLOGLOG)
  {
    for (i=0; i<n; i++)
      y[i] = -exp(-exp(x[i])) + 1.0;
  }
  else if (link == GLM_LINK_IDENTITY)
  {
    for (i=0; i<n; i++)
      y[i] = x[i];
  }
  else if (link == GLM_LINK_INVERSE)
  {
    for (i=0; i<n; i++)
    {
      if (x[i] > 0.0)
        y[i] = 1.0 / x[i];
      else
        y[i] = 0.0;
    }
  }
  else if (link == GLM_LINK_LOG)
  {
    for (i=0; i<n; i++)
      y[i] = exp(x[i]);
  }
  else if (link == GLM_LINK_LOGIT)
  {
    for (i=0; i<n; i++)
    {
      const double tmp = exp(x[i]);
      y[i] = tmp / (1.0 + tmp);
    }
  }
  else if (link == GLM_LINK_SQRT)
  {
    for (i=0; i<n; i++)
      y[i] = x[i] * x[i];
  }
  else if (link == GLM_LINK_PROBIT)
  {
    for (i=0; i<n; i++)
      y[i] = pnorm(x[i], 0.0, 1.0, true, false);
  }
  else if (link == GLM_LINK_CAUCHIT)
  {
    for (i=0; i<n; i++)
      y[i] = pcauchy(x[i], 0.0, 1.0, true, false);
      end do
    !$omp end do
  }
  else if (link == GLM_LINK_INVERSESQUARE)
  {
    for (i=0; i<n; i++)
      y[i] = 1.0 / sqrt(x[i]);
  }
}



static inline void glm_linkinv_deriv(const int link, const int n, const double *const restrict x, double *const restrict y)
{
  int i;
  
  if (link == GLM_LINK_CLOGLOG)
  {
    for (i=0; i<n; i++)
    {
      const double tmp = exp(x[i]);
      y[i] = tmp * exp(-tmp);
    }
  }
  else if (link == GLM_LINK_IDENTITY)
  {
    for (i=0; i<n; i++)
      y[i] = 1.0;
  }
  else if (link == GLM_LINK_INVERSE)
  {
    for (i=0; i<n; i++)
      y[i] = -1.0 / (x[i] * x[i]);
  }
  else if (link == GLM_LINK_LOG)
  {
    for (i=0; i<n; i++)
      y[i] = exp(x[i]);
  }
  else if (link == GLM_LINK_LOGIT)
  {
    for (i=0; i<n; i++)
    {
      const double tmp = exp(x[i]);
      y[i] = tmp / ((1.0 + tmp)*(1.0 + tmp));
    }
  }
  else if (link == GLM_LINK_SQRT)
  {
    for (i=0; i<n; i++)
      y[i] = 2.0 * x[i];
  }
  else if (link == GLM_LINK_PROBIT)
  {
    for (i=0; i<n; i++)
      y[i] = dnorm(x[i], 0.0, 1.0, false);
  }
  else if (link == GLM_LINK_CAUCHIT)
  {
    for (i=0; i<n; i++)
      y[i] = dcauchy(x[i], 0.0, 1.0, false);
  }
  else if (link == GLM_LINK_INVERSESQUARE)
  {
    for (i=0; i<n; i++)
    {
      const double tmp = sqrt(x[i]);
      y[i] = -0.5 / (tmp * sqrt(tmp));
    }
  }
}



// "working" residuals
static inline void glm_residuals(const int link, const int n, const double *const restrict y, const double *const restrict mu, const double *const restrict eta, double *const restrict resids)
{
  int i;
  
  if (link == GLM_LINK_CLOGLOG)
  {
    for (i=0; i<n; i++)
    {
      const double tmp = exp(eta[i]);
      resids[i] = (y[i] - mu[i] / (tmp * exp(-tmp));
    }
  }
  else if (link == GLM_LINK_IDENTITY)
  {
    for (i=0; i<n; i++)
      resids[i] = y[i] - mu[i];
  }
  else if (link == GLM_LINK_INVERSE)
  {
    for (i=0; i<n; i++)
      resids[i] = -1.0 * eta[i]*eta[i] * (y[i] - mu[i]);
  }
  else if (link == GLM_LINK_LOG)
  {
    for (i=0; i<n; i++)
      resids[i] = y[i]/mu[i] - 1.0;
  }
  else if (link == GLM_LINK_LOGIT)
  {
    for (i=0; i<n; i++)
      resids[i] = (y[i]/mu[i] - 1.0) / (1.0 - mu[i]);
  }
  else if (link == GLM_LINK_SQRT)
  {
    for (i=0; i<n; i++)
      resids[i] = (y[i] - mu[i]) / (2.0 * eta[i]);
  }
  else if (link == GLM_LINK_PROBIT)
  {
    resids[i] = (y[i] - mu[i]) / dnorm(et[i], 0.0, 1.0, false);
  }
  else if (link == GLM_LINK_PROBIT)
  {
    for (i=0; i<n; i++)
      resids[i] = (y[i] - mu[i]) / dcauchy(eta[i], 0.0, 1.0, false);
  }
  else if (link == GLM_LINK_INVERSESQUARE)
  {
    for (i=0; i<n; i++)
      resids(i) = -0.5 * (y[i] - mu[i]) / (eta[i] * sqrt(eta[i]));
  }
}


#endif
