#ifndef __LINMOD2_LM__
#define __LINMOD2_LM__

#include "internal/cdefs.h"
#include "internal/lapack.h"
#include "internal/types.h"


typedef struct linmodel
{
  int m;
  int n;
  int max_mn;
  int nrhs;
  double *restrict x;
  const double *restrict y;
  double *restrict coef;
  double *restrict resid;
  double *restrict fttd;
  double *restrict eff;
  int info;
} linmodel_t;


static inline void copymat(cint m_in, cint m_out, cint n, cdbl_r in, dbl_r out)
{
  const int min = MIN(m_in, m_out);
  
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<min; i++)
      out[i + m_out*j] = in[i + m_in*j];
    
    for (int i=min; i<m_out; i++)
      out[i + m_out*j] = 0.0;
  }
}



static inline int lm_fit_workspace(linmodel_t *const restrict lm, double **restrict work)
{
  int lwork;
  const char side = 'L';
  const char trans = 'N';
  double tmp = 0;
  
  dgels_(&trans, &lm->m, &lm->n, &lm->nrhs, &(double){0}, &lm->m, &(double){0}, &lm->max_mn, &tmp, &(int){-1}, &lm->info);
  
  if (lm->info)
    goto err_handler;
  
  lwork = (int) tmp;
  
  if (lm->m > lm->n)
    dormqr_(&side, &trans, &lm->m, &lm->nrhs, &lm->n, &(double){0}, &lm->m, &(double){0}, &(double){0}, &lm->m, &tmp, &(int){-1}, &lm->info);
  else
    dormlq_(&side, &trans, &lm->m, &lm->nrhs, &lm->m, &(double){0}, &lm->m, &(double){0}, &(double){0}, &lm->m, &tmp, &(int){-1}, &lm->info);
  
  if (lm->info)
    goto err_handler;
  
  lwork = MAX((int) tmp, lwork);
  *work = malloc(lwork * sizeof(**work));
  
  if (work == NULL)
    goto err_handler;
  
  return lwork;
  
  
err_handler:
  
  free(lm->x);
  if (lm->info != 0)
    THROW_LAPACKERR(lm->info);
  if (work == NULL)
    THROW_MEMERR;
  
  return -1;
}



static inline void lm_compute_coefficients(linmodel_t *const restrict lm, dbl_r work, cint lwork)
{
  const int m = lm->m;
  const int n = lm->n;
  double *restrict y;
  // const char trans = 'N';
  
  // keep y const, so store it in eff or coef (whichever is larger) for now
  if (m >= n)
    y = lm->eff;
  else
    y = lm->coef;
  
  for (int j=0; j<lm->nrhs; j++)
  {
    for (int i=0; i<m; i++)
      y[i + lm->max_mn*j] = lm->y[i + m*j];
  }
  
  // dgels_(&trans, &lm->m, &lm->n, &lm->nrhs, lm->x, &lm->m, y, &lm->max_mn, work, &lwork, &lm->info);
  if (m >= n)
  {
    // x = Q*R
    dgeqrf_(&m, &n, lm->x, &m, work, work+n, &lwork, &lm->info);
    // coef = R^(-1) * Q^T * y
    dormqr_(&(char){'L'}, &(char){'T'}, &m, &lm->nrhs, &n, lm->x, &m, work, y, &m, work+n, &lwork, &lm->info);
    dtrtrs_(&(char){'U'}, &(char){'N'}, &(char){'N'}, &n, &lm->nrhs, lm->x, &m, y, &m, &lm->info);
    
    for (int j=0; j<lm->nrhs; j++)
    {
      for (int i=0; i<n; i++)
        lm->coef[i + m*j] = lm->eff[i + m*j];
    }
  }
  else
  {
    // x = L*Q
    dgelqf_(&m, &n, lm->x, &m, work, work+m, &lwork, &lm->info);
    // coef = Q^T * L^(-1) * y
    dtrtrs_(&(char){'L'}, &(char){'N'}, &(char){'N'}, &m, &lm->nrhs, lm->x, &m, y, &n, &lm->info);
    dormlq_(&(char){'L'}, &(char){'T'}, &n, &lm->nrhs, &m, lm->x, &m, work, y, &n, work+m, &lwork, &lm->info);
  }
}



// TODO:
// m>=n : R*coef
// m<n  : 
static inline void lm_compute_effects(linmodel_t *const restrict lm, dbl_r work, cint lwork)
{
  const char side = 'L';
  const char trans = 'N';
  const int min_mn = MIN(lm->m, lm->n);
  
  // the solution is stored in 'eff'; copy over the first rank elements
  copymat(lm->m, lm->n, lm->nrhs, lm->coef, lm->eff);
  
  dormqr_(&side, &trans, &lm->m, &lm->nrhs, &lm->n, lm->x, &lm->m, work, lm->eff, &lm->m, work+min_mn, &lwork, &lm->info);
}



static inline void lm_compute_fitted(linmodel_t *const restrict lm, dbl_r work, cint lwork)
{
  const char side = 'L';
  const char trans = 'N';
  const char diag = 'N';
  const double alpha = 1.0;
  
  const int m = lm->m;
  const int n = lm->n;
  
  
  // we don't have x anymore so use the factorization in its place
  if (m >= n)
  {
    const char uplo = 'U';
    
    copymat(m, n, lm->nrhs, lm->coef, lm->fttd);
    
    // fttd = x*coef = Q*R*coef
    dtrmm_(&side, &uplo, &trans, &diag, &n, &lm->nrhs, &alpha, lm->x, &m, lm->fttd, &m);
    dormqr_(&side, &trans, &m, &lm->nrhs, &n, lm->x, &m, work, lm->fttd, &m, work+n, &lwork, &lm->info);
  }
  else
  {
    const char uplo = 'L';
    
    // fttd = x*coef = L*Q*coef
    //FIXME lwork probably too small
    for (int j=0; j<lm->nrhs; j++)
    {
      memcpy(work+lm->m, lm->coef, lm->n*sizeof(*work));
      dormlq_(&side, &trans, &n, &lm->nrhs, &m, lm->x, &m, work, work+m, &lm->n, work+n+m, &lwork, &lm->info);
      memcpy(lm->fttd + lm->m * j, work+lm->m, lm->m*sizeof(*work));
      dtrmm_(&side, &uplo, &trans, &diag, &m, &lm->nrhs, &alpha, lm->x, &m, lm->fttd, &m);
    }
  }
}



static inline void lm_compute_residuals(linmodel_t *const restrict lm)
{
  for (int j=0; j<lm->nrhs; j++)
  {
    for (int i=0; i<lm->m; i++)
      lm->resid[i + lm->m*j] = lm->y[i + lm->m*j] - lm->fttd[i + lm->m*j];
  }
}



// the main driver
static inline int lm_fit(linmodel_t *const restrict lm)
{
  double *work;
  
  int lwork = lm_fit_workspace(lm, &work);
  if (lwork == -1)
    return -1;
  
  // fit y ~ x
  lm_compute_coefficients(lm, work, lwork);
  
  // effects = R*coef
  // lm_compute_effects(lm, work, lwork);
  
  // fitted values = X*coef
  lm_compute_fitted(lm, work, lwork);
  
  // residuals = y - fitted
  lm_compute_residuals(lm);
  
  
  free(work);
  
  if (lm->info != 0)
    THROW_LAPACKERR(lm->info);
  
  return 0;
}


#endif
