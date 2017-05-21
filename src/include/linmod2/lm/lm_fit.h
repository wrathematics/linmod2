#ifndef __LINMOD2_LM_FIT__
#define __LINMOD2_LM_FIT__

#include "types.h"

#include "../internal/cdefs.h"
#include "../internal/lapack.h"
#include "../internal/mmult.h"
#include "../internal/returns.h"
#include "../internal/types.h"


static inline int lm_init(const int ret_qraux, const int m, const int n, const int nrhs, double *restrict x, double *restrict y, linmodel_t *restrict lm)
{
  double *coef, *resid, *fttd, *eff, *qraux;
  
  LINMOD_TRY_MALLOC(1, lm);
  LINMOD_TRY_MALLOC(n*nrhs, coef);
  LINMOD_TRY_MALLOC(m*nrhs, resid);
  LINMOD_TRY_MALLOC(m*nrhs, fttd);
  LINMOD_TRY_MALLOC(m*nrhs, eff);
  
  if (ret_qraux)
    LINMOD_TRY_MALLOC(n, qraux);
  else
    qraux = NULL;
  
  lm->oldx = NULL;
  
  lm->m = m;
  lm->n = n;
  lm->max_mn = MAX(m, n);
  lm->nrhs = nrhs;
  lm->x = x;
  lm->y = y;
  lm->coef = coef;
  lm->resid = resid;
  lm->fttd = fttd;
  lm->eff = eff;
  lm->qraux = qraux;
  lm->info = 0;
  
  return LINMOD2_OK;
  
  
  
OOM:
  FREE(lm);
  FREE(coef);
  FREE(resid);
  FREE(fttd);
  FREE(eff);
  
  return LINMOD2_ERR_MALLOC;
}



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
    return LINMOD2_ERR_LAPACK;
  
  lwork = (int) tmp;
  
  if (lm->m > lm->n)
    dormqr_(&side, &trans, &lm->m, &lm->nrhs, &lm->n, &(double){0}, &lm->m, &(double){0}, &(double){0}, &lm->m, &tmp, &(int){-1}, &lm->info);
  else
    dormlq_(&side, &trans, &lm->m, &lm->nrhs, &lm->m, &(double){0}, &lm->m, &(double){0}, &(double){0}, &lm->m, &tmp, &(int){-1}, &lm->info);
  
  if (lm->info)
    return LINMOD2_ERR_LAPACK;
  
  lwork = MAX((int) tmp, lwork);
  *work = malloc(lwork * sizeof(**work));
  
  if (work == NULL)
    return LINMOD2_ERR_MALLOC;
  
  return lwork;
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
        lm->coef[i + n*j] = y[i + m*j];
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
  const int nrhs = lm->nrhs;
  
  
  if (lm->oldx == NULL)
  {
    // we don't have x anymore so use the X=QR or X=LQ factorization in its place
    if (m >= n)
    {
      const char uplo = 'U';
      
      copymat(m, n, nrhs, lm->coef, lm->fttd);
      memset(lm->fttd, 0, m*nrhs*sizeof(*lm->fttd));
      for (int j=0; j<nrhs; j++)
      {
        for (int i=0; i<n; i++)
        lm->fttd[i + m*j] = lm->coef[i + n*j];
      }
      
      // fttd = x*coef = Q*R*coef
      dtrmm_(&side, &uplo, &trans, &diag, &n, &nrhs, &alpha, lm->x, &m, lm->fttd, &m);
      dormqr_(&side, &trans, &m, &nrhs, &n, lm->x, &m, work, lm->fttd, &m, work+n, &lwork, &lm->info);
    }
    else
    {
      const char uplo = 'L';
      
      // fttd = x*coef = L*Q*coef
      //FIXME lwork probably too small
      for (int j=0; j<nrhs; j++)
      {
        memcpy(work+lm->m, lm->coef, lm->n*sizeof(*work));
        dormlq_(&side, &trans, &n, &nrhs, &m, lm->x, &m, work, work+m, &lm->n, work+n+m, &lwork, &lm->info);
        memcpy(lm->fttd + lm->m * j, work+lm->m, lm->m*sizeof(*work));
        dtrmm_(&side, &uplo, &trans, &diag, &m, &nrhs, &alpha, lm->x, &m, lm->fttd, &m);
      }
    }
  }
  else
    matmult(false, false, 1.0, lm->m, lm->n, lm->oldx, lm->n, lm->nrhs, lm->coef, lm->fttd);
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
  double *work = NULL;
  
  int lwork = lm_fit_workspace(lm, &work);
  if (lwork == LINMOD2_ERR_MALLOC)
    return LINMOD2_ERR_MALLOC;
  
  // fit y ~ x
  lm_compute_coefficients(lm, work, lwork);
  
  // effects = R*coef
  // lm_compute_effects(lm, work, lwork);
  
  // fitted values = X*coef
  lm_compute_fitted(lm, work, lwork);
  
  // residuals = y - fitted
  lm_compute_residuals(lm);
  
  if (lm->qraux != NULL)
  {
    for (int i=0; i<lm->n; i++)
      (lm->qraux)[i] = work[i];
  }
  
  free(work);
  
  if (lm->info != 0)
    THROW_LAPACKERR(lm->info);
  
  return 0;
}


#endif
