#include "linmod2.h"


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
    dgeqrf_(&m, &n, lm->x, &m, work, work+n, &lwork, &lm->info);
    dormqr_(&(char){'L'}, &(char){'T'}, &m, &lm->nrhs, &n, lm->x, &m, work, y, &m, work+n, &lwork, &lm->info);
    dtrtrs_(&(char){'U'}, &(char){'N'}, &(char){'N'}, &n, &lm->nrhs, lm->x, &m, y, &m, &lm->info);
    
    for (int j=0; j<lm->nrhs; j++)
    {
      for (int i=0; i<n; i++)
        lm->coef[i + m*j] = lm->eff[i + m*j];
    }
  }
  // else // TODO
  // {
  //   
  // }
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
  
  
  copymat(lm->m, lm->n, lm->nrhs, lm->coef, lm->fttd);
  
  
  // we don't have x anymore so use the factorization in its place
  if (lm->m >= lm->n)
  {
    char uplo = 'U';
    
    work[0] = 0;
    // fttd = x*coef = Q*R*coef
    dtrmm_(&side, &uplo, &trans, &diag, &lm->n, &lm->nrhs, &alpha, lm->x, &lm->m, lm->fttd, &lm->m);
    dormqr_(&side, &trans, &lm->m, &lm->nrhs, &lm->n, lm->x, &lm->m, work, lm->fttd, &lm->m, work+lm->n, &lwork, &lm->info);
  }
  else
  {
    char uplo = 'L';
    
    // fttd = x*coef = L*Q*coef
    dormlq_(&side, &trans, &lm->m, &lm->nrhs, &lm->m, lm->x, &lm->m, work, lm->fttd, &lm->m, work+lm->m, &lwork, &lm->info);
    dtrmm_(&side, &uplo, &trans, &diag, &lm->m, &lm->nrhs, &alpha, lm->x, &lm->m, lm->fttd, &lm->m);
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
  lm_compute_effects(lm, work, lwork);
  
  // fitted values = X*coef
  lm_compute_fitted(lm, work, lwork);
  
  // residuals = y - fitted
  lm_compute_residuals(lm);
  
  
  free(work);
  
  if (lm->info != 0)
    THROW_LAPACKERR(lm->info);
  
  return 0;
}



SEXP R_lm_fit(SEXP x, SEXP y, SEXP intercept)
{
  SEXP ret, ret_names;
  SEXP coef, resid, fttd, eff, incpt_;
  linmodel_t lm;
  double *x_cp;
  
  CHECK_IS_MATRIX(x);
  CHECK_IS_NUMERIC(x);
  CHECK_IS_NUMERIC(y);
  CHECK_IS_FLAG(intercept);
  
  const int incpt = INT(intercept);
  const int m = nrows(x);
  const int n = incpt ? ncols(x)+1 : ncols(x);
  const int nrhs = ncols(y);
  
  newRmat(coef, n, nrhs, "dbl");
  newRmat(resid, m, nrhs, "dbl");
  newRmat(fttd, m, nrhs, "dbl");
  newRmat(eff, m, nrhs, "dbl");
  newRvec(incpt_, 1, "int");
  
  INT(incpt_) = incpt;
  
  x_cp = malloc(m*n*sizeof(*x_cp));
  CHECKMALLOC(x_cp);
  if (incpt)
  {
    for (int i=0; i<m; i++)
      x_cp[i] = 1.0;
    
    memcpy(x_cp + m, DBLP(x), m*(n-1)*sizeof(*x_cp));
  }
  else
    memcpy(x_cp, DBLP(x), m*n*sizeof(*x_cp));
  
  lm.m = m;
  lm.n = n;
  lm.max_mn = MAX(m, n);
  lm.nrhs = nrhs;
  lm.x = x_cp;
  lm.y = DBLP(y);
  lm.coef = DBLP(coef);
  lm.resid = DBLP(resid);
  lm.fttd = DBLP(fttd);
  lm.eff = DBLP(eff);
  
  int check = lm_fit(&lm);
  free(x_cp);
  if (check == -1)
    THROW_MEMERR;
  
  ret_names = make_list_names(5, "coefficients", "fitted.values", "residuals", "effects", "intercept");
  ret = make_list(ret_names, 5, coef, fttd, resid, eff, incpt_);
  
  R_END;
  return ret;
}
