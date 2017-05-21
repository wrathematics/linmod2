#include "linmod2.h"
#include "include/linmod2/lm.h"

SEXP R_dot_lm_fit_minimal(SEXP x, SEXP y, SEXP intercept)
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
  
  
  memset(DBLP(coef), 0, n*nrhs*sizeof(double));
  // memset(DBLP(resid), 0, m*nrhs*sizeof(double));
  memset(DBLP(fttd), 0, m*nrhs*sizeof(double));
  // memset(DBLP(eff), 0, m*nrhs*sizeof(double));
  
  
  INT(incpt_) = incpt;
  
  x_cp = malloc(m*n*sizeof(*x_cp));
  if (x_cp == NULL)
    THROW_MEMERR;
  
  if (incpt)
  {
    for (int i=0; i<m; i++)
      x_cp[i] = 1.0;
    
    memcpy(x_cp + m, DBLP(x), m*(n-1)*sizeof(*x_cp));
  }
  else
    memcpy(x_cp, DBLP(x), m*n*sizeof(*x_cp));
  
  lm.oldx = REAL(x);
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
  lm.qraux = NULL;
  
  int check = lm_fit(&lm);
  free(x_cp);
  if (check == LINMOD2_ERR_MALLOC)
    THROW_MEMERR;
  CHECKINFO(lm.info);
  
  ret_names = make_list_names(5, "coefficients", "fitted.values", "residuals", "effects", "intercept");
  ret = make_list(ret_names, 5, coef, fttd, resid, eff, incpt_);
  
  R_END;
  return ret;
}




SEXP R_dot_lm_fit(SEXP x, SEXP y, SEXP intercept)
{
  SEXP ret, ret_names;
  SEXP qr, coef, resid, eff;
  SEXP rank, pivot, qraux, tol, pivoted;
  double *fttd;
  linmodel_t lm;
  
  CHECK_IS_MATRIX(x);
  CHECK_IS_NUMERIC(x);
  CHECK_IS_NUMERIC(y);
  CHECK_IS_FLAG(intercept);
  
  const int incpt = INT(intercept);
  const int m = nrows(x);
  const int n = incpt ? ncols(x)+1 : ncols(x);
  const int nrhs = ncols(y);
  
  newRmat(qr, m, n, "dbl");
  newRmat(coef, n, nrhs, "dbl");
  newRmat(resid, m, nrhs, "dbl");
  newRmat(eff, m, nrhs, "dbl");
  fttd = (double*) R_alloc(m*n, sizeof(double));
  
  
  newRvec(rank, 1, "int");
  newRvec(pivot, n, "int");
  newRvec(qraux, n, "dbl");
  newRvec(tol, 1, "dbl");
  newRvec(pivoted, 1, "lgl");
  
  INT(rank) = n;
  DBL(tol) = 1e-8;
  INT(pivoted) = 0;
  
  for (int i=0; i<n; i++)
    INT(pivot, i) = i+1;
  
  memset(DBLP(coef), 0, n*nrhs*sizeof(double));
  // memset(DBLP(resid), 0, m*nrhs*sizeof(double));
  memset(fttd, 0, m*nrhs*sizeof(double));
  // memset(DBLP(eff), 0, m*nrhs*sizeof(double));
  
  double *const restrict x_cp = REAL(qr);
  
  if (incpt)
  {
    for (int i=0; i<m; i++)
      x_cp[i] = 1.0;
    
    memcpy(x_cp + m, DBLP(x), m*(n-1)*sizeof(*x_cp));
  }
  else
    memcpy(x_cp, DBLP(x), m*n*sizeof(*x_cp));
  
  lm.oldx = REAL(x);
  lm.m = m;
  lm.n = n;
  lm.max_mn = MAX(m, n);
  lm.nrhs = nrhs;
  lm.x = x_cp;
  lm.y = DBLP(y);
  lm.coef = DBLP(coef);
  lm.resid = DBLP(resid);
  lm.fttd = fttd;
  lm.eff = DBLP(eff);
  lm.qraux = DBLP(qraux);
  
  int check = lm_fit(&lm);
  if (check == LINMOD2_ERR_MALLOC)
    THROW_MEMERR;
  CHECKINFO(lm.info);
  
  ret_names = make_list_names(9, "qr", "coefficients", "residuals", "effects", "rank", "pivot", "qraux", "tol", "pivoted");
  ret = make_list(ret_names, 9, qr, coef, resid, eff, rank, pivot, qraux, tol, pivoted);
  
  R_END;
  return ret;
}



#define setDimNames(DIMNAMES, ROWNAMES, X) \
  newRlist(DIMNAMES, 2); \
  SET_VECTOR_ELT(DIMNAMES, 0, ROWNAMES); \
  SET_VECTOR_ELT(DIMNAMES, 1, RNULL); \
  setAttrib(X, R_DimNamesSymbol, DIMNAMES);

SEXP R_lm_fit(SEXP x, SEXP y, SEXP intercept)
{
  SEXP ret, ret_names;
  SEXP coef, resid, eff, fttd;
  SEXP coef_names, coef_dimnames, eff_names, eff_dimnames;
  SEXP QR, QR_names, QR_class;
  SEXP class = install("class");
  SEXP qr, qraux, pivot, tol, rank;
  SEXP df_resid;
  SEXP assign = RNULL;
  linmodel_t lm;
  
  CHECK_IS_MATRIX(x);
  CHECK_IS_NUMERIC(x);
  CHECK_IS_NUMERIC(y);
  CHECK_IS_FLAG(intercept);
  
  const int incpt = INT(intercept);
  const int m = nrows(x);
  const int n = incpt ? ncols(x)+1 : ncols(x);
  const int nrhs = ncols(y);
  
  newRmat(qr, m, n, "dbl");
  newRmat(coef, n, nrhs, "dbl");
  newRmat(fttd, m, nrhs, "dbl");
  newRmat(resid, m, nrhs, "dbl");
  newRmat(eff, m, nrhs, "dbl");
  
  newRvec(rank, 1, "int");
  newRvec(pivot, n, "int");
  newRvec(qraux, n, "dbl");
  newRvec(tol, 1, "dbl");
  
  newRvec(df_resid, 1, "int");
  newRvec(QR_class, 1, "str");
  
  INT(rank) = MIN(m, n);
  DBL(tol) = 1e-8;
  INT(df_resid) = m - INT(rank);
  
  for (int i=0; i<n; i++)
    INT(pivot, i) = i+1;
  
  memset(DBLP(coef), 0, n*nrhs*sizeof(double));
  // memset(DBLP(resid), 0, m*nrhs*sizeof(double));
  memset(REAL(fttd), 0, m*nrhs*sizeof(double));
  // memset(DBLP(eff), 0, m*nrhs*sizeof(double));
  
  double *const restrict x_cp = REAL(qr);
  
  if (incpt)
  {
    for (int i=0; i<m; i++)
      x_cp[i] = 1.0;
    
    memcpy(x_cp + m, DBLP(x), m*(n-1)*sizeof(*x_cp));
  }
  else
    memcpy(x_cp, DBLP(x), m*n*sizeof(*x_cp));
  
  lm.oldx = REAL(x);
  lm.m = m;
  lm.n = n;
  lm.max_mn = MAX(m, n);
  lm.nrhs = nrhs;
  lm.x = x_cp;
  lm.y = DBLP(y);
  lm.coef = DBLP(coef);
  lm.resid = DBLP(resid);
  lm.fttd = REAL(fttd);
  lm.eff = DBLP(eff);
  lm.qraux = DBLP(qraux);
  
  int check = lm_fit(&lm);
  if (check == LINMOD2_ERR_MALLOC)
    THROW_MEMERR;
  CHECKINFO(lm.info);
  
  QR_names = make_list_names(5, "qr", "qraux", "pivot", "tol", "rank");
  QR = make_list(QR_names, 5, qr, qraux, pivot, tol, rank);
  SET_STRING_ELT(QR_class, 0, mkChar("qr"));
  setAttrib(QR, class, QR_class);
  
  // set_NA_coef(n, nrhs, INT(rank), REAL(coef));
  
  coef_names = make_lmfit_default_coefnames(n);
  eff_names = make_lmfit_default_effectnames(m, n, INT(rank), INTP(pivot));
  setDimNames(coef_dimnames, coef_names, coef);
  setDimNames(eff_dimnames, eff_names, eff);
  
  ret_names = make_list_names(8, "coefficients", "residuals", "effects", "rank", "fitted.values", "assign", "qr", "df.residual");
  ret = make_list(ret_names, 8, coef, resid, eff, rank, fttd, assign, QR, df_resid);
  
  R_END;
  return ret;
}
