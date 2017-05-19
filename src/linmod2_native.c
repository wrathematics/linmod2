/* Automatically generated. Do not edit by hand. */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

extern SEXP R_dot_lm_fit_minimal(SEXP x, SEXP y, SEXP intercept);
extern SEXP R_dot_lm_fit(SEXP x, SEXP y, SEXP intercept);
extern SEXP R_lm_fit(SEXP x, SEXP y, SEXP intercept);

static const R_CallMethodDef CallEntries[] = {
  {"R_dot_lm_fit_minimal", (DL_FUNC) &R_dot_lm_fit_minimal, 3},
  {"R_dot_lm_fit", (DL_FUNC) &R_dot_lm_fit, 3},
  {"R_lm_fit", (DL_FUNC) &R_lm_fit, 3},
  {NULL, NULL, 0}
};

void R_init_lqr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
