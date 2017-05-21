#ifndef R__LINMOD2_H__
#define R__LINMOD2_H__


#include "include/cdefs.h"
#include "include/reactor.h"
#include "include/RNACI.h"


// utils.c
SEXP make_lmfit_default_coefnames(const int n);
SEXP make_lmfit_default_effectnames(const int m, const int n, const int rank, const int *pvt);
void set_NA_coef(const int n, const int nrhs, const int rank, double *const restrict coef);


#endif
