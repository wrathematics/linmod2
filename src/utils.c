// Copyright 2014, 2017 Drew Schmidt

#include "linmod2.h"
#include <math.h>

// NOTE: ceil(log10(2^31-1)) = 10
#define BUFLEN 16

// buflen = ((int) ceil(log10((double)n))) + 1;
static char buf[BUFLEN];

SEXP make_lmfit_default_coefnames(const int n)
{
  int i;
  int wordlen;
  SEXP ret;
  
  buf[0] = 'x';
  
  newRvec(ret, n, "str");
  
  for (i=0; i<n; i++)
  {
    sprintf(buf+1, "%d", i+1);
    wordlen = ((int) ceil(log10((double)i+2))) + 1;
    wordlen = MAX(wordlen, 2);
    SET_STRING_ELT(ret, i, mkCharLen(buf, wordlen));
  }
  
  return ret;
}



SEXP make_lmfit_default_effectnames(const int m, const int n, const int rank, const int *pvt)
{
  int i, j;
  int wordlen;
  int maxpvt = rank; //MIN(m, n);
  SEXP ret;
  
  for (i=1; i<n; i++)
  {
    if (pvt[i] < pvt[i-1])
    {
      maxpvt = i;
      break;
    }
  }
  
  buf[0] = 'x';
  
  newRvec(ret, m, "str");
  
  for (i=0; i<maxpvt; i++)
  {
    j = pvt[i];
    sprintf(buf+1, "%d", j);
    wordlen = ((int) ceil(log10((double)j+1))) + 1;
    wordlen = MAX(wordlen, 2);
    SET_STRING_ELT(ret, i, mkCharLen(buf, wordlen));
  }
  
  for (i=maxpvt; i<m; i++)
    SET_STRING_ELT(ret, i, mkChar(""));
  
  return ret;
}



void set_NA_coef(const int n, const int nrhs, const int rank, double *const restrict coef)
{
  int i, j;
  
  if (rank < n)
  {
    for (j=0; j<nrhs; j++)
    {
      for (i=rank; i<n; i++)
        coef[i + n*j] = NA_REAL;
    }
  }
}
