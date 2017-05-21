#ifndef __LINMOD2_LM_TYPES__
#define __LINMOD2_LM_TYPES__


typedef struct linmodel
{
  int m;
  int n;
  int max_mn;
  int nrhs;
  double *restrict x;
  double *restrict oldx; // keeping the copy around improves performance a bit
  const double *restrict y;
  double *restrict coef;
  double *restrict resid;
  double *restrict fttd;
  double *restrict eff;
  double *restrict qraux; // NOTE set to NULL if you don't want the return
  int info;
} linmodel_t;


#endif
