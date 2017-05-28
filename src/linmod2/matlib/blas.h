#ifndef __LINMOD2_MATLIB_BLAS__
#define __LINMOD2_MATLIB_BLAS__

#include "../internal/types.h"


void dgemm_(cchar_r transa, cchar_r transb, cint_r m, cint_r n, cint_r k,
  cdbl_r alpha, cdbl_r a, cint_r lda, cdbl_r b, cint_r ldb, cdbl_r beta,
  dbl_r c, cint_r ldc);


#endif
