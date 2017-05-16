#ifndef __LINMOD2_CDEFS_H__
#define __LINMOD2_CDEFS_H__


#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define CHECKMALLOC(ptr) if (ptr == NULL) return -1;
#define FREE(ptr) if(ptr!=NULL) free(ptr)

#define LINMOD_TRY_MALLOC(n, x) \
  x = malloc(n * sizeof(*x)); \
  if (x == NULL) goto OOM


#endif
