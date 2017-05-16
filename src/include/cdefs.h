#ifndef R_LINMOD2_CDEFS_H_
#define R_LINMOD2_CDEFS_H_


#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define THROW_MEMERR error("unable to allocate necessary memory")
#define CHECKMALLOC(ptr) if (ptr == NULL) return -1;

#define THROW_LAPACKERR(info) error("LAPACK returned error code %d", info)
#define CHECKINFO(info) if (info!=0) THROW_LAPACKERR(info)

#define FREE(ptr) if(ptr!=NULL) free(ptr)


#endif
