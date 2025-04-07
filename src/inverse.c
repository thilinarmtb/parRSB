#include <gslib.h>
#include <stdlib.h>

#include "types.h"

#if defined(PARRSB_BLAS)

#define fdgetrf GS_FORTRAN_UNPREFIXED(dgetrf, DGETRF)
void fdgetrf(int *M, int *N, double *A, int *lda, int *IPIV, int *info);

#define fdgetri GS_FORTRAN_UNPREFIXED(dgetri, DGETRI)
void fdgetri(int *N, double *A, int *lda, int *ipiv, double *work, int *lwork,
             int *info);

#define fsgetrf GS_FORTRAN_UNPREFIXED(sgetrf, SGETRF)
void fsgetrf(int *M, int *N, float *A, int *lda, int *IPIV, int *info);

#define fsgetri GS_FORTRAN_UNPREFIXED(sgetri, SGETRI)
void fsgetri(int *N, float *A, int *lda, int *ipiv, float *work, int *lwork,
             int *info);

int matrix_inverse(int N, scalar *A) {
  int size = N * N;
  int info;

  int *ipiv = (int *)calloc(N, sizeof(int));
  scalar *work = (scalar *)calloc(N * N, sizeof(scalar));

  if (sizeof(scalar) == sizeof(double))
    fdgetrf(&N, &N, A, &N, ipiv, &info);
  else if (sizeof(scalar) == sizeof(float))
    fsgetrf(&N, &N, (float *)A, &N, ipiv, &info);
  else
    info = N + 1;

  if (info != 0) goto cleanup;

  if (sizeof(scalar) == sizeof(double))
    fdgetri(&N, A, &N, ipiv, work, &size, &info);
  else if (sizeof(scalar) == sizeof(float))
    fsgetri(&N, (float *)A, &N, ipiv, (float *)work, &size, &info);
  else
    info = N + 1;

cleanup:
  free(ipiv);
  free(work);
  return info;
}

#else
int matrix_inverse(int N, scalar *A) {
  (void)A;
  fprintf(stderr, "Error: Compile parRSB with BLAS enabled!\n");
  return N + 1;
}
#endif // PARRSB_BLAS
