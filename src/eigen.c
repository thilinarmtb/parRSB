#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "gslib.h"

#define FDGETRF GS_FORTRAN_UNPREFIXED(dgetrf, DGETRF)
#define FDGETRI GS_FORTRAN_UNPREFIXED(dgetri, DGETRI)

#if defined(PARRSB_BLAS)
void FDGETRF(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
void FDGETRI(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork,
             int *INFO);

void matrix_inverse(int N, double *A) {
  int size = N * N;
  int info;

  int *ipiv = (int *)calloc(N, sizeof(int));
  double *work = (double *)calloc(N * N, sizeof(double));

  FDGETRF(&N, &N, A, &N, ipiv, &info);
  if (info != 0) printf("dgetrf: %d\n", info);

  FDGETRI(&N, A, &N, ipiv, work, &size, &info);
  if (info != 0) printf("dgetri: %d\n", info);

  free(ipiv);
  free(work);
}
#else
void matrix_inverse(int N, double *A) {
  fprintf(stderr, "matrix_inverse: Build with BLAS support.\n");
  exit(EXIT_FAILURE);
}
#endif

#undef FDGETRF
#undef FDGETRI
