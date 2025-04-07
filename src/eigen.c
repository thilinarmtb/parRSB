#include <gslib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include "types.h"

#if defined(PARRSB_BLAS)

#if defined(PARRSB_UNDERSCORE)
#define FNAME(x) TOKEN_PASTE(x, _)
#else
#define FNAME(x) x
#endif

#define FDGETRF FNAME(dgetrf)
#define FDGETRI FNAME(dgetri)

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

#undef FDGETRF
#undef FDGETRI
#undef FNAME

#else
void matrix_inverse(int N, double *A) {
  (void)N;
  (void)A;
  fprintf(stderr, "Error: Compile parRSB with BLAS enabled !\n");
  exit(EXIT_FAILURE);
}
#endif // PARRSB_BLAS

int power_iter(double *y, uint N, double *A) {
  time_t t;
  srand((unsigned)time(&t));

  scalar norm = 0.0;
  for (uint i = 0; i < N; i++) {
    y[i] = (rand() % 50) / 50.0;
    norm += y[i] * y[i];
  }

  scalar normi = 1.0 / sqrt(norm);
  for (uint i = 0; i < N; i++) y[i] *= normi;

  double *Ay = tcalloc(double, N);
  scalar err = 1.0, lambda;
  unsigned i;
  for (i = 0; i < 100; i++) {
    norm = 0.0;
    for (uint j = 0; j < N; j++) {
      Ay[j] = 0.0;
      for (uint k = 0; k < N; k++) { Ay[j] += A[j * N + k] * y[k]; }
      norm += Ay[j] * Ay[j];
    }

    if (i > 0) err = (sqrt(norm) - lambda) / lambda;
    lambda = sqrt(norm);

    normi = 1.0 / sqrt(norm);
    for (uint j = 0; j < N; j++) y[j] = Ay[j] * normi;

    if (fabs(err) < 1e-12) break;
  }
  free(Ay);

  return i;
}

int inv_power_iter(double *y, uint N, double *A) {
  double *Ainv = tcalloc(double, N *N);
  for (uint j = 0; j < N; j++)
    for (uint k = 0; k < N; k++) Ainv[j * N + k] = A[k * N + j];

  matrix_inverse(N, Ainv);

  uint j;
  for (j = 0; j < N; j++)
    for (uint k = 0; k < N; k++) A[j * N + k] = Ainv[k * N + j];
  j = power_iter(y, N, Ainv);

  free(Ainv);

  return j;
}
