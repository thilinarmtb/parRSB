#include "parrsb_impl.h"

#include <math.h>
#include <stdio.h>
#include <time.h>

int power_iter(double *y, uint N, double *A) {
  srand((unsigned)time(NULL));

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

  return 0;
}

int inv_power_iter(double *y, uint N, double *A) {
  double *Ainv = tcalloc(double, N *N);
  for (uint j = 0; j < N; j++)
    for (uint k = 0; k < N; k++) Ainv[j * N + k] = A[k * N + j];

  int err = matrix_inverse(N, Ainv);
  if (err != 0) {
    fprintf(stderr, "Matrix inverse failed!\n");
    fflush(stderr);
    return err;
  }

  for (uint j = 0; j < N; j++)
    for (uint k = 0; k < N; k++) A[j * N + k] = Ainv[k * N + j];
  power_iter(y, N, Ainv);

  free(Ainv);

  return 0;
}

inline static double sign(scalar a, scalar b) {
  scalar m = b >= 0.0 ? 1.0 : -1.0;
  return fabs(a) * m;
}

int tqli(scalar *evectors, scalar *evalues, sint n, scalar *diagonal,
         scalar *upper, int id) {
  if (n == 0) return 0;

  scalar *d = tcalloc(scalar, 2 * n), *e = d + n;
  sint i;
  for (i = 0; i < n; i++) d[i] = diagonal[i];
  for (i = 0; i < n - 1; i++) e[i] = upper[i];
  e[n - 1] = 0.0;

  for (i = 0; i < n; i++) {
    for (sint j = 0; j < n; j++) evectors[i * n + j] = 0;
    evectors[i * n + i] = 1;
  }

  sint j, k, l, iter, m;
  for (l = 0; l < n; l++) {
    iter = 0;
    do {
      for (m = l; m < n - 1; m++) {
        scalar dd = fabs(d[m]) + fabs(d[m + 1]);
        /* Should use a tolerance for this check */
        if (fabs(e[m]) / dd < SCALAR_TOL) break;
      }

      if (m != l) {
        if (iter++ == 30) {
          if (id == 0) printf("Too many iterations.\n");
          // vec_copy(*evalues, d);
          for (i = 0; i < n; i++) evalues[i] = d[i];
          return 1;
        }

        scalar g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        scalar r = sqrt(g * g + 1.0);

        g = d[m] - d[l] + e[l] / (g + sign(r, g));
        scalar s = 1.0, c = 1.0, p = 0.0;

        for (i = m - 1; i >= l; i--) {
          scalar f = s * e[i];
          scalar b = c * e[i];

          if (fabs(f) >= fabs(g)) {
            c = g / f;
            r = sqrt(c * c + 1.0);
            e[i + 1] = f * r;
            s = 1.0 / r;
            c = c * s;
          } else {
            s = f / g;
            r = sqrt(s * s + 1.0);
            e[i + 1] = g * r;
            c = 1.0 / r;
            s = s * c;
          }

          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          /* Find eigenvectors */
          for (k = 0; k < n; k++) {
            f = evectors[k * n + i + 1];
            evectors[k * n + i + 1] = s * evectors[k * n + i] + c * f;
            evectors[k * n + i] = c * evectors[k * n + i] - s * f;
          }
          /* Done with eigenvectors */
        }

        if (r < SCALAR_TOL && i >= l) continue;

        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }

  /* Orthnormalize eigenvectors -- Just normalize? */
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      scalar tmp = evectors[i * n + j];
      evectors[i * n + j] = evectors[j * n + i];
      evectors[j * n + i] = tmp;
    }
  }

  for (k = 0; k < n; k++) {
    e[k] = 0;
    for (sint i = 0; i < n; i++)
      e[k] += evectors[k * n + i] * evectors[k * n + i];
    if (e[k] > 0.0) e[k] = sqrt(fabs(e[k]));
    scalar scale = 1.0 / e[k];
    for (sint i = 0; i < n; i++) evectors[k * n + i] *= scale;
  }

  // vec_copy(*evalues, d);
  for (i = 0; i < n; i++) evalues[i] = d[i];

  free(d);

  return 0;
}
