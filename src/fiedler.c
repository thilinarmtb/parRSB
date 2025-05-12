#include "metrics.h"
#include "parrsb_impl.h"
#include "sort.h"

#include <math.h>
#include <time.h>

inline static scalar dot(const scalar *y, const scalar *x, uint n) {
  scalar result = 0.0;
  for (uint i = 0; i < n; i++) result += x[i] * y[i];
  return result;
}

inline static scalar norm2(const scalar *r, uint n, const struct comm *c) {
  scalar wrk;
  scalar norm = dot(r, r, n);
  comm_allreduce(c, gs_scalar, gs_add, &norm, 1, &wrk);
  return sqrt(norm);
}

inline static void ortho(scalar *q, uint n, ulong ng, const struct comm *c) {
  if (ng == 0) return;

  scalar sum = 0.0;
  for (uint i = 0; i < n; i++) sum += q[i];

  scalar wrk;
  comm_allreduce(c, gs_scalar, gs_add, &sum, 1, &wrk);
  sum /= ng;

  for (uint i = 0; i < n; i++) q[i] -= sum;
}

static int lanczos_aux(scalar *diag, scalar *upper, scalar *rr, laplacian l,
                       const scalar *f, int niter, double tol,
                       const struct comm *c) {
  uint n = laplacian_size(l);
  slong ng = n, wrk;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, &wrk);

  scalar *r = tcalloc(scalar, n);
  scalar *p = tcalloc(scalar, n);
  scalar *w = tcalloc(scalar, n);

  for (uint i = 0; i < n; i++) r[i] = f[i];
  ortho(r, n, ng, c);

  scalar rtz1 = 1, pap = 0, alpha, beta, rtz2, pap0;
  scalar rtr = dot(r, r, n), buf[2];
  comm_allreduce(c, gs_scalar, gs_add, &rtr, 1, buf);
  scalar rnorm = sqrt(rtr), rtol = rnorm * tol;

  metric_set(TOL_INIT, rnorm);
  metric_set(TOL_TGT, rtol);

  scalar rni = 1.0 / rnorm;
  for (uint i = 0; i < n; i++) rr[0 * n + i] = r[i] * rni;

  int iter = 0;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1, rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0) beta = 0.0;

    for (uint i = 0; i < n; i++) p[i] = beta * p[i] + r[i];
    ortho(p, n, ng, c);

    laplacian_op(w, l, p);

    pap0 = pap;
    pap = dot(w, p, n);
    comm_allreduce(c, gs_scalar, gs_add, &pap, 1, buf);

    alpha = rtz1 / pap;
    for (uint i = 0; i < n; i++) r[i] = r[i] - alpha * w[i];

    rtr = dot(r, r, n);
    comm_allreduce(c, gs_scalar, gs_add, &rtr, 1, buf);
    rnorm = sqrt(rtr), rni = 1.0 / rnorm;

    for (uint i = 0; i < n; i++) rr[(iter + 1) * n + i] = r[i] * rni;

    if (iter == 0) {
      diag[iter] = pap / rtz1;
    } else {
      diag[iter] = (beta * beta * pap0 + pap) / rtz1;
      upper[iter - 1] = -beta * pap0 / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      iter++;
      break;
    }
  }
  metric_set(TOL_FNL, rnorm);

  free(r), free(p), free(w);
  return iter;
}

int lanczos(scalar *f, laplacian l, scalar *vi, const struct comm *c,
            const parrsb_options opts) {
  uint n = laplacian_size(l);
  slong ng = n, wrk;
  comm_allreduce(c, gs_long, gs_add, &ng, 1, &wrk);
  if (ng == 0) return 0;

  uint miter = opts->rsb_max_iter;
  if (miter > ng) miter = ng;

  scalar *alpha = tcalloc(scalar, miter);
  scalar *beta = tcalloc(scalar, miter - 1);
  scalar *evals = tcalloc(scalar, miter);
  scalar *evecs = tcalloc(scalar, miter * miter);
  scalar *rr = tcalloc(scalar, (miter + 1) * n);

  uint iter = miter, ipass = 0;
  for (; (iter == miter) && (ipass < (uint)opts->rsb_max_passes); ipass++) {
    metric_tic(c, RSB_LANCZOS);
    iter = lanczos_aux(alpha, beta, rr, l, vi, miter, opts->rsb_tol, c);
    metric_toc(c, RSB_LANCZOS);

    // Use TQLI and find the minimum eigenvalue and associated vector
    metric_tic(c, RSB_LANCZOS_TQLI);

    tqli(evecs, evals, iter, alpha, beta, c->id);
    scalar eval_min = fabs(evals[0]);
    uint eval_min_idx = 0;
    for (uint i = 1; i < iter; i++) {
      if (fabs(evals[i]) < eval_min) {
        eval_min = fabs(evals[i]);
        eval_min_idx = i;
      }
    }

    for (uint i = 0; i < n; i++) {
      f[i] = 0.0;
      for (uint j = 0; j < iter; j++)
        f[i] += rr[j * n + i] * evecs[eval_min_idx * iter + j];
    }
    ortho(f, n, ng, c);
    for (uint i = 0; i < n; i++) vi[i] = f[i];

    metric_toc(c, RSB_LANCZOS_TQLI);
  }

  free(alpha), free(beta), free(evals), free(evecs), free(rr);
  return (ipass - 1) * miter + iter;
}

static void set_rhs(scalar *r, uint n, const struct comm *c) {
  slong out[2][1], wrk[2][1], in = n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], ng = out[1][0];

  for (uint i = 0; i < n; i++) {
    r[i] = start + i + 1.0;
    if ((start + i) < ng / 2) r[i] += 1e3 * ng;
  }

  ortho(r, n, ng, c);

  scalar normi = 1.0 / norm2(r, n, c);
  for (uint i = 0; i < n; i++) r[i] *= normi;
}

int fiedler(scalar *f, laplacian l, const parrsb_options opts,
            const struct comm *c) {
  // Return if the number of processes is equal to 1.
  if (c->np == 1) return 0;
  // Or if rsb algorithm is set to something other than lanczos.
  if (opts->rsb_algo > 0) return 1;

  metric_tic(c, RSB_FIEDLER);

  uint n = laplacian_size(l);
  scalar *r = tcalloc(scalar, n);
  set_rhs(r, n, c);
  int iter = lanczos(f, l, r, c, opts);
  free(r);

  metric_acc(RSB_FIEDLER_CALC_NITER, iter);

  scalar normi = 1.0 / norm2(f, n, c);
  for (uint i = 0; i < n; i++) f[i] *= normi;

  metric_toc(c, RSB_FIEDLER);

  return 0;
}
