#include "metrics.h"
#include "parrsb_impl.h"
#include "sort.h"

#include <math.h>
#include <time.h>

inline static scalar dot(scalar *y, scalar *x, uint n) {
  scalar result = 0.0;
  for (uint i = 0; i < n; i++) result += x[i] * y[i];
  return result;
}

inline static void ortho(scalar *q, uint lelt, ulong n, const struct comm *c) {
  scalar sum = 0.0;
  for (uint i = 0; i < lelt; i++) sum += q[i];

  scalar wrk;
  comm_allreduce(c, gs_scalar, gs_add, &sum, 1, &wrk);
  sum /= n;

  for (uint i = 0; i < lelt; i++) q[i] -= sum;
}

static int lanczos_aux(scalar *diag, scalar *upper, scalar *rr, uint n,
                       ulong nelg, int niter, double tol, scalar *f,
                       laplacian gl, const struct comm *c) {
  scalar *r = tcalloc(scalar, 3 * n);
  scalar *p = r + n;
  scalar *w = p + n;

  for (uint i = 0; i < n; i++) r[i] = f[i];
  ortho(r, n, nelg, c);

  scalar rtz1 = 1, pap = 0, alpha, beta, rtz2, pap_old;
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

    // add2s1(p,r,beta,n)
    for (uint i = 0; i < n; i++) p[i] = beta * p[i] + r[i];

    scalar pp = dot(p, p, n);
    comm_allreduce(c, gs_scalar, gs_add, &pp, 1, buf);

    // vec_ortho(c, p, nelg);
    ortho(p, n, nelg, c);

    laplacian_op(w, gl, p);

    pap_old = pap, pap = dot(w, p, n);
    comm_allreduce(c, gs_scalar, gs_add, &pap, 1, buf);

    alpha = rtz1 / pap;
    // vec_axpby(r, r, 1.0, w, -1.0 * alpha);
    for (uint i = 0; i < n; i++) r[i] = r[i] - alpha * w[i];

    rtr = dot(r, r, n);
    comm_allreduce(c, gs_scalar, gs_add, &rtr, 1, buf);
    rnorm = sqrt(rtr), rni = 1.0 / rnorm;

    // vec_scale(rr[iter + 1], r, rni);
    for (uint i = 0; i < n; i++) rr[(iter + 1) * n + i] = r[i] * rni;

    if (iter == 0) {
      diag[iter] = pap / rtz1;
    } else {
      diag[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      iter++;
      break;
    }
  }

  metric_set(TOL_FNL, rnorm);

  free(r);

  return iter;
}

static int lanczos(scalar *fiedler, laplacian l, scalar *initv,
                   const struct comm *c, const parrsb_options opts, slong nelg,
                   buffer *bfr) {
  uint miter = opts->rsb_max_iter;
  if (miter > nelg) miter = nelg;
  uint n = laplacian_size(l);

  scalar *alpha = tcalloc(scalar, 2 * miter - 1);
  scalar *beta = alpha + miter;
  scalar *rr = tcalloc(scalar, (miter + 1) * n);
  scalar *evecs = tcalloc(scalar, miter * miter);
  scalar *evals = tcalloc(scalar, miter);

  uint iter = miter, ipass = 0;
  for (; (iter == miter) && (ipass < (uint)opts->rsb_max_passes); ipass++) {
    metric_tic(c, RSB_LANCZOS);
    iter = lanczos_aux(alpha, beta, rr, n, nelg, miter, opts->rsb_tol, initv, l,
                       c);
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
      fiedler[i] = 0.0;
      for (uint j = 0; j < iter; j++)
        fiedler[i] += rr[j * n + i] * evecs[eval_min_idx * iter + j];
    }
    ortho(fiedler, n, nelg, c);
    for (uint i = 0; i < n; i++) initv[i] = fiedler[i];

    metric_toc(c, RSB_LANCZOS_TQLI);
  }

  free(alpha), free(rr), free(evecs), free(evals);

  return (ipass - 1) * miter + iter;
}

int fiedler(scalar *f, laplacian l, const parrsb_options opts,
            const struct comm *c, buffer *bfr) {
  // Return if the number of processes is equal to 1 or rsb algorithm is set to
  // something other than lanczos.
  if (c->np == 1 || opts->rsb_algo > 0) return 1;

  metric_tic(c, RSB_FIEDLER);

  uint n = laplacian_size(l);
  slong out[2][1], wrk[2][1], in = n;
  comm_scan(out, c, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], ng = out[1][0];

  scalar *vi = tcalloc(scalar, n);
  for (uint i = 0; i < n; i++) {
    vi[i] = start + i + 1.0;
    if (start + i < ng / 2) vi[i] += 1000 * ng;
  }

  ortho(vi, n, ng, c);
  scalar norm = dot(vi, vi, n);
  comm_allreduce(c, gs_scalar, gs_add, &norm, 1, wrk);
  scalar normi = 1.0 / sqrt(norm);
  for (uint i = 0; i < n; i++) vi[i] *= normi;

  int iter = 0;
  switch (opts->rsb_algo) {
  case 0: iter = lanczos(f, l, vi, c, opts, ng, bfr); break;
  default: break;
  }
  metric_acc(RSB_FIEDLER_CALC_NITER, iter);

  norm = dot(f, f, n);
  comm_allreduce(c, gs_scalar, gs_add, &norm, 1, wrk);
  normi = 1.0 / sqrt(norm);
  for (uint i = 0; i < n; i++) f[i] *= normi;

  free(vi);

  metric_toc(c, RSB_FIEDLER);

  return 0;
}
