#include "metrics.h"
#include "multigrid.h"
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
  comm_allreduce(c, gs_double, gs_add, &sum, 1, &wrk);
  sum /= n;

  for (uint i = 0; i < lelt; i++) q[i] -= sum;
}

static int lanczos_aux(scalar *diag, scalar *upper, scalar *rr, uint lelt,
                       ulong nelg, int niter, double tol, scalar *f,
                       laplacian gl, const struct comm *gsc, buffer *bfr) {
  scalar *r = tcalloc(scalar, 3 * lelt), *p = r + lelt, *w = p + lelt;
  // vec_copy(r, f);
  uint i;
  for (i = 0; i < lelt; i++) r[i] = f[i];

  // vec_ortho(gsc, r, nelg);
  ortho(r, lelt, nelg, gsc);

  scalar rtz1 = 1, pap = 0, alpha, beta, rtz2, pap_old;
  scalar rtr = dot(r, r, lelt), buf[2];
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
  scalar rnorm = sqrt(rtr), rtol = rnorm * tol;

  metric_set(TOL_INIT, rnorm);
  metric_set(TOL_TGT, rtol);

  // vec_scale(rr[0], r, rni);
  scalar rni = 1.0 / rnorm;
  for (i = 0; i < lelt; i++) rr[0 * lelt + i] = r[i] * rni;

  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1, rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0) beta = 0.0;

    // add2s1(p,r,beta,n)
    for (i = 0; i < lelt; i++) p[i] = beta * p[i] + r[i];

    scalar pp = dot(p, p, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &pp, 1, buf);

    // vec_ortho(gsc, p, nelg);
    ortho(p, lelt, nelg, gsc);

    laplacian_op(w, gl, p, bfr);

    pap_old = pap, pap = dot(w, p, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, buf);

    alpha = rtz1 / pap;
    // vec_axpby(r, r, 1.0, w, -1.0 * alpha);
    for (i = 0; i < lelt; i++) r[i] = r[i] - alpha * w[i];

    rtr = dot(r, r, lelt);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, buf);
    rnorm = sqrt(rtr), rni = 1.0 / rnorm;

    // vec_scale(rr[iter + 1], r, rni);
    for (i = 0; i < lelt; i++) rr[(iter + 1) * lelt + i] = r[i] * rni;

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

static int lanczos(scalar *fiedler, laplacian wl, scalar *initv,
                   const struct comm *gsc, const parrsb_options opts,
                   slong nelg, buffer *bfr) {
  metric_tic(gsc, RSB_LANCZOS_SETUP);

  uint miter = opts->rsb_max_iter;
  uint mpass = opts->rsb_max_passes;
  double tol = opts->rsb_tol;
  uint lelt = laplacian_get_size(wl);

  if (nelg < miter) miter = nelg;

  scalar *alpha = tcalloc(scalar, 2 * miter - 1);
  scalar *beta = alpha + miter;
  scalar *rr = tcalloc(scalar, (miter + 1) * lelt);

  scalar *eVectors = tcalloc(scalar, miter * miter);
  scalar *eValues = tcalloc(scalar, miter);

  metric_toc(gsc, RSB_LANCZOS_SETUP);

  uint iter = miter, ipass;
  for (ipass = 0; iter == miter && ipass < mpass; ipass++) {
    metric_tic(gsc, RSB_LANCZOS);
    iter = lanczos_aux(alpha, beta, rr, lelt, nelg, miter, tol, initv, wl, gsc,
                       bfr);
    metric_toc(gsc, RSB_LANCZOS);

    metric_tic(gsc, RSB_LANCZOS_TQLI);

    // Use TQLI and find the minimum eigenvalue and associated vector
    tqli(eVectors, eValues, iter, alpha, beta, gsc->id);
    scalar eValMin = fabs(eValues[0]);
    uint eValMinI = 0;
    for (uint i = 1; i < iter; i++) {
      if (fabs(eValues[i]) < eValMin) {
        eValMin = fabs(eValues[i]);
        eValMinI = i;
      }
    }

    for (uint i = 0; i < lelt; i++) {
      fiedler[i] = 0.0;
      for (uint j = 0; j < iter; j++)
        fiedler[i] += rr[j * lelt + i] * eVectors[eValMinI * iter + j];
    }
    ortho(fiedler, lelt, nelg, gsc);
    for (uint i = 0; i < lelt; i++) initv[i] = fiedler[i];

    metric_toc(gsc, RSB_LANCZOS_TQLI);
  }

  free(alpha), free(rr), free(eVectors), free(eValues);

  return (ipass - 1) * miter + iter;
}

static int fiedler_mesh(struct array *elems, const element_info ei,
                        const parrsb_options opts, const struct comm *gsc,
                        buffer *buf) {
  metric_tic(gsc, RSB_FIEDLER);

  // Return if the number of processes is equal to 1.
  if (gsc->np == 1) goto early_exit;

  metric_tic(gsc, RSB_FIEDLER_SETUP);

  uint lelt = elems->n;
  slong out[2][1], wrk[2][1], in = lelt;
  comm_scan(out, gsc, gs_long, gs_add, &in, 1, wrk);
  slong start = out[0][0], nelg = out[1][0];

  scalar *initv = tcalloc(scalar, lelt);
  for (uint i = 0; i < lelt; i++) {
    initv[i] = start + i + 1.0;
    if (start + i < nelg / 2) initv[i] += 1000 * nelg;
  }

  ortho(initv, lelt, nelg, gsc);
  scalar rtr = dot(initv, initv, lelt), rni;
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);
  rni = 1.0 / sqrt(rtr);
  for (uint i = 0; i < lelt; i++) initv[i] *= rni;

  metric_toc(gsc, RSB_FIEDLER_SETUP);

  metric_tic(gsc, RSB_LAPLACIAN_SETUP);
  laplacian wl;
  laplacian_init(&wl, elems, ei, gsc, buf);
  metric_toc(gsc, RSB_LAPLACIAN_SETUP);

  metric_tic(gsc, RSB_FIEDLER_CALC);

  int iter = 0;
  scalar *f = tcalloc(scalar, lelt);
  switch (opts->rsb_algo) {
  case 0: iter = lanczos(f, wl, initv, gsc, opts, nelg, buf); break;
  case 1:
    iter = inverse(f, wl, elems, ei->nv, initv, gsc, opts, nelg, buf);
    break;
  }
  metric_acc(RSB_FIEDLER_CALC_NITER, iter);

  scalar norm = 0;
  for (uint i = 0; i < lelt; i++) norm += f[i] * f[i];
  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, wrk);
  scalar normi = 1.0 / sqrt(norm);
  for (uint i = 0; i < lelt; i++) f[i] *= normi;

  metric_toc(gsc, RSB_FIEDLER_CALC);

  laplacian_free(&wl);

  struct rsb_element *pe = (struct rsb_element *)elems->ptr;
  for (uint i = 0; i < lelt; i++) pe[i].fiedler = f[i];

  if (initv) free(initv);
  if (f) free(f);

early_exit:
  metric_toc(gsc, RSB_FIEDLER);
  return 0;
}

int fiedler(scalar *f, laplacian l, const parrsb_options opts,
            const struct comm *c, buffer *buf) {
  // Return if the number of processes is equal to 1 or rsb algorithm is set to
  // something other than lanczos.
  if (c->np == 1 || opts->rsb_algo > 0) return 1;

  metric_tic(c, RSB_FIEDLER);

  uint n = laplacian_get_size(l);
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
  comm_allreduce(c, gs_double, gs_add, &norm, 1, wrk);
  scalar normi = 1.0 / sqrt(norm);
  for (uint i = 0; i < n; i++) vi[i] *= normi;

  int iter = 0;
  switch (opts->rsb_algo) {
  case 0: iter = lanczos(f, l, vi, c, opts, ng, buf); break;
  default: break;
  }
  metric_acc(RSB_FIEDLER_CALC_NITER, iter);

  norm = dot(f, f, n);
  comm_allreduce(c, gs_double, gs_add, &norm, 1, wrk);
  normi = 1.0 / sqrt(norm);
  for (uint i = 0; i < n; i++) f[i] *= normi;

  free(vi);

  metric_toc(c, RSB_FIEDLER);

  return 0;
}
