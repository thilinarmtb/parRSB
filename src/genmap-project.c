#include <genmap-iterative.h>

#define MM 505

int project(genmap_handle h, struct comm *gsc, struct precond *d,
            genmap_vector ri, int max_iter, genmap_vector x) {
  assert(x->size == ri->size);
  uint lelt = x->size;

  slong bfr, nelg;
  nelg = lelt;
  comm_allreduce(gsc, gs_long, gs_add, &nelg, 1, &bfr);

  genmap_vector z0, z, dz, w, p, r;
  genmap_vector_create(&z, lelt);
  genmap_vector_create(&w, lelt);
  genmap_vector_create(&r, lelt);
  genmap_vector_create(&p, lelt);
  genmap_vector_create(&z0, lelt);
  genmap_vector_create(&dz, lelt);

  assert(max_iter < MM);
  double *P, *W;
  GenmapCalloc(lelt * MM, &P);
  GenmapCalloc(lelt * MM, &W);

  uint i;
  for (i = 0; i < lelt; i++) {
    x->data[i] = 0.0;
    r->data[i] = ri->data[i];
  }

  genmap_vector_copy(z, r);
  genmap_vector_copy(p, z);

  GenmapScalar rz1, buf;
  rz1 = genmap_vector_dot(r, z);
  comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, &buf);

  GenmapScalar rr = genmap_vector_dot(r, r);
  comm_allreduce(gsc, gs_double, gs_add, &rr, 1, &buf);

  GenmapScalar alpha, beta, rz0, rz2, scale;

  double tol = 1e-3;
  double res_tol = rr * tol;

  uint j, k;
  i = 0;
  while (i < max_iter) {
    genmap_laplacian(h, p->data, w->data);

    GenmapScalar den = genmap_vector_dot(p, w);
    comm_allreduce(gsc, gs_double, gs_add, &den, 1, &buf);
    alpha = rz1 / den;

    scale = 1.0 / sqrt(den);
    for (j = 0; j < lelt; j++) {
      W[i * lelt + j] = scale * w->data[j];
      P[i * lelt + j] = scale * p->data[j];
    }

    genmap_vector_axpby(x, x, 1.0, p, alpha);
    genmap_vector_axpby(r, r, 1.0, w, -alpha);

    rr = genmap_vector_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rr, 1, &buf);

    if (gsc->id == 0)
      printf("i = %u, rr = %lf, res_tol = %lf, tol = %lf\n", i, rr, res_tol,
             tol * tol);
    if (rr < res_tol || rr < tol * tol)
      break;

    GenmapScalar norm0 = genmap_vector_dot(z, z);
    comm_allreduce(gsc, gs_double, gs_add, &norm0, 1, &buf);

    genmap_vector_copy(z0, z);

    metric_tic(gsc, VCYCLE);
    precond_apply(z->data, r->data, d, &h->buf);
    metric_toc(gsc, VCYCLE);

    GenmapScalar norm1 = genmap_vector_dot(z, z);
    comm_allreduce(gsc, gs_double, gs_add, &norm1, 1, &buf);

    rz0 = rz1;
    genmap_vector_ortho_one(gsc, z, nelg);
    rz1 = genmap_vector_dot(r, z);
    comm_allreduce(gsc, gs_double, gs_add, &rz1, 1, &buf);

    genmap_vector_axpby(dz, z, 1.0, z0, -1.0);
    rz2 = genmap_vector_dot(r, dz);
    comm_allreduce(gsc, gs_double, gs_add, &rz2, 1, &buf);

    beta = rz2 / rz0;
    genmap_vector_axpby(p, z, 1.0, p, beta);

    i++;

    metric_tic(gsc, PROJECT);
    for (k = 0; k < lelt; k++)
      P[(MM - 1) * lelt + k] = 0.0;

    for (j = 0; j < i; j++) {
      double a = 0.0;
      for (k = 0; k < lelt; k++)
        a += W[j * lelt + k] * p->data[k];
      comm_allreduce(gsc, gs_double, gs_add, &a, 1, &buf);
      for (k = 0; k < lelt; k++)
        P[(MM - 1) * lelt + k] += a * P[j * lelt + k];
    }

    for (k = 0; k < lelt; k++)
      p->data[k] -= P[(MM - 1) * lelt + k];
    metric_toc(gsc, PROJECT);
  }

  genmap_destroy_vector(z);
  genmap_destroy_vector(w);
  genmap_destroy_vector(p);
  genmap_destroy_vector(r);
  genmap_destroy_vector(z0);
  genmap_destroy_vector(dz);

  GenmapFree(P);
  GenmapFree(W);

  return i + 1;
}
