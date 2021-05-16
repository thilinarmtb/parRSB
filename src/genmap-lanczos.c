#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>

int genmap_lanczos(genmap_handle h, struct comm *gsc, genmap_vector f,
                   GenmapInt niter, genmap_vector **rr, genmap_vector diag,
                   genmap_vector upper) {
  assert(diag->size == niter);
  assert(diag->size == upper->size + 1);
  assert(f->size == genmap_get_nel(h));

  if (genmap_get_partition_nel(h) < niter) {
    niter = genmap_get_partition_nel(h);
    diag->size = niter;
    upper->size = niter - 1;
  }

  GenmapScalar eps = 1.e-5;
  GenmapScalar alpha, beta;
  GenmapScalar rnorm, rtol, rni, rtr, rtz1, rtz2, pap, pap_old;
  genmap_vector r, p, w;

  rtz1 = 1.0;
  pap = 0.0;
  GenmapInt lelt = genmap_get_nel(h);

  genmap_vector_create_zeros(&p, lelt);
  genmap_vector_create(&w, lelt);
  genmap_vector_create(&r, lelt);

  genmap_vector_copy(r, f);

  genmap_vector_ortho_one(gsc, r, genmap_get_partition_nel(h));
  rtr = genmap_vector_dot(r, r);
  comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);
  rnorm = sqrt(rtr);
  rtol = rnorm * eps;
  rni = 1.0 / rnorm;

  if (*rr == NULL) {
    GenmapMalloc((size_t)(niter + 1), rr);
    GenmapInt i;
    for (i = 0; i < niter + 1; ++i)
      (*rr)[i] = NULL;
  }
  genmap_vector_create(&(*rr)[0], lelt);

  genmap_vector_scale((*rr)[0], r, rni);

  int iter;
  for (iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1 / rtz2;
    if (iter == 0)
      beta = 0.0;

    genmap_vector_axpby(p, p, beta, r, 1.0);
    genmap_vector_ortho_one(gsc, p, genmap_get_partition_nel(h));

    genmap_laplacian(h, p->data, w->data);

    genmap_vector_scale(w, w, -1.0);

    pap_old = pap;
    pap = genmap_vector_dot(w, p);
    comm_allreduce(gsc, gs_double, gs_add, &pap, 1, &rni);

    alpha = rtz1 / pap;
    genmap_vector_axpby(r, r, 1.0, w, -1.0 * alpha);

    rtr = genmap_vector_dot(r, r);
    comm_allreduce(gsc, gs_double, gs_add, &rtr, 1, &rni);
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;

    genmap_vector_create(&(*rr)[iter + 1], lelt);
    genmap_vector_scale((*rr)[iter + 1], r, rni);

    if (iter == 0) {
      diag->data[iter] = pap / rtz1;
    } else {
      diag->data[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper->data[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }

    if (rnorm < rtol) {
      diag->size = iter + 1;
      upper->size = iter;
      iter = iter + 1;
      break;
    }
  }

  metric_acc(LANCZOSTOLFINAL, rnorm);
  metric_acc(LANCZOSTOLTARGET, rtol);

  genmap_destroy_vector(p);
  genmap_destroy_vector(w);
  genmap_destroy_vector(r);

  return iter;
}
