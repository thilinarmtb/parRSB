#include <limits.h>
#include <math.h>

#include <genmap-impl.h>
#include <genmap-iterative.h>
#include <genmap-partition.h>

static int fiedler_rqi(genmap_handle h, struct comm *gsc, int max_iter,
                       genmap_vector init_vec) {
  GenmapInt lelt = genmap_get_nel(h);
  genmap_vector y;
  genmap_vector_create_zeros(&y, lelt);

  metric_tic(gsc, RQI);
  int iter = rqi(h, gsc, init_vec, max_iter, y);
  metric_toc(gsc, RQI);
  metric_acc(NRQI, iter);

  GenmapScalar norm = 0, normi;
  sint i;
  for (i = 0; i < lelt; i++)
    norm += y->data[i] * y->data[i];
  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
  normi = 1.0 / sqrt(norm);

  genmap_vector_scale(y, y, normi);
  struct rsb_element *elements = genmap_get_elements(h);
  for (i = 0; i < lelt; i++)
    elements[i].fiedler = y->data[i];

  genmap_destroy_vector(y);

  return iter;
}

static int fiedler_lanczos(genmap_handle h, struct comm *gsc, int max_iter,
                           genmap_vector init_vec) {
  genmap_vector alphaVec, betaVec;
  genmap_vector_create(&alphaVec, max_iter);
  genmap_vector_create(&betaVec, max_iter - 1);
  genmap_vector *q = NULL;

  metric_tic(gsc, LANCZOS);
  int iter = genmap_lanczos(h, gsc, init_vec, max_iter, &q, alphaVec, betaVec);
  metric_toc(gsc, LANCZOS);
  metric_acc(NLANCZOS, iter);

  genmap_vector evLanczos, evTriDiag;
  genmap_vector_create(&evTriDiag, iter);

  /* Use TQLI and find the minimum eigenvalue and associated vector */
  genmap_vector *eVectors, eValues;
  metric_tic(gsc, TQLI);
  GenmapTQLI(h, alphaVec, betaVec, &eVectors, &eValues);
  metric_toc(gsc, TQLI);

  GenmapScalar eValMin = fabs(eValues->data[0]);
  GenmapInt eValMinI = 0;
  int i;
  for (i = 1; i < iter; i++) {
    if (fabs(eValues->data[i]) < eValMin) {
      eValMin = fabs(eValues->data[i]);
      eValMinI = i;
    }
  }
  genmap_vector_copy(evTriDiag, eVectors[eValMinI]);

  GenmapInt lelt = genmap_get_nel(h);
  genmap_vector_create_zeros(&evLanczos, lelt);
  GenmapInt j;
  for (i = 0; i < lelt; i++) {
    for (j = 0; j < iter; j++)
      evLanczos->data[i] += q[j]->data[i] * evTriDiag->data[j];
  }

  GenmapScalar norm = 0, normi;
  for (i = 0; i < lelt; i++)
    norm += evLanczos->data[i] * evLanczos->data[i];

  comm_allreduce(gsc, gs_double, gs_add, &norm, 1, &normi);
  genmap_vector_scale(evLanczos, evLanczos, 1. / sqrt(norm));
  struct rsb_element *elements = genmap_get_elements(h);
  for (i = 0; i < lelt; i++)
    elements[i].fiedler = evLanczos->data[i];

  genmap_destroy_vector(alphaVec);
  genmap_destroy_vector(betaVec);
  genmap_destroy_vector(evLanczos);
  genmap_destroy_vector(evTriDiag);
  genmap_destroy_vector(eValues);
  for (i = 0; i < iter; i++)
    genmap_destroy_vector(eVectors[i]);
  GenmapFree(eVectors);

  for (i = 0; i < iter + 1; i++)
    genmap_destroy_vector(q[i]);
  GenmapFree(q);

  return iter;
}

int GenmapFiedler(genmap_handle h, struct comm *lc, int max_iter,
                  genmap_vector init_vec) {
  genmap_laplacian_init(h, lc);

  if (h->options->rsb_algo == 0)
    return fiedler_lanczos(h, lc, max_iter, init_vec);
  else if (h->options->rsb_algo == 1)
    return fiedler_rqi(h, lc, max_iter, init_vec);

  if (h->M != NULL)
    csr_mat_free(h->M);
}
