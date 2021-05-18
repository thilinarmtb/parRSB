#include <limits.h>

#include <genmap-impl.h>
#include <genmap-sort.h>

static int check_convergence(struct comm *gc, int max_pass, int max_iter) {
  int max_levels = log2(gc->np);

  int i;
  for (i = 0; i < max_levels; i++) {
    sint converged = 1;
    int val = (int)metric_get_value(i, NFIEDLER);
    if (val >= max_pass * max_iter) {
      converged = 0;
      break;
    }

    sint ibfr;
    comm_allreduce(gc, gs_int, gs_min, &converged, 1, &ibfr);

    double dbfr;
    double final = (double)metric_get_value(i, LANCZOSTOLFINAL);
    comm_allreduce(gc, gs_double, gs_min, &final, 1, &dbfr);

    double target = (double)metric_get_value(i, LANCZOSTOLTARGET);

    if (converged == 0 && gc->id == 0) {
      printf("\tWarning: Partitioner only reached a tolerance of %lf given %lf "
             "after %d x %d iterations in Level=%d!\n",
             final, target, max_pass, max_iter, i);
      fflush(stdout);
    }
  }
}

int genmap_init_vector(genmap_vector *ivec_, int global, struct comm *lc,
                       genmap_handle h) {
  sint lelt = genmap_get_nel(h);
  genmap_vector_create(ivec_, lelt);
  genmap_vector ivec = *ivec_;

  slong out[2][1], buf[2][1];
  slong in = lelt;
  comm_scan(out, lc, gs_long, gs_add, &in, 1, buf);
  h->nel = out[1][0];

  struct rsb_element *elements = genmap_get_elements(h);
  sint i;
  if (global > 0)
    for (i = 0; i < lelt; i++)
      ivec->data[i] = out[0][0] + i + 1;
  else
    for (i = 0; i < lelt; i++)
      ivec->data[i] = elements[i].fiedler;

  genmap_vector_ortho_one(lc, ivec, genmap_get_partition_nel(h));
  GenmapScalar rtr = genmap_vector_dot(ivec, ivec);
  GenmapScalar rni;
  comm_allreduce(lc, gs_double, gs_add, &rtr, 1, &rni);
  rni = 1.0 / sqrt(rtr);
  genmap_vector_scale(ivec, ivec, rni);

  return 0;
}

int genmap_rsb(genmap_handle h) {
  int max_iter = 50;
  int max_pass = 50;

  struct comm *lc = h->local;
  struct comm *gc = h->global;

  genmap_number_faces_and_edges(h, gc);
  genmap_comm_scan(h, gc);

  uint nelt = genmap_get_nel(h);
  struct rsb_element *e = genmap_get_elements(h);

  int nv = h->nv;
  int ndim = (nv == 8) ? 3 : 2;

  int level = 0;

  while (lc->np > 1) {
    /* Run RCB, RIB pre-step or just sort by global id */
    metric_tic(lc, RCB);
    if (h->options->rsb_prepartition == 1) // RCB
      rcb(lc, h->elements, ndim, &h->buf);
    else if (h->options->rsb_prepartition == 2) // RIB
      rib(lc, h->elements, ndim, &h->buf);
    metric_toc(lc, RCB);

    /* Run fiedler */
    metric_tic(lc, FIEDLER);
    int ipass = 0, iter;
    do {
      genmap_vector ivec;
      genmap_init_vector(&ivec, ipass == 0, lc, h);

      iter = GenmapFiedler(h, lc, max_iter, ivec);
      metric_acc(NFIEDLER, iter);

      genmap_destroy_vector(ivec);
    } while (++ipass < max_pass && iter == max_iter);
    metric_toc(lc, FIEDLER);

    /* Sort by Fiedler vector */
    metric_tic(lc, FIEDLERSORT);
    parallel_sort(struct rsb_element, h->elements, fiedler, gs_double, 0, 1, lc,
                  &h->buf);
    metric_toc(lc, FIEDLERSORT);

    /* Bisect */
    int bin = 1;
    if (lc->id < (lc->np + 1) / 2)
      bin = 0;
    repair_partitions(h, bin, level, lc, gc);

    struct comm tc;
    genmap_comm_split(lc, bin, lc->id, &tc);
    comm_free(lc);
    comm_dup(lc, &tc);
    comm_free(&tc);

    genmap_comm_scan(h, lc);
    metric_push_level();
    level++;
  }

  check_convergence(gc, max_pass, max_iter);

  return 0;
}
