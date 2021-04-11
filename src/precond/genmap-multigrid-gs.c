#include <genmap-impl.h>
#include <genmap-multigrid-gs.h>

/* Levels internal to processor */
static int setup_internal_levels(genmap_handle h, struct comm *c,
                                 struct mg_data_gs *d) {
  sint nelt = genmap_get_nel(h);
  int ninternal = log2i(nelt) + 1;
  // Global reduction to find maximum internal levels
  slong buf[2][2];
  comm_allreduce(c, gs_int, gs_max, &ninternal, 1, buf);

  int nexternal = log2i(c->np) + 1;
  d->nlevels = ninternal + nexternal;

  // Allocate data structures
  GenmapMalloc(d->nlevels + 1, &d->level_off);
  GenmapMalloc(d->nlevels, &d->levels);

  sint off_size = nelt;
  d->level_off[0] = 0;
  int i;
  for (i = 1; i < ninternal + 1; i++) {
    d->level_off[i] = d->level_off[i - 1] + off_size;
    off_size /= 2;
  }

  sint k = 1;
  for (i = ninternal + 1; i < d->nlevels + 1; i++) {
    if ((c->id + 1) % k == 0)
      d->level_off[i] = d->level_off[i - 1] + 1;
    else
      d->level_off[i] = d->level_off[i - 1];
    k *= 2;
  }

  slong *wrk;
  size_t wrk_size = 2 * nelt + c->np;
  GenmapMalloc(wrk_size, &wrk);

  slong in[2];
  slong out[2][2];
  uint j, cur_size, next_size;
  int nlevels = d->nlevels;

  // Setup R, last internal level is the first external level.
  // It would be taken care at external R setup
  for (i = ninternal - 2; i >= 0; i--) {
    // Setup global ids for gs_setup
    cur_size = d->level_off[i + 1] - d->level_off[i];
    next_size = d->level_off[i + 2] - d->level_off[i + 1];

    in[0] = next_size;
    comm_scan(out, &d->c, gs_long, gs_add, in, 1, buf);
    slong eid = out[0][0];

    for (j = 0; j < cur_size; j++)
      wrk[j] = -(eid + j / 2 + 1);
    if (j > 2 && j % 2 == 0)
      wrk[j] += 1;

    for (j = 0; j < next_size; j++)
      wrk[cur_size + j] = eid + j + 1;

    // TODO: See if I can call gs_setup without c and get rid of eid
    d->levels[i].R =
        gs_setup(wrk, next_size + cur_size, c, 0, gs_crystal_router, 0);
  }

  // External R setup. Last external level doesn't need R
  for (i = nlevels - 2; i >= ninternal - 1; i--) {
    // Both next_size and cur_size are equal to 1 or  0
    cur_size = d->level_off[i + 1] - d->level_off[i];
    next_size = d->level_off[i + 2] - d->level_off[i + 1];

    in[0] = cur_size;
    in[1] = next_size;
    comm_scan(out, &d->c, gs_long, gs_add, in, 2, buf);
    slong cur_id = out[0][0];
    slong next_id = out[0][1];

    if (cur_size == 1) {
      wrk[0] = -(cur_id / 2 + 1);
      if (cur_id + 1 == out[1][0] && cur_id % 2 == 0)
        wrk[0] += 1;
    }

    if (next_size == 1) // next_size == 1 iff cur_size == 1
      wrk[1] = next_id + 1;

    d->levels[i].R =
        gs_setup(wrk, next_size + cur_size, c, 0, gs_crystal_router, 0);
  }

  // Setup R0, last level and zeroth level doesn't need R0
  for (i = nlevels - 2; i > 0; i--) {
    cur_size = d->level_off[i + 1] - d->level_off[i];

    in[0] = cur_size;
    comm_scan(out, &d->c, gs_long, gs_add, in, 1, buf);
    slong eid = out[0][0];

    uint cur_off = d->level_off[i];

    for (j = 0; j < cur_size; j++)
      wrk[cur_off + j] = eid + j;
    // Apply R^T till we get to Level0
    for (j = i - 1; j >= 0; j--)
      gs(wrk + d->level_off[j], gs_double, gs_add, 0, d->levels[j].R, &h->buf);

    // Setup R0 now that we know the corresponding ids at level0
    for (j = 0; j < nelt; j++)
      wrk[j] *= -1;
    for (j = 0; j < cur_size; j++)
      wrk[nelt + j] = wrk[cur_off + j];

    // TODO: See if I can call gs_setup without c and get rid of eid
    d->levels[i].R0 =
        gs_setup(wrk, nelt + cur_size, c, 0, gs_crystal_router, 0);
  }

  GenmapFree(wrk);

  return nlevels;
}

void mg_setup_gs(genmap_handle h, struct comm *c, struct mg_data *dd) {
  struct mg_data_gs *d = dd->data;

  comm_dup(&d->c, c);

  slong in = genmap_get_nel(h);
  slong out[2][1], bf[2][1];
  comm_scan(out, &d->c, gs_long, gs_add, &in, 1, bf);
  slong nelg = out[1][0];

  d->nlevels = log2i(nelg) + 1;
  GenmapMalloc(d->nlevels, &d->levels);
  GenmapMalloc(d->nlevels + 1, &d->level_off);
}

void mg_free_gs(struct mg_data *dd) {
  struct mg_data_gs *d = dd->data;

  int i;
  struct mg_level_gs *l = d->levels;
  for (i = 0; i < d->nlevels; i++) {
    if (i < d->nlevels - 1)
      gs_free(l[i].R);
    GenmapFree(&l[i]);
    if (i > 0 && i < d->nlevels - 1)
      gs_free(l[i].R0);
  }
  GenmapFree(l);

  comm_free(&d->c);
  GenmapFree(d->level_off);
  GenmapFree(d);
}
