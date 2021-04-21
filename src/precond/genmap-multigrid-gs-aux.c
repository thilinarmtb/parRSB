// FIXME: genmap-impl.h is only for GenmapFree
#include <genmap-impl.h>
#include <genmap-multigrid-gs.h>

void mg_setup_aux_gs(struct mg_data_gs *d, slong *wrk, buffer *buf) {
  struct comm *c = &d->c;

  uint nelt = d->level_off[1] - d->level_off[0];

  slong in;
  slong out[2][1], bf[2][1];

  uint size0, size1;
  size0 = nelt;

  uint j, lvl;
  for (lvl = 1; lvl < d->nlevels; lvl++) {
    in = size0;
    comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
    slong start0 = out[0][0];
    slong ng0 = out[1][0];

    size1 = (start0 + size0 + 1) / 2 - (start0 + 1 + 1) / 2;

    in = size1;
    comm_scan(out, c, gs_long, gs_add, &in, 1, buf);
    slong start1 = out[0][0];
    slong ng1 = out[1][0];

    for (j = 0; j < size0; j++)
      wrk[j] = -((start0 + j + 1) + 1) / 2;
    if (ng0 > 2 && (ng0 % 2) == 1 && (start0 + size0) == ng0)
      wrk[j - 1] += 1;
    for (j = 0; j < size1; j++)
      wrk[size0 + j] = start1 + j + 1;

    d->J[lvl - 1] = gs_setup(wrk, size0 + size1, c, 0, gs_pairwise, 0);

    d->nsmooth[lvl] = 2;
    d->sigma[lvl] = 0.6;
    d->J[lvl] = NULL;

    d->level_off[lvl + 1] = d->level_off[lvl] + size1;

    size0 = size1;
  }

  for (lvl = d->nlevels - 2; lvl > 0; lvl--) {
    size1 = d->level_off[lvl + 1] - d->level_off[lvl];

    in = size1;
    comm_scan(out, &d->c, gs_long, gs_add, &in, 1, buf);
    slong eid = out[0][0];

    uint cur_off = d->level_off[lvl];

    for (j = 0; j < size1; j++)
      wrk[cur_off + j] = eid + j + 1;
    // Apply R^T till we get to Level0
    for (j = lvl - 1; j >= 0; j--)
      gs(wrk + d->level_off[j], gs_long, gs_add, 0, d->J[j], buf);

    // Setup R0 now that we know the corresponding ids at level0
    for (j = 0; j < nelt; j++)
      wrk[j] *= -1;
    for (j = 0; j < size1; j++)
      wrk[nelt + j] = wrk[cur_off + j];

    // TODO: See if I can call gs_setup without c and get rid of eid
    d->R0[lvl] = gs_setup(wrk, nelt + size1, c, 0, gs_crystal_router, 0);
  }
}
