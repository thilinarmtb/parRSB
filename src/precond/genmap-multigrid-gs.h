#ifndef _GENMAP_MULTIGRID_GS_H_
#define _GENMAP_MULTIGRID_GS_H_

#include <genmap-multigrid.h>

struct mg_data_gs {
  struct comm *c;

  int nlevels;
  uint *level_off;

  int *nsmooth;
  GenmapScalar *sigma;

  struct gs_data **J;

  GenmapScalar *buf;

  struct gs_data **J0;
  GenmapScalar *diagonal; // Level 0 diagonal

  int nv;
};

void mg_setup_aux_gs(struct mg_data_gs *d, slong *wrk, buffer *buf);

#endif
