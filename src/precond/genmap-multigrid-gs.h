#ifndef _GENMAP_MULTIGRID_GS_H_
#define _GENMAP_MULTIGRID_GS_H_

#include <genmap-multigrid.h>
#include <genmap.h>

struct mg_level_gs {
  int nsmooth;
  GenmapScalar sigma;
  struct gs_data *R;  // Restriction operator to level below
  struct gs_data *R0; // Restriction operartor from level 0
};

struct mg_data_gs {
  struct comm c;
  int nlevels;
  struct mg_level_gs *levels;
  uint *level_off;
};

void mg_setup_gs(genmap_handle h, struct comm *c, struct mg_data *d);

void mg_free_gs(struct mg_data *d);

#endif
