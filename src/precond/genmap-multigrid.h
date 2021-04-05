#ifndef _GENMAP_MULTIGRID_H_
#define _GENMAP_MULTIGRID_H_

#include <genmap-multigrid-csr.h>
#include <genmap-multigrid-gs.h>

struct mg_data;

struct mg_level {
  int (*get_nsmooth)(struct mg_data *data);
  GenmapScalar (*get_sigma)(struct mg_data *data);
  void (*G)(GenmapScalar *v, GenmapScalar *u);
  void (*restrict_down)(GenmapScalar *v, GenmapScalar *u);
  void (*interpolate_up)(GenmapScalar *v, GenmapScalar *u);
  void *data;
};

struct mg_data {
  int (*get_nlevels)(struct mg_data *data);
  uint *(*get_level_off)(struct mg_data *data);
  struct mg_level *(*get_levels)(struct mg_data *data);
  void (*setup)(genmap_handle h, struct comm *c, struct mg_data *d);
  void (*free)(struct mg_data *d);
  void *data;
};

int log2i(sint i);

void mg_setup(genmap_handle h, struct comm *c, struct mg_data *d);

void mg_vcycle(GenmapScalar *u, GenmapScalar *rhs, struct mg_data *d);

void mg_free(struct mg_data *d);

int project(genmap_handle h, struct comm *c, struct mg_data *d, genmap_vector r,
            int max_iter, genmap_vector x);

int flex_cg(genmap_handle h, struct comm *c, struct mg_data *d, genmap_vector r,
            int max_iter, genmap_vector x);

int rqi(genmap_handle h, struct comm *c, genmap_vector z, int max_iter,
        genmap_vector fiedler);

#endif
