#ifndef _GENMAP_MULTIGRID_H_
#define _GENMAP_MULTIGRID_H_

#include <genmap.h>

/* mg_level */
struct mg_level {
  int (*get_nsmooth)(struct mg_level *lvl);
  GenmapScalar (*get_sigma)(struct mg_level *lvl);
  void (*G)(struct mg_level *lvl, GenmapScalar *v, GenmapScalar *u);
  void (*rstrct)(struct mg_level *lvl, GenmapScalar *v, GenmapScalar *u);
  void (*interp)(struct mg_level *lvl, GenmapScalar *v, GenmapScalar *u);
  void *data;
};

int mg_level_get_nsmooth(struct mg_level *lvl);

GenmapScalar mg_level_get_sigma(struct mg_level *lvl);

/* mg_data */
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

int mg_get_nlevels(struct mg_data *d);

uint *mg_get_level_off(struct mg_data *d);

void mg_free(struct mg_data *d);

/* Algorithms */
void vcycle(GenmapScalar *u, GenmapScalar *rhs, struct mg_data *d, buffer *buf);

int project(genmap_handle h, struct comm *c, struct mg_data *d, genmap_vector r,
            int max_iter, genmap_vector x);

int flex_cg(genmap_handle h, struct comm *c, struct mg_data *d, genmap_vector r,
            int max_iter, genmap_vector x);

int rqi(genmap_handle h, struct comm *c, genmap_vector z, int max_iter,
        genmap_vector fiedler);

#endif
