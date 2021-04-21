#ifndef _GENMAP_MULTIGRID_H_
#define _GENMAP_MULTIGRID_H_

#include <genmap.h>

int log2i(sint i);

/* mg_data */
struct mg_data {
  int (*get_nlevels)(struct mg_data *data);
  uint *(*get_level_off)(struct mg_data *data);
  void (*free)(struct mg_data *d);

  int (*get_nsmooth)(struct mg_data *data, int lvl);
  GenmapScalar (*get_sigma)(struct mg_data *data, int lvl);
  GenmapScalar *(*get_diagonal)(struct mg_data *data, int lvl);
  void (*G)(struct mg_data *data, int lvl, GenmapScalar *v, GenmapScalar *u,
            buffer *buf);
  void (*rstrct)(struct mg_data *data, int lvl, GenmapScalar *v, buffer *buf);
  void (*intrp)(struct mg_data *data, int lvl, GenmapScalar *v, buffer *buf);

  void *data;
};

void mg_setup_gs(genmap_handle h, struct comm *c, struct mg_data *d);
void mg_setup_csr(genmap_handle h, struct comm *c, struct mg_data *d);
void mg_setup(genmap_handle h, struct comm *c, struct mg_data *d);

int mg_get_nlevels(struct mg_data *d);

uint *mg_get_level_off(struct mg_data *d);

int mg_get_nsmooth(struct mg_data *d, int level);

GenmapScalar mg_get_sigma(struct mg_data *d, int level);

GenmapScalar *mg_get_diagonal(struct mg_data *d, int level);

void mg_operator(GenmapScalar *v, GenmapScalar *u, int level, struct mg_data *d,
                 buffer *buf);

void mg_restrict(GenmapScalar *v, int level, struct mg_data *d, buffer *buf);

void mg_interpolate(GenmapScalar *v, int level, struct mg_data *d, buffer *buf);

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
