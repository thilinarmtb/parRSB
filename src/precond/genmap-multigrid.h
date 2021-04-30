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
  void (*diagonal_scaling)(struct mg_data *data, int lvl, GenmapScalar *v,
                           GenmapScalar *u, GenmapScalar sigma);
  void (*G)(struct mg_data *data, int lvl, GenmapScalar *v, GenmapScalar *u,
            buffer *buf);
  void (*rstrct)(struct mg_data *data, int lvl, GenmapScalar *v, buffer *buf);
  void (*intrp)(struct mg_data *data, int lvl, GenmapScalar *v, buffer *buf);
  void (*coarse)(struct mg_data *data, GenmapScalar *u, GenmapScalar *r);

  genmap_handle h;
  void *data;
};

void mg_check(genmap_handle h, struct comm *c);

void mg_setup_gs(genmap_handle h, struct comm *c, struct mg_data *d);
void mg_setup_csr(genmap_handle h, struct comm *c, struct mg_data *d);
void mg_setup(genmap_handle h, struct comm *c, int type, struct mg_data *d);

int mg_get_nlevels(struct mg_data *d);

uint *mg_get_level_off(struct mg_data *d);

int mg_get_nsmooth(struct mg_data *d, int level);

GenmapScalar mg_get_sigma(struct mg_data *d, int level);

void mg_diagonal_scaling(GenmapScalar *v, GenmapScalar *u, GenmapScalar sigma,
                         int level, struct mg_data *d);

void mg_operator(GenmapScalar *v, GenmapScalar *u, int level, struct mg_data *d,
                 buffer *buf);

void mg_restrict(GenmapScalar *v, int level, struct mg_data *d, buffer *buf);

void mg_interpolate(GenmapScalar *v, int level, struct mg_data *d, buffer *buf);

void mg_coarse_solve(GenmapScalar *u, GenmapScalar *r, struct mg_data *d);

void mg_free(struct mg_data *d);

// TODO: Move this to a separate file
/* Algorithms */
void vcycle(GenmapScalar *u, GenmapScalar *rhs, struct mg_data *d, buffer *buf);

int project(genmap_handle h, struct comm *c, struct mg_data *d, genmap_vector r,
            int max_iter, genmap_vector x);

int flex_cg(genmap_handle h, struct comm *c, struct mg_data *d, genmap_vector r,
            int max_iter, genmap_vector x);

int rqi(genmap_handle h, struct comm *c, genmap_vector z, int max_iter,
        genmap_vector fiedler);

#endif
