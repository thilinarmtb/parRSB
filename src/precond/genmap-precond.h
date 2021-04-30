#ifndef _GENMAP_PRECOND_H_
#define _GENMAP_PRECOND_H_

#include <genmap.h>

struct precond {
  void (*apply)(GenmapScalar *u, GenmapScalar *v, void *data, buffer *buf);
  void *data;
};

struct precond *precond_setup(int type, genmap_handle h, struct comm *c);

void precond(GenmapScalar *u, GenmapScalar *v, struct precond *p, buffer *buf);

int precond_free(struct precond *p);

#endif
