#include <genmap-ilu.h>
#include <genmap-multigrid.h>
#include <genmap-precond.h>

struct precond *precond_setup(int type, genmap_handle h, struct comm *c) {
  // Only does MG for now
  struct precond *p = tmalloc(struct precond, 1);
  p->type = type;

  if (type < 2) { /* MG (0 = CSR, 1 = GS) */
    p->data = (void *)tmalloc(struct mg_data, 1);
    mg_setup(h, c, type, (struct mg_data *)p->data);
  } else if (type == 2) {
    p->data = (void *)tmalloc(struct ilu_data, 1);
    ilu_setup(h, c, (struct ilu_data *)p->data);
  }

  return p;
}

void precond_apply(GenmapScalar *u, GenmapScalar *v, struct precond *p,
                   buffer *buf) {
  if (p->type < 2)
    mg_apply(u, v, (struct mg_data *)p->data, buf);
  else if (p->type == 2)
    ilu_apply(u, v, (struct mg_data *)p->data, buf);
}

void precond_check(genmap_handle h, struct comm *c) { mg_check(h, c); }

int precond_free(struct precond *p) {
  if (p != NULL) {
    if (p->type < 2)
      mg_free((struct mg_data *)p->data);
    else
      ilu_free((struct ilu_data *)p->data);
    free(p);
  }

  return 0;
}
