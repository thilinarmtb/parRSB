#include <genmap-multigrid.h>
#include <genmap-precond.h>

struct precond *precond_setup(int type, genmap_handle h, struct comm *c) {
  // Only does MG for now
  struct precond *p = tmalloc(struct precond, 1);
  p->data = (void *)tmalloc(struct mg_data, 1);
  mg_setup(h, c, 0, (struct mg_data *)p->data);
  return p;
}

int precond_free(struct precond *p) {
  // Only does MG for now
  // TODO: mg_free();
  if (p->data != NULL)
    free(p->data);
  free(p);

  return 0;
}
