#include "parrsb_impl.h"

static void buffer_initialize_fields(buffer *bfr) {
  bfr->n = bfr->max = 0, bfr->ptr = 0;
}

static int buffer_initialized(buffer *bfr) {
  return !(bfr->n == 0 && bfr->ptr == 0 && bfr->max == 0);
}

void arena_init(arena_t *arena) {
  arena_t a = *arena = tcalloc(struct arena, 1);
  a->n = 0;
  for (uint i = 0; i < MAXDEPTH; i++) buffer_initialize_fields(&(a->bfr[i]));
}

void *arena_start(arena_t a, size_t capacity) {
  assert(a->n < MAXDEPTH && "MAXDEPTH is too small!");

  buffer *bfr = &(a->bfr[a->n++]);
  if (!buffer_initialized(bfr))
    buffer_init(bfr, capacity);
  else
    buffer_reserve(bfr, capacity);

  return bfr->ptr;
}

void *arena_alloc(arena_t a, size_t size) {
  assert(a->n > 0 && "arena_start() must be called before arena_talloc()!");

  buffer *bfr = &(a->bfr[a->n - 1]);
  assert((bfr->n + size) <= bfr->max && "arena capacity is not enough!");
  size_t offset = bfr->n;
  bfr->n += size;

  return (void *)((char *)bfr->ptr + offset);
}

void arena_stop(arena_t a) {
  assert(a->n > 0 && "arena_stop() must be called after arena_start()!");
  a->bfr[--(a->n)].n = 0;
}

void arena_free(arena_t *arena) {
  if (!arena || !(*arena)) return;

  arena_t a = *arena;
  for (uint i = 0; i < MAXDEPTH; i++) {
    buffer *bfr = &(a->bfr[i]);
    if (buffer_initialized(bfr)) buffer_free(bfr);
  }

  free(a), a = 0;
}
