#include "parrsb_impl.h"

void arena_init(arena_t *arena, size_t size) {
  arena_t a = *arena = tcalloc(struct arena, 1);

  buffer_init(&(a->bfr), size);

  a->num_offset = 0;
  a->current_offset = 0;
  for (uint i = 0; i < MAXDEPTH; i++) a->previous_offset[i] = 0;
}

void arena_start(arena_t arena) {
  assert(arena->num_offset < MAXDEPTH && "MAXDEPTH is too small.");

  arena->previous_offset[arena->num_offset++] = arena->current_offset;
}

void *arena_alloc(arena_t arena, size_t size) {
  size_t offset = arena->current_offset;
  arena->current_offset += size;

  buffer *bfr = &(arena->bfr);
  buffer_reserve(bfr, arena->current_offset);

  return (void *)((char *)bfr->ptr + offset);
}

void arena_stop(arena_t arena) {
  if (arena->num_offset == 0) return;

  arena->current_offset = arena->previous_offset[--(arena->num_offset)];
}

void arena_free(arena_t *a) {
  if (!a || !(*a)) return;

  buffer_free(&((*a)->bfr));
  free(*a), *a = 0;
}
