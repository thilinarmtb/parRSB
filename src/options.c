#include "parrsb_impl.h"

static struct parrsb_options default_options = {
    .partitioner = 0,
    .levels = 2,
    .verbose_level = 1,
    .profile_level = 0,
    .find_cc = 1,
    .repair = 0,
    .rsb_pre = 1,
    .rsb_max_iter = 50,
    .rsb_max_passes = 50,
    .rsb_tol = 1e-5,
};

int parrsb_options_get_default(parrsb_options_t *options) {
  *options = tcalloc(struct parrsb_options, 1);
  memcpy(*options, &default_options, sizeof(struct parrsb_options));
  return 0;
}

void parrsb_options_print(const parrsb_options_t options) {
#define PRINT_OPTION(OPT, STR, FMT) printf("%s = " FMT "\n", STR, options->OPT)

  PRINT_OPTION(partitioner, "PARRSB_PARTITIONER", "%u");
  PRINT_OPTION(levels, "PARRSB_LEVELS", "%u");
  PRINT_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL", "%d");
  PRINT_OPTION(profile_level, "PARRSB_PROFILE_LEVEL", "%d");
  PRINT_OPTION(find_cc, "PARRSB_FIND_CONNECTED_COMPONENTS", "%d");
  PRINT_OPTION(repair, "PARRSB_REPAIR", "%d");
  PRINT_OPTION(rsb_pre, "PARRSB_RSB_PRE", "%d");
  PRINT_OPTION(rsb_max_iter, "PARRSB_RSB_MAX_ITER", "%d");
  PRINT_OPTION(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", "%d");
  PRINT_OPTION(rsb_tol, "PARRSB_RSB_TOL", "%lf");

#undef PRINT_OPTION
}

int parrsb_options_update(parrsb_options_t options) {
#define OPT_UPDATE(opt, var, is_int)                                           \
  do {                                                                         \
    const char *val = getenv(var);                                             \
    if (val != NULL) {                                                         \
      if (is_int)                                                              \
        options->opt = atoi(val);                                              \
      else                                                                     \
        options->opt = atof(val);                                              \
    }                                                                          \
  } while (0)

  OPT_UPDATE(partitioner, "PARRSB_PARTITIONER", 1);
  OPT_UPDATE(levels, "PARRSB_LEVELS", 1);
  OPT_UPDATE(verbose_level, "PARRSB_VERBOSE_LEVEL", 1);
  OPT_UPDATE(profile_level, "PARRSB_PROFILE_LEVEL", 1);
  OPT_UPDATE(find_cc, "PARRSB_FIND_CONNECTED_COMPONENTS", 1);
  OPT_UPDATE(repair, "PARRSB_REPAIR", 1);
  OPT_UPDATE(rsb_pre, "PARRSB_RSB_PRE", 1);
  OPT_UPDATE(rsb_max_iter, "PARRSB_RSB_MAX_ITER", 1);
  OPT_UPDATE(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", 1);
  OPT_UPDATE(rsb_tol, "PARRSB_RSB_TOL", 0);

#undef OPT_UPDATE

  return 0;
}

int parrsb_options_set_partitioner(parrsb_options_t options, int partitioner) {
  if (partitioner < 0 || partitioner > 2) return 1;
  options->partitioner = partitioner;
  return 0;
}

int parrsb_options_free(parrsb_options_t *options) {
  if (!options) return 1;
  if (*options) free(*options), *options = 0;
  return 0;
}
