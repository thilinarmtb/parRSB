#include "parrsb_impl.h"

static struct parrsb_options default_options = {
    // General options
    .partitioner = 0,
    .tagged = 0,
    .levels = 2,
    .find_disconnected_comps = 1,
    .repair = 0,
    .verbose_level = 1,
    .profile_level = 0,
    // RSB common (Lanczos and MG) options
    .rsb_algo = 0,
    .rsb_pre = 1,
    .rsb_max_iter = 50,
    .rsb_max_passes = 50,
    .rsb_tol = 1e-5,
    .rsb_dump_stats = 0,
    // RSB MG specific options
    .rsb_mg_grammian = 0,
    .rsb_mg_factor = 2};

int parrsb_options_get_default(parrsb_options *options) {
  *options = tcalloc(struct parrsb_options, 1);
  memcpy(*options, &default_options, sizeof(struct parrsb_options));
  return 0;
}

int parrsb_options_copy(parrsb_options *dest, const parrsb_options src) {
  *dest = tcalloc(struct parrsb_options, 1);
  memcpy(*dest, src, sizeof(struct parrsb_options));
  return 0;
}

void parrsb_options_print(const parrsb_options options) {
#define PRINT_OPTION(OPT, STR, FMT) printf("%s = " FMT "", STR, options->OPT)

  PRINT_OPTION(partitioner, "PARRSB_PARTITIONER", "%d");
  PRINT_OPTION(tagged, "PARRSB_TAGGED", "%d");
  PRINT_OPTION(levels, "PARRSB_LEVELS", "%d");
  PRINT_OPTION(find_disconnected_comps, "PARRSB_FIND_DISCONNECTED_COMPONENTS",
               "%d");
  PRINT_OPTION(repair, "PARRSB_REPAIR", "%d");
  PRINT_OPTION(verbose_level, "PARRSB_VERBOSE_LEVEL", "%d");
  PRINT_OPTION(profile_level, "PARRSB_PROFILE_LEVEL", "%d");
  PRINT_OPTION(rsb_algo, "PARRSB_RSB_ALGO", "%d");
  PRINT_OPTION(rsb_pre, "PARRSB_RSB_PRE", "%d");
  PRINT_OPTION(rsb_max_iter, "PARRSB_RSB_MAX_ITER", "%d");
  PRINT_OPTION(rsb_max_passes, "PARRSB_RSB_MAX_PASSES", "%d");
  PRINT_OPTION(rsb_tol, "PARRSB_RSB_TOL", "%lf");
  PRINT_OPTION(rsb_mg_grammian, "PARRSB_RSB_MG_GRAMMIAN", "%d");
  PRINT_OPTION(rsb_mg_factor, "PARRSB_RSB_MG_FACTOR", "%d");

#undef PRINT_OPTION
}

int parrsb_options_free(parrsb_options *options) {
  if (!options) return 1;
  if (*options) free(*options);
  *options = 0;
  return 0;
}
