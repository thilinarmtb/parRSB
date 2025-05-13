#ifndef _PARRSB_EXAMPLE_H_
#define _PARRSB_EXAMPLE_H_

#include <parrsb.h>

#include <getopt.h>
#include <math.h>

static void parrsb_print_part_stat(long long *vtx, unsigned nelt, unsigned nv,
                                   MPI_Comm ce) {
  int id, np;
  MPI_Comm_rank(ce, &id);
  MPI_Comm_size(ce, &np);

  int nc[3], ns[3], nss[3], nel[3];
  parrsb_get_part_stat(&nc[0], &ns[0], &nss[0], &nel[0], vtx, nelt, nv, ce);

  if (id == 0) goto flush;
  printf("Min neighbors: %d | Max neighbors: %d | Avg neighbors: %lf\n", nc[0],
         nc[1], (double)nc[2] / np);
  printf("Min nvolume: %d | Max nvolume: %d | Avg nvolume: %lf\n", ns[0], ns[1],
         (double)ns[2] / np);
  printf("Min volume: %d | Max volume: %d | Avg volume: %lf\n", nss[0], nss[1],
         (double)nss[2] / np);
  printf("Min elements: %d | Max elements: %d\n", nel[0], nel[1]);
flush:
  fflush(stdout);
}

typedef struct {
  char *mesh;  // Mesh name, required.
  double tol;  // gencon tolerance, default: 0.2
  int test;    // run tests, default: 0
  int dump;    // dump the connectivity or map file, default: 1
  int nactive; // # of active MPI ranks, default: INT_MAX
  int verbose; // Verbosity, default: 0
} parrsb_cmd_line_opts;

static parrsb_cmd_line_opts *parrsb_parse_cmd_opts(int argc, char *argv[]) {
  parrsb_cmd_line_opts *in = tcalloc(parrsb_cmd_line_opts, 1);

  in->mesh = NULL, in->tol = 2e-1;
  in->test = 0, in->dump = 0, in->verbose = 0, in->nactive = INT_MAX;

  static struct option long_options[] = {{"mesh", required_argument, 0, 0},
                                         {"tol", optional_argument, 0, 10},
                                         {"test", optional_argument, 0, 20},
                                         {"dump", optional_argument, 0, 30},
                                         {"nactive", optional_argument, 0, 40},
                                         {"verbose", optional_argument, 0, 50},
                                         {"help", optional_argument, 0, 99},
                                         {0, 0, 0, 0}};

  size_t len;
  for (;;) {
    int c = getopt_long(argc, argv, "", long_options, NULL);
    if (c == -1) break;

    switch (c) {
    case 0:
      len = strnlen(optarg, PATH_MAX);
      in->mesh = tcalloc(char, len + 1);
      strncpy(in->mesh, optarg, len);
      break;
    case 10: in->tol = atof(optarg); break;
    case 20: in->test = 1; break;
    case 30: in->dump = 1; break;
    case 40: in->nactive = atoi(optarg); break;
    case 50: in->verbose = atoi(optarg); break;
    case 99:
      printf("`--mesh (required), Name of the input mesh (.re2 file)\n");
      printf("--tol (optional, default = 0.2), Tolerance used for mesh "
             "connectivity.\n");
      printf("--test (optional), Run tests in `genmap`/`gencon` examples.\n");
      printf("--dump (optional), Dump `.co2`/`.ma2` file(s) after running "
             "`gencon`/`genmap`.\n");
      printf("--nactive (optional, default: `INT_MAX`), Number of active "
             "processes.\n");
      exit(EXIT_SUCCESS);
      break;
    default: exit(EXIT_FAILURE); break;
    }
  }

  if (in->mesh == NULL) {
    fprintf(stderr, "Required argument `--mesh` was not found !\n");
    exit(EXIT_FAILURE);
  }

  return in;
}

static void parrsb_cmd_opts_free(parrsb_cmd_line_opts *opts) {
  if (opts) {
    if (opts->mesh) free(opts->mesh);
    free(opts);
  }
}

void parrsb_check_error_(int err, char *file, int line, MPI_Comm comm);
#define parrsb_check_error(err, comm)                                          \
  parrsb_check_error_(err, __FILE__, __LINE__, comm)

#endif
