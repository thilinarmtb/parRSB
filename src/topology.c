#include <parrsb_impl.h>

INTERN void mpi_topology(uint *num_levels, struct comm *const comms,
                         const struct comm *const c, const int verbose) {
  // Level 1 communicator is the global communicator.
  comm_dup(&comms[0], c);

  // Level 2 communicator is the node level communicator.
  MPI_Comm node;
  MPI_Comm_split_type(c->c, MPI_COMM_TYPE_SHARED, c->id, MPI_INFO_NULL, &node);
  comm_init(&comms[1], node);
  MPI_Comm_free(&node);

  // Set num levels to 2.
  *num_levels = 2;

  // Sanity check: rpn should be the same across all the ranks.
  sint max = comms[1].np, min = comms[1].np, wrk;
  comm_allreduce(&comms[0], gs_int, gs_max, &max, 1, &wrk);
  comm_allreduce(&comms[0], gs_int, gs_min, &min, 1, &wrk);
  parrsb_error(max == min, c, verbose,
               "ranks per node is not the same across all the nodes.");

  // Sanity check: ranks per node must be positive.
  parrsb_error(min > 0, c, verbose, "ranks per node must be positive.");

  // Find the number of nodes under the global communicator.
  sint num_nodes = (comms[1].id == 0);
  comm_allreduce(c, gs_int, gs_add, &num_nodes, 1, &wrk);
  parrsb_info(c, verbose, "parRSB: nodes = %u, ranks per node = %u", num_nodes,
              min);
  parrsb_info(c, verbose, "parRSB: levels enabled = %u", *num_levels);
}

INTERN void slingshot_topology(uint *num_levels, const char *const template,
                               struct comm *const comms,
                               const struct comm *const c, const int verbose) {
  char *modified = calloc(3 * strlen(template) + 1, sizeof(char));
  const char *sub = "%u%n";

  // Find rack, cabinet, router, blade, node, etc. based on format string.
  uint size = 0, sidx = 0, didx = 0;
  while (template[sidx] != '\0') {
    // replace "%d" by "%u%n\0"
    if (strncmp(&template[sidx], "%d", 2) == 0) {
      strcpy(&modified[didx], sub);
      size++, sidx += 2, didx += strlen(sub) + 1 /* for '\0' */;
    } else {
      modified[didx++] = template[sidx++];
    }
  }
  modified[didx] = '\0';

  int hname_len;
  char hname[MPI_MAX_PROCESSOR_NAME];
  MPI_Get_processor_name(hname, &hname_len);

  uint hidx = 0, idx = 0;
  uint *elements = calloc(size + 1, sizeof(uint));
  for (uint i = 0; i < size; i++) {
    uint offset;
    sscanf(&hname[hidx], &modified[idx], &elements[i], &offset);
    hidx += offset;

    while (modified[idx] != '\0') idx++;
    idx++;
  }
  free(modified);

  // Level 1 communicator is the global communicator.
  comm_dup(&comms[0], c);
  for (uint i = 0; i < size; i++)
    comm_split(&comms[i], elements[i], elements[i + 1], &comms[i + 1]);
  free(elements);

  *num_levels = size + 1;
  parrsb_info(c, verbose, "parRSB: levels enabled = %u", *num_levels);
}
