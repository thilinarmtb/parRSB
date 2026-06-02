#include <parrsb_impl.h>

INTERN void mpi_topology(unsigned *num_levels, struct comm *const comms,
                         const struct comm *const c, const int verbose) {
  uint requested = *num_levels;

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
  parrsb_info(c, verbose, "parRSB: levels requested %u, enabled = %u",
              requested, *num_levels);
}
