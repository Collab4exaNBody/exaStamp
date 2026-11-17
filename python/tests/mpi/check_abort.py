"""Run under mpirun with 2+ processes: an uncaught Python exception on one rank
must abort every process (non-zero exit) instead of leaving the others blocked."""

import pyexastamp as xs

ctx = xs.init(xs.command_line(xs.find_example("regression_new/potentials/eam/eam_vniitf/single_specy.msp"),
                              program="check_abort", overrides={"global": {"max_iteration": 1}}))
if ctx.mpi_rank == 1:
    raise RuntimeError("deliberate failure on rank 1")
xs.run(ctx)   # rank 0 waits for rank 1 in collective operations
xs.end(ctx)
