"""Run under mpirun with 2+ processes: results must be consistent across ranks.

A failed check raises an exception, which aborts every process (non-zero exit).
"""

import numpy as np

import pyexastamp as xs

MSP = xs.find_example("regression_new/potentials/eam/eam_vniitf/single_specy.msp")


def run(gather):
    ctx = xs.init(xs.command_line(MSP, program="check_mpi", overrides={
        "global": {"max_iteration": 3},
        "simulation_epilog": [{"extract_particle_fields": {"fields": ["rx", "virial"], "gather": gather}}],
        "grid_flavor": "grid_flavor_full",
    }))
    assert ctx.error_code < 0
    assert ctx.mpi_nprocs > 1, "run with mpirun -n 2 or more"
    xs.run(ctx)
    thermo = xs.read_thermo(ctx)
    p = xs.particles(ctx)
    rank, nprocs, n = ctx.mpi_rank, ctx.mpi_nprocs, thermo["particle_count"]
    xs.end(ctx)
    return rank, nprocs, n, thermo["temperature"], p


# local copies: every rank holds a part of the particles
rank, nprocs, n, temperature, p = run(gather=False)
local = len(p["id"])
assert 0 < local < n, (rank, local, n)
assert p["virial"].shape == (local, 3, 3)
print(f"rank {rank}/{nprocs}: {local} of {n} particles, T={temperature:.6f}", flush=True)

# gathered copy: rank 0 holds every particle once, the others nothing
rank, nprocs, n, temperature, p = run(gather=True)
if rank == 0:
    assert len(p["id"]) == n, (len(p["id"]), n)
    assert len(np.unique(p["id"])) == n
    assert p["rx"].shape == (n,) and p["virial"].shape == (n, 3, 3)
else:
    assert len(p["id"]) == 0 and p["virial"].shape == (0, 3, 3)
print(f"rank {rank}/{nprocs}: gather ok", flush=True)
