#!/usr/bin/env python3
"""Run an exaStamp EAM sample and inspect the results from Python.

Usage: python run_eam_sample.py [input.msp]
       mpirun -n 4 python run_eam_sample.py [input.msp]

With several MPI processes, particles are gathered on rank 0, which prints.
"""

import sys

import numpy as np
import pyexastamp as xs

msp = sys.argv[1] if len(sys.argv) > 1 else xs.find_example("regression_new/potentials/eam/eam_vniitf/single_specy.msp")

# copy every exportable per-particle field at the end of the run, gathered on rank 0
ctx = xs.init(xs.command_line(msp, overrides={
    "simulation_epilog": [{"extract_particle_fields": {"gather": True}}],
}))
if ctx.error_code >= 0:
    sys.exit(ctx.error_code)
xs.run(ctx)

if ctx.mpi_rank != 0:   # results are printed once, by rank 0
    xs.end(ctx)
    sys.exit(0)

thermo = xs.read_thermo(ctx)
print("\n=== thermodynamic state ===")
for key in ("particle_count", "temperature", "pressure", "total_energy", "potential_energy", "volume", "density"):
    print(f"  {key:18s} {thermo[key]:.6g} {thermo['units'].get(key, '')}")

print("\n=== species ===")
for s in xs.read_species(ctx):
    print(f"  {s['type']}: {s['name']}  mass={s['mass']} Da  z={s['z']}  charge={s['charge']} e")

atoms = xs.particles(ctx)
if atoms is not None:
    print("\n=== particles ===")
    for name, a in atoms.items():
        print(f"  {name:8s} dtype={a.dtype}  shape={a.shape}")
    r = np.stack([atoms["rx"], atoms["ry"], atoms["rz"]], axis=1)
    print(f"  centre of mass: {r.mean(axis=0)}")

xs.end(ctx)
