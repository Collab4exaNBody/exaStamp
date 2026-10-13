# Physics regression collection

One case per interatomic potential (`<case>/run.msp`). Each case builds a small system (lattice or configuration
file), adds deterministic noise to positions and velocities, runs 10 steps and compares the final state to two
**committed** references:

- `ref_values.dat`: positions, velocities and accelerations of sampled particles (`check_values`);
- `ref_thermo.yaml`: potential, kinetic, rotational and total energies, virial tensor, volume and particle count
  (`check_thermodynamic_state`).

A test fails when the physics changes: different forces or energies, or a missing reference.

## Shared settings: `common.yaml`

- `deterministic_noise: true`: the noise is drawn per domain cell, so seq, mt and par start from the same state and
  compare against the same references.
- `physics_regression_noise`: noise block appended to each case's `setup_system` (gaussian on positions, uniform on
  velocities).
- `simulation_epilog`: recomputes forces and energies, then runs both checks. Thresholds are the globals
  `pr_pos_threshold`, `pr_vel_threshold`, `pr_acc_threshold` (absolute L2 errors of `check_values`, internal units) and
  `pr_thermo_threshold` (relative). A case can override them in its `+global:` block.

## Running

```bash
ctest -R physics_regression            # all cases, seq / mt / par
ctest -R physics_regression_eam_alloy  # one family
```

Each case runs in `<build>/regression_runs/physics_regression/<case>/`.

## (Re)generating references

Regenerate references only for an intentional change of the physics, and commit the reference diff together with
the change that explains it:

```bash
./make_references.sh <exaStamp launcher> <case ...|all>
# e.g. ./make_references.sh ~/local/exaStamp/bin/exaStamp pair_lj_single
```

The script runs each case once (1 MPI process, 1 thread) in a temporary copy with `create_references.msp`, then copies
both reference files back into the case folder.

## Adding a case

Create `<case>/run.msp`: include `../common.yaml`, define species, the potential, `compute_force` and `setup_system`
(domain, lattice or configuration, then `physics_regression_noise`), plus `rcut_inc` and `dt` in `+global:`. Potential
files are taken from `data/potentials`. Then run `make_references.sh` for the case.

## Disabled variants

The first line of a `run.msp` may hold `no-seq`, `no-mt` or `no-mpi`. The second line explains why a variant is
disabled; these are known problems still to be fixed (see the comment in each case).
