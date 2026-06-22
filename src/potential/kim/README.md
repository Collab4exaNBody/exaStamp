# exaStampKIM plugin

Wraps [OpenKIM](https://openkim.org/) models so they can be used as a pair
potential in exaNBody pipelines. The model is initialized once
(`kim_init`, `kim_init.cu`) and then evaluated per-particle, per-timestep
(`kim_force`, `kim_force.cu` / `kim_force_op.h`).

- `kim_init`: creates a temporary `KIM::Model`, checks that the model
  supports every species declared in `species`, records the per-species
  `kim_particle_codes`, reads the model's neighbor-list cutoffs (diagnostic
  `kim_ctx->rcut`) and influence distance (used as `parameters->rcut`, the
  cutoff that actually drives neighbor search), then destroys the temporary
  model.
- `kim_force`: lazily creates one `KIM::Model` per OpenMP thread
  (`KIMThreadContext::kim_model`, stored in `kim_ctx->m_thread_ctx`) and runs
  `compute_cell_particle_pairs` with `KimForceOp` (`kim_force_op.h`), which
  builds a small "central atom + neighbors" cluster and calls
  `KIM::Model::Compute` on it for every particle.

## Recent fixes

- **Central-atom species code** (`kim_force_op.h`): the central particle's
  KIM species code was set directly from exaStamp's internal species index
  (`species_codes[0] = type;`), while neighbor species were correctly
  translated through `kim_particle_codes[]`. The central atom now uses the
  same translation (`species_codes[0] = kim_particle_codes[type];`), so it is
  evaluated as the right species when a model supports multiple species with
  non-trivial KIM species codes.

- **`kim_local_model` leak** (`kim_init.cu`): the temporary `KIM::Model`
  created to query species support, cutoffs, and influence distance was never
  destroyed. `KIM::Model::Destroy(&kim_local_model)` is now called at the end
  of `execute()`.

- **Discarded `ComputeArguments` objects** (`kim_init.cu`, `kim_force.cu`):
  both files created a `KIM::ComputeArguments` object via
  `ComputeArgumentsCreate` that was never used and never destroyed (one leak
  in `kim_init`, one per OpenMP thread in `kim_force`'s thread-context setup).
  Both have been removed; thread-local `ComputeArguments` are created and
  destroyed where they are actually used, in `KimForceOp::operator()`.

- **Dead/commented-out code** (`kim.h`, `kim_init.cu`, `kim_force.cu`):
  removed the unused `KIMContext::kim_model` member (only
  `KIMThreadContext::kim_model` is used), commented-out `m_test`
  declarations, a commented-out `rcut` YAML decode, and ~60 lines of
  commented-out validation/diagnostic code in `kim_force.cu` that duplicated
  checks already done in `kim_init`.

- **Misc**: fixed a signed/unsigned comparison and a duplicate semicolon in
  the species-code printing loop in `kim_init.cu`. Added a comment in
  `kim.h` clarifying that `KIMContext::rcut` (max of the model's neighbor-list
  cutoffs) is diagnostic-only, while `KIMParams::rcut` (the influence
  distance) is the cutoff that drives `compute_cell_particle_pairs`.

## Known limitations (not yet addressed)

- **`modelWillNotRequestNeighborsOfNoncontributingParticles`** is fetched in
  `kim_init.cu` but never checked. `KimForceOp` builds a "central + direct
  neighbors only" cluster and marks only the central particle as
  contributing; if a model needs neighbor-of-neighbor information for
  non-contributing particles, this per-atom cluster approach may silently
  produce incorrect results.

- **Per-atom allocations and `ComputeArguments` churn** (`kim_force_op.h`):
  `KimForceOp::operator()` allocates several `std::vector`s (`coords`,
  `species_codes`, `contributing`, `energies`, `forces`, `virials`, the
  `CentralOnlyNL` neighbor list) and calls
  `ComputeArgumentsCreate`/`ComputeArgumentsDestroy` for every particle, every
  timestep. This is correct but costly; a future optimization is to give each
  `KIMThreadContext` persistent, pre-sized scratch buffers and a single
  `KIM::ComputeArguments*` set up once per thread and reused across particles.

- **Mixed error-handling style**: `MY_ERROR` (defined in `kim.h`, calls
  `exit(1)` directly) is used alongside onika's `fatal_error()`/`lerr`
  streams.
