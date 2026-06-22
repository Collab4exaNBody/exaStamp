# exaStampMDStressLab plugin

Wraps [MDStressLab++](https://github.com/AmitAcc/MDStressLab) to compute Hardy
(Cauchy) stress fields from an exaStamp simulation, using `process_dedr` and/or
projected-force methods evaluated through OpenKIM models.

The `compute_stress_grid_mdstresslab` operator (`compute_stress_grid_mdstresslab.cpp`)
has two modes, selected by the `read_from_file` slot (default `false`):

- **External file mode** (`read_from_file: true`, `filename` slot set): reads
  an MDStressLab-format configuration file and computes the stress field on a
  fixed analysis grid.
- **Internal grid mode** (default, `read_from_file: false`): builds an
  MDStressLab `BoxConfiguration` directly from the current exaStamp
  `grid`/`domain` (current configuration) and `grid_t0`/`xform_t0` (reference
  configuration), then computes the stress field over the local sub-domain.

## Recent fixes

- **Mode selection decoupled from `filename`**: the operator used to switch
  to external-file mode whenever `filename.has_value()` was true. In a
  pipeline, the `filename` input slot can end up auto-connected to an
  unrelated operator's output slot of the same name, which forced file mode
  even when internal-grid mode was intended. Mode is now controlled
  explicitly via the `read_from_file` boolean slot (default `false`, i.e.
  internal-grid mode).

- **Reference-configuration copy loop bound**: the loop copying particle
  positions from `grid_t0` into `body.coordinates[Reference]` iterated using
  the *current* cell's particle count (`n`) instead of the *reference* cell's
  particle count (`n0`). When the two differ (the normal case once particles
  have moved between cells), this could read out of bounds from `grid_t0` or
  silently drop reference coordinates. The loop now iterates `0..n0`.

- **Hardcoded particle species**: `body.species[iloc]` was unconditionally set
  to `"Ta"`, regardless of the actual particle species, which produced
  incorrect stress fields for any non-Ta or multi-species system. It is now
  set from the particle's `field::_type` via the `species` slot
  (`(*species)[type[j]].name()`).

- **`grid_t0` / `xform_t0` slots**: these were declared `OPTIONAL` but
  dereferenced unconditionally in internal-grid mode, which would crash if
  unset. They are now `REQUIRED`, matching how `grid`/`domain` are already
  declared.

## Known limitations (not yet addressed)

- External-file mode hardcodes the analysis grid resolution (`60x5x60`) and
  bounding box (`(0,0,0)`-`(217.8,16.5,198.0)`), which match a specific test
  case rather than being general parameters.
- The Hardy sphere radius (`MethodSphere hardy(6, "hardy")`) is hardcoded in
  both modes and does not feed into `rcut_max` / ghost-layer sizing.
- Several declared slots (`grid_subdiv`, `mesg`, `rcut_max`) are currently
  unused in `execute()`.
