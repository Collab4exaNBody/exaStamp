# Delaunay / PTM / DXA pipeline

Goal: build towards a DXA-style dislocation extraction pipeline in exaStamp, following
Stukowski, Bulatov, Arsenlis, *"Automated identification and indexing of dislocations in
crystal interfaces"*, Model. Simul. Mater. Sci. Eng. 20 (2012) 085007
(local copy: `/home/lafourcadep/Bureau/DXA_ref.pdf`).

DXA's published pipeline has 9 steps. Status:

| Step | What it does | Status |
|---|---|---|
| (i) | Atomic structure identification: crystal type **and** local lattice orientation per atom | crystal type: done (`compute_slcsa`, `src/analysis_particle/`). Orientation: **done** via PTM, see below |
| (ii) | Space-filling Delaunay tessellation | **done** — this directory |
| (iii) | Assign an ideal lattice vector to each tessellation edge | **done** — `compute_dxa_edge_vectors`, `src/ptm/`, see below |
| (iv) | Classify each tetrahedron good/bad from edge compatibility | **done** — `compute_dxa_tet_classification`, `src/ptm/`, see below |
| (v) | Build the interface mesh (2D manifold separating good/bad regions, half-edge structure) | **done** — `compute_interface_mesh`, `src/delaunay/`, see below |
| (vi)-(vii) | Burgers circuit construction, real-vs-noise defect filtering | **done** — `compute_dxa_burgers_circuits`, `src/delaunay/`, see below |
| (viii)-(ix) | Sweep along the interface mesh to extract dislocation line geometry, junction detection | not started |

## This directory: Delaunay tessellation (done, verified)

`compute_delaunay.cpp` builds a 3D Delaunay tessellation of each MPI rank's owned+ghost
particles using [Geogram](https://github.com/BrunoLevy/geogram)'s `GEO::PeriodicDelaunay3d`,
constructed in **non-periodic** mode. exaStamp's own ghost particles are already literal
duplicated position copies for both periodic images and cross-rank neighbors, so Geogram never
needs to handle periodicity itself — this sidesteps `PeriodicDelaunay3d`'s own periodicity
machinery, which is designed for a single complete point set, not domain-decomposed data.

**Parallelism**: MPI parallelism comes from exaStamp's own domain decomposition (one independent
tessellation per rank, nothing new needed). OMP parallelism comes from Geogram's own internal
threading — confirmed the installed `libgeogram.so` links `libgomp` and calls
`omp_get_max_threads()` internally.

**Ownership / ghost-fringe-trust rule**: a tetrahedron is kept by a rank only if its centroid
falls in one of that rank's own **owned** (non-ghost) cells. This single rule handles two
concerns at once:
- No gaps or duplicates once every rank's output is combined (every point in space belongs to
  exactly one rank's owned region).
- Discards exactly the tetrahedra that could be geometrically unreliable — Delaunay has no fixed
  radius of influence the way a physical interaction cutoff does, so a tet near the outer fringe
  of the ghost halo isn't guaranteed correct; those can only occur far from any owned cell, which
  this rule already excludes.

**Real bug found and fixed**: for a rank whose ghost-inclusive local block straddles a periodic
domain edge (e.g. a 2-rank split along i: rank 0 owns global cells `[0,5)`, `ghost_layers=2`
gives `block_start.i=-2`), exaNBody's `domain_periodic_location()` always wraps into the domain's
*canonical* `[0,domain_dims)` range — the wrong frame to subtract a possibly-negative
`block_start` from. Fixed with a custom per-axis resolver (`resolve_local_axis` in
`compute_delaunay.cpp`) that computes the unwrapped global cell coordinate directly and shifts by
`±domain_dims` only as needed to land inside this rank's own local extent.

**Verified**: total tetrahedron count is exactly **96000**, identically, whether run with 1, 2,
or 4 MPI ranks (`data/regression_new/delaunay/compute_delaunay.msp`, 16000-atom noisy BCC Ta
lattice). No gaps, no duplicates.

`write_delaunay_vtk.cpp` writes each rank's tessellation as a `.vtu` piece plus a `.pvtu` master
file (same convention as exaNBody's `write_grid_vtk.cpp`), with a per-tetrahedron `rank` CellData
field so Paraview can color by MPI rank. Open the `.pvtu`, not the individual `.vtu` pieces.

## `src/ptm/`: PTM strategy

DXA step (iii) needs, per atom, more than a structure-type label — it needs a **local lattice
orientation**: a one-to-one mapping from the atom's actual neighbor bonds to the ideal lattice's
template directions. `compute_slcsa`'s bispectrum-based classifier is deliberately
**rotation-invariant** (that's what makes it a good structure-type classifier), so it structurally
cannot recover orientation — that's a separate registration/fitting problem.

First attempt (now deleted) hand-rolled this: iteratively assign each bond to its nearest ideal
BCC `<111>` template direction under a running rotation estimate, recompute the rotation as the
orthogonal polar factor of the resulting cross-covariance matrix (reusing the same Newton
iteration as `compute_polar_decomposition`'s deformation-gradient polar decomposition). **This was
the wrong call** — it conflated two unrelated numerical problems that only happen to share the
same underlying polar-decomposition math, and re-derived from scratch something that already has
an established, validated implementation: **PTM (Polyhedral Template Matching)**, Larsen, Schmidt
& Schiøtz, Model. Simul. Mater. Sci. Eng. 24 (2016) 055007. PTM uses Theobald's QCP
(quaternion-characteristic-polynomial) method for the rotation fit, plus proper combinatorial
template/graph matching for correspondence — not a greedy heuristic.

**PTM is vendored, not linked.** OVITO bundles PTM's source directly (no separate installable
library exists), so it's copied from `/home/lafourcadep/CODES/VISU/ovito/src/3rdparty/ptm/`
(MIT license, preserved as `src/ptm/PTM_LICENSE`) into:
- `src/ptm/lib/` — 17 vendored `.cpp` files, auto-compiled by `onika_add_plugin` into an internal
  `exaStampPTM` shared library.
- `src/ptm/include/` — 22 vendored headers, flat (matches upstream's own flat layout).

**PTM does NOT need Voro++**, despite OVITO's own `CMakeLists.txt` linking `VoroPlusPlus` — every
PTM source file was grepped for real Voro++ includes / `voro::` namespace usage (as opposed to
`ptm_voro::`, PTM's own self-contained internal Voronoi-cell implementation): zero hits. That
link in OVITO's build is vestigial. (Voro++'s files were copied once, then deleted after
confirming this — don't re-add them.)

`src/ptm/compute_ptm.cu` is the new operator, scoped to **BCC only** for now (`PTM_CHECK_BCC`,
matches this session's only test lattice — other `PTM_CHECK_*` flags for FCC/HCP/ICO/SC would
work the same way). Per particle: gather neighbors via the normal pair-compute loop, sort by
distance, build a `ptm_atomicenv_t` (central atom at the origin, neighbors as relative vectors —
confirmed via `ptm_normalize_vertices.cpp` that PTM re-centers on the barycentre of all points
internally, so relative vectors are fine), call `ptm_index()` through a trivial callback bridge
(valid because PTM's neighbour-ordering for BCC/FCC/HCP/ICO/SC calls the callback exactly once
per atom), and write out structure type, the fitted orientation quaternion (as-is, w,x,y,z), and
RMSD, to flat per-particle output buffers (not named grid fields — see below). One
`ptm_local_handle_t` scratch object per OpenMP thread (not
thread-safe to share) — same per-thread-context pattern the old `compute_local_structural_metrics.cpp`
used for its `SnapLegacyBS` instances. PTM itself is host-only (combinatorial graph/template
matching, not GPU-portable) — the functor is intentionally `CudaCompatible=false`.

## `compute_ptm.cu`: compiles, verified

Three compile errors were fixed getting here (all in `compute_ptm.cu`'s git history):
1. The bufferless `ComputePairParticleContextStart/Stop` pattern doesn't work with
   `CudaCompatible=false` + chunk-neighbor iteration. Switched to the buffer-based pattern
   (`operator()(jnum, buf, cells)`, matching `eam_force_op_singlemat.h`'s `CudaCompatible=false`
   precedent).
2. A real bug: the functor's `operator()` had 4 params `(jnum, buf, itype, cells)`, copying
   `BispectrumOpRealT`'s shape without also copying its `FieldSet<field::_type>` requirement —
   with an empty field set the caller only ever passes 3 arguments. Fixed by dropping `itype`.
3. Writing named grid fields via `cells[cell][field][part]` from the buffer-based pattern needs
   `compute_cell_particle_pairs2`'s deep `cells_accessor()` (`force_use_cells_accessor=true`) —
   but forcing that exposed the real, deeper problem below, so this operator doesn't write named
   grid fields at all in the end (see next point).

**The actual blocker, found and fixed**: `compute_cell_particle_pairs2`'s own `CS=1`/`CS=VARIMPL`
dispatch macro (`compute_cell_particle_pairs.h`) unconditionally compiles *both* branches for any
chunk-neighbor caller. The `CS=VARIMPL` (runtime chunk size) branch cannot compile for a
`CudaCompatible=false` functor on this exaNBody version — none of `compute_cell_particle_pairs_
cell`'s three overloads (`chunk.h`/`chunk_cs1.h`/`chunk_scb.h`) accept a genuinely runtime chunk
size, only the compile-time `onika::UIntConst<CS>`. This is apparently untested for
`CudaCompatible=false` on this GPU build specifically: `eam_force_op_singlemat.h`'s own
`CudaCompatible=false` branch is dead code here (`USTAMP_POTENTIAL_CUDA_COMPATIBLE` resolves true
in a GPU build), and `compute_bispectrum.cu`'s functor is `CudaCompatible=true` so never exercises
this path either — confirmed by rebuilding both, one clean, and by confirming the error persists
identically whether or not `force_use_cells_accessor` is set (ruling out the cells-accessor
question as the actual cause). PTM is the first genuinely `CudaCompatible=false` chunk-neighbor
buffer-based operator to compile in this GPU build.

**Fix**: bypass `compute_cell_particle_pairs2`'s broken dispatch macro entirely — call
`compute_cell_particle_pairs_cell` directly with a hardcoded `onika::UIntConst<1>{}`, in our own
`#pragma omp parallel for collapse(3)` loop over owned cells (mirrors exactly what
`ComputeParticlePairFunctor::operator()` does internally on the host path). This requires
`chunk_neighbors` to actually be built with `chunk_size==1`; `execute()` checks this at runtime
(`fatal_error()` otherwise) — set `chunk_neighbors: { config: { chunk_size: 1 } }` in the `.msp`.

**Output shape changed along the way** (user's call, and the right one): instead of named grid
fields, `compute_ptm` writes three flat `onika::memory::CudaMMVector<double>` OUTPUT slots
(`ptm_struct_type`, `ptm_orientation` — quaternion w,x,y,z, `ptm_rmsd`), indexed via
`grid->cell_particle_offset_data()`, exactly `compute_bispectrum.cu`'s own pattern. Sidesteps the
cells-accessor question entirely, and matches how `compute_slcsa` already consumes bispectrum's
flat buffer directly rather than via a named field — DXA step (iii) can do the same with PTM's
orientation buffer. Note: only *owned* cells are touched by the loop, so ghost-particle slots in
these buffers are left at their `CudaMMVector::resize()` default, not a real PTM result.

**Verified** (`data/regression_new/delaunay/compute_ptm.msp`, same 16000-atom BCC Ta lattice as
the Delaunay test, unrotated/unnoised): all 16000 owned particles matched `PTM_MATCH_BCC` with
`mean rmsd=0.000000` and `mean |angle from identity|=0.000000 deg` — exact match, as expected for
a perfect unrotated lattice. (An operator whose only outputs are non-grid-field OUTPUT slots gets
silently pruned from the graph if nothing downstream consumes them — the `.msp` wires
`ptm_n_matched` through `print_int` via the `rebind`+`body` "implicit batch" trick to keep it
alive; watch for a `*SUPPRESS* ...compute_ptm` startup line if this ever goes quiet again.)

**Two real bugs found via rcut-sensitivity and a cross-check against OVITO's own PTM** (user
observation: changing `rcut` on the BCC test shouldn't matter once it's generous enough, and a
separate `compute_ptm_fcc.msp` FCC lattice was getting 0% matches when OVITO's PTM modifier gets
100% FCC / max rmsd~0.05 on the same generated `fcc_sample.xyz`):
1. The neighbor-selection sort's inner loop was bounded by `n` (the target count), not `jnum` (the
   real number of candidates) — so it only ever sorted whichever `n` neighbors happened to land
   first in the buffer (traversal order, not distance order), not the true nearest ones. A larger
   `rcut` changed which neighbors filled those slots, making results silently depend on `rcut`.
   Fixed: inner loop now scans the full `jnum` range.
2. No RMSD acceptance gate: PTM always returns its *best* combinatorial match, however bad — it
   never refuses on its own, that's the caller's job. Without a cutoff, genuinely non-BCC atoms
   (a real FCC lattice, checked only against the BCC template) still got *some* correspondence,
   just a terrible one (rmsd~0.26, orientation ~44° off) — and got silently accepted as "BCC".
   Added `rmsd_cutoff` (default `0.1`, same as OVITO's own PTM modifier default): a match above it
   is reported as `PTM_MATCH_NONE` (rmsd/orientation still written, for diagnostics).

**Structure coverage widened**: `compute_ptm` now checks FCC, HCP, BCC, ICO and SC in one pass
(`PTM_MATCH_*` 1–5) — all of PTM's *single-shot* neighbour-ordering structures (`ptm_index.cpp`'s
`num_outer==0` branch, the one `ptm_get_neighbours_from_prebuilt_env` actually implements).
DCUB/DHEX/GRAPHENE are deliberately excluded: `ptm_index.cpp` calls those through a genuinely
different *two-shell* ordering (separate inner/outer-shell requests back to the neighbour
callback), and our callback always hands back the same single prebuilt env regardless of what's
requested — enabling those flags as-is would silently produce wrong results, not just "unsupported
right now". Needs the callback taught the two-shell protocol first.

One subtlety worth remembering if the candidate-neighbor cap is ever revisited: `ptm_index.cpp`
normalizes (barycentre + mean-length scale) *all* `env.num` supplied points once, up front, before
any per-structure check runs — so an overly generous env size can, in principle, bias every
structure's fit through that shared normalization, not just add harmless unused points. In
practice this build's tests (BCC checked with a 14-point vs. an 18-point pool) showed no
measurable difference, so the cap was left at `PTM_MAX_INPUT_POINTS-1=18` (OVITO's own value,
`PTMAlgorithm.cpp`) rather than re-deriving a tighter one — flag if a future structure/lattice
shows the opposite.

**Verified again** with both bugs fixed and FCC/HCP/BCC/ICO/SC all enabled: BCC lattice unchanged
(16000/16000, rmsd=0.0257), FCC lattice (`compute_ptm_fcc.msp`, 32000 owned particles) now matches
32000/32000 as `PTM_MATCH_FCC` with mean rmsd=0.033 — consistent with OVITO's own max~0.05 on the
identical generated `fcc_sample.xyz`. A separate SC test (`compute_ptm_sc.msp`) initially got 0
matches too — root cause there was purely `.msp` wiring, not code: `compute_force: nop` means
nothing establishes `rcut_max` before the neighbor list is built, and `compute_ptm` only runs in
`simulation_epilog`, after that list already exists, so its own `rcut` can't retroactively widen a
candidate pool that was never gathered. Fixed by declaring `global: { rcut_max: 8.0 ang }`
explicitly; then SC matched 63055/64000 owned (~98.5%, the rest lost to noise breaking SC's strict
degree-4 topology check — expected, not a bug).

## `ptm_fields` + `write_delaunay_vtk`'s `color_field`: materialize PTM into fields, visualize on the tessellation

`compute_ptm` writes flat `CudaMMVector` buffers, not grid fields (deliberately, to dodge the whole
cells-accessor/CS-dispatch mess above) — but for anything that expects a real per-particle grid
field (`write_xyz`, `write_grid_vtk`, or coloring the Delaunay tessellation), that's the wrong
shape. Two new pieces close this gap:

- **`ptm_fields`** (new operator, `src/ptm/compute_ptm.cu`): copies `ptm_struct_type`/
  `ptm_orientation`/`ptm_rmsd` into named grid fields (`field::mk_generic_real`/`mk_generic_mat3`,
  default names `ptm_type`/`ptm_orientation`/`ptm_rmsd`) — quaternion is converted to a rotation
  tensor here. Purely pointwise, no neighbor search, so it uses `exanb::compute_cell_particles`
  (`compute_cell_particles.h`) instead of the pair-compute machinery `compute_ptm.cu` itself needed
  — a completely different, much simpler code path with none of the CS/VARIMPL baggage, since
  there's no chunk-neighbor iteration involved at all. Run any time after `compute_ptm`.

- **`write_delaunay_vtk`'s new `color_fields` slot** (list, originally a single `color_field`,
  generalized on request so several fields can be visualized from the same output without
  re-running the pipeline): generic, not PTM-specific — takes any number of per-particle
  `field::mk_generic_real` field names and writes one extra `Float64` `CellData` array per field,
  each tet's value averaged over its 4 vertices. Needed `DelaunayTessellation` to grow a
  `vertex_particle_index` array (`compute_delaunay.cpp`, populated via
  `grid->cell_particle_offset_data()` in the same point-building loop that already exists) mapping
  each compacted vertex back to a flat particle index — without it there was no way to look up any
  per-particle field for a given tetrahedron's corners. `write_delaunay_vtk` had to become
  grid-variant (`make_grid_variant_operator`, was `make_simple_operator`) to gain a `grid` input at
  all. Leaving `color_fields` empty preserves the exact old behavior (verified: `.pvtu` has no
  extra `PDataArray` when empty).

  **Boundary-fluctuation caveat, found and fixed** (user-reported: `write_xyz`'s own `ptm_type`
  output was correct everywhere, but the Paraview/VTK tessellation showed fluctuations right at
  rank boundaries — exactly the ghost-copy issue predicted above, since `write_xyz` only ever
  writes *owned* particles while a boundary tet's vertices can include ghost copies). exaNBody
  already has a generic mechanism for this: `ghost_update_opt: { opt_fields: [ "ptm_type" ] }`
  added to the `.msp` pipeline (after `ptm_fields`, before `compute_delaunay`) explicitly
  synchronizes that dynamic field's ghost copies from their owning rank, same as any other
  ghost-communicated field — confirmed this removes the boundary artifacts. Add every
  `ptm_fields`-produced field name that ends up used by `color_fields` to `opt_fields` the same way.

**Verified end-to-end** (`compute_ptm_delaunay_color.msp`, BCC lattice): pipeline
`compute_ptm → ptm_fields → ghost_update_opt: { opt_fields: [ptm_type, ptm_rmsd] } →
compute_delaunay → write_delaunay_vtk: { color_fields: [ptm_type, ptm_rmsd] }` runs clean —
`ptm_type`/`ptm_rmsd` both appear as separate selectable `CellData` arrays in the same `.pvtu`,
`ptm_type` is `3` (pure BCC) uniformly including at rank boundaries — no more fractional values
from unsynced ghost copies.

## `compute_dxa_edge_vectors`: DXA step (iii)

Assigns an ideal lattice vector to each Delaunay tessellation edge, against a single **user-chosen
`target_structure`** (`FCC`/`HCP`/`BCC`/`ICO`/`SC`) — not whatever PTM happened to locally match.
User's requirement: this lets a specific structure be probed even in a system with mixed regions
(e.g. classify only the BCC matrix around a precipitate, or only the FCC phase in a two-phase
system), rather than the classification silently changing definition from one atom to the next.

Per edge `(u,v)`:
- both endpoints must have `PTM_MATCH_<target_structure>` as their own `ptm_type` (from
  `ptm_fields`) — otherwise unresolved.
- rotate the actual bond vector into each endpoint's own ideal-lattice frame using its
  `ptm_orientation` rotation tensor (`R` maps ideal→actual, so actual→ideal is `R^T` for a
  rotation — same convention as the "grid space vs physics space" question earlier: `domain->
  xform()` is applied to the raw bond vector first, matching the space PTM's own orientation fit
  was computed in), then find the target structure's own ideal neighbor direction — PTM's own
  reference templates (`ptm_initialize_data.h`'s `refdata_t::points`, exactly what PTM matched
  atoms against, not re-derived) closest to it.
- resolved only if both endpoints agree on the **same** ideal direction (expected within a single
  grain — PTM already returns near-identical orientations atom-to-atom there, confirmed by
  `compute_ptm`'s own BCC validation, 0.0deg mean deviation from identity) within
  `angle_tolerance` (default 15°). Disagreement is exactly the defect/grain-boundary signal step
  (iv) needs.

Edges are deduplicated across the tessellation (`DXAEdgeVectors`, new struct alongside
`DelaunayTessellation`) — many tetrahedra share an edge, so collect all 6 vertex-pairs per tet,
sort+unique, and keep a `(v0,v1)->index` map for step (iv)'s later per-tetrahedron lookups.

**File placement note**: despite being conceptually a "delaunay pipeline" step, the operator lives
in `src/ptm/` (`compute_dxa_edge_vectors.cpp`), not `src/delaunay/` — it needs `ptm_constants.h`/
`ptm_initialize_data.h` for the reference templates, and those live in a flat (non-namespaced)
include dir wired only into the `ptm` plugin. Cross-plugin `#include <exaStamp/...>` headers are
**not** globally visible by default in this build (learned the hard way: neither direction worked
until `exaStampPTM_INCLUDE_DIRS` was set in `src/ptm/CMakeLists.txt` pointing at `../delaunay/
include` — the per-plugin variable `onika_add_plugin` actually reads, mirroring the existing
`exaStampPTM_LINK_LIBRARIES` convention). Header-only, no link dependency on Geogram or the
delaunay plugin's own `.so`.

**Verified** (`compute_dxa_edges.msp`, same unrotated BCC Ta lattice): **120795/120795 edges
(100%) resolved** against `target_structure: BCC` — expected for a single perfect grain, and only
achieved once `ghost_update_opt` also synchronized `ptm_orientation` (not just `ptm_type` as in
the earlier color_fields fix), confirming boundary edges need every field they read synced. Sanity
check: the same lattice classified against `target_structure: FCC` resolved **0/120881** — confirms
the gate is real, not a no-op.

## `compute_dxa_tet_classification`: DXA step (iv)

Classifies each Delaunay tetrahedron as "good" (an undistorted patch of `target_structure`'s ideal
lattice) or "bad" (a defect: dislocation core, grain boundary, stacking fault, second phase, ...):
good means all six of a tet's edges resolved in step (iii) *and* those six ideal vectors close
consistently — going vertex 0→1→2 must match going 0→2 directly (`e12 == e02-e01`, similarly for
the other two independent relations; checking all three redundantly is cheap and catches any
bookkeeping mistake directly). No grid access needed at all — purely a post-process on
`DelaunayTessellation` + `DXAEdgeVectors` — so it's a plain (non-grid-variant) operator, unlike
every other operator in this pipeline.

Output is a new `DXATetClassification` struct (`good`: `1.0`/`0.0` per tetrahedron, parallel to
`DelaunayTessellation::tetrahedra`) — plain `double`, not `bool`, so `write_delaunay_vtk` can write
it straight out as CellData with no conversion, same convention as `ptm_type`. `write_delaunay_vtk`
gained a matching `dxa_tet_classification` INPUT slot (auto-wires in if the classification
operator ran earlier in the pipeline) that writes it as a direct `"dxa_good"` CellData array — no
per-vertex averaging needed, unlike `color_fields`, since it's already per-tetrahedron.

**Pipeline ordering matters**: `compute_dxa_edge_vectors`/`compute_dxa_tet_classification` must
run *before* the `write_delaunay_vtk` call that's meant to pick up `dxa_good` — auto-wiring only
sees a slot's current value at the point the consuming operator actually executes, not "eventually
in this pipeline". Caught this in `compute_dxa_edges.msp` (`write_delaunay_vtk` was originally
sequenced before the DXA steps, from when the file only did `color_fields` coloring); reordered.

**Verified** (`compute_dxa_edges.msp`, extended to a realistic mixed-phase system — BCC matrix
with a cylindrical FCC inclusion and a cylindrical SC inclusion, same three-lattice setup as
`analysis_particle/compute_slcsa.msp`'s own test): classifying against `target_structure: BCC`
gives 188649/217850 edges resolved and **146570/172559 tetrahedra (85%) classified good** — matches
physical expectation exactly: the BCC bulk is overwhelmingly good, while the FCC/SC inclusions
(wrong `ptm_type` there, so every one of their edges is unresolved) and the interfaces around them
come out bad. Confirmed the `dxa_good` CellData array lands correctly in the `.vtu` output with
matching counts (146570 values `== 1.0` out of 172559).

## `compute_interface_mesh` + `write_interface_mesh`: DXA step (v)

Builds the 2D surface separating "good" from "bad" tetrahedra: one triangle per tetrahedron face
shared by exactly one good and one bad tet (`InterfaceMesh`, new struct, `src/delaunay/include/
exaStamp/delaunay/interface_mesh.h`). No PTM dependency at all here (pure `DelaunayTessellation` +
`DXATetClassification` geometry), so — unlike steps (iii)/(iv) — this lives cleanly in
`src/delaunay/`, no cross-plugin include-path issue this time.

**Finding interface faces**: hash every tetrahedron's 4 faces (sorted 3-vertex key, a small custom
`FaceKeyHash` since `std::array<uint32_t,3>` has no built-in `std::hash`) and group by key. A face
referenced by exactly 2 tets is interior — check their good/bad status; if they differ, that face
becomes one interface triangle (vertex order taken from the *bad* tet's own face-vertex listing, a
fixed convention, not normal-oriented — noted as a known simplification, not needed yet). A face
referenced by only 1 tet means the neighboring tet wasn't kept in this rank's own
`DelaunayTessellation` — same "ghost-fringe-trust" concern `compute_delaunay.cpp` already reasons
about for tets themselves — so its far side's status is genuinely unknown and it's excluded rather
than guessed at.

**Edge adjacency, not full half-edge (yet)**: `InterfaceMesh::edge_triangles` maps each interface
mesh edge to the triangle indices sharing it — normally exactly 2 (the surface is topologically a
closed 2-manifold, since a dislocation line can't end inside a perfect crystal), an edge count of
1 marks where this rank's own domain-decomposition boundary cuts the surface off (not a real
defect edge), and >2 would mark a branch point. Deliberately stopped short of a strict half-edge
(twin/next/prev) structure — that specific traversal shape isn't scoped out yet, steps (vi)-(ix)
(circuit search / sweep) will dictate what's actually needed there.

**`write_interface_mesh`** exports it as a VTK triangle mesh, same per-rank-piece + `.pvtu` master
convention as `write_delaunay_vtk` — compacts to only the vertices an interface triangle actually
references (interface triangles are a small fraction of the full tessellation's vertices).

**Verified** (`compute_dxa_edges.msp`, same mixed BCC/FCC/SC system): **3510 interface triangles**,
366 cutoff edges (~7%, plausible for a periodic single-rank run — even one rank has a
self-wraparound ghost boundary), 1938 points in the compacted `.vtu` output (matches the triangle
count's own vertex references) — and the emitted triangle coordinates cluster geometrically right
around the FCC/SC cylindrical inclusions, exactly where the good/bad boundary should be.

## `compute_dxa_burgers_circuits`: DXA steps (vi)-(vii), and filtering "bad candidates"

User's ask: some "bad" regions found in step (iv) aren't real dislocations at all (isolated
misclassified tetrahedra, or in the earlier mixed-phase test, entire phase boundaries) — need a
way to tell those apart from genuine dislocation lines before trusting the interface mesh.

**First attempt didn't work — a real methodological gap, not a bug.** The natural first idea: for
every tessellation edge touching a "bad" tet, trace a small Burgers circuit — the ring of
tetrahedra sharing that edge (a fan, hinged on the edge like pages of a book) — and sum the ring's
already-resolved (step iii) ideal vectors going around it. Physically sound (a real closed loop's
holonomy is zero unless it encircles a dislocation), but tested against a real periodic BCC
dislocation quadrupole (`compute_dxa_real_case.msp`, `quadrupole_dislo.xyz`, 128000 atoms) it found
**zero** confirmed core edges. Added debug counters to find out why rather than guess:
`qualifying=18838 incomplete_ring=0 unresolved_ring=17577 zero_residual=1261 confirmed=0
max_norm_seen=0`. 93% of candidate rings touched an atom PTM couldn't classify at all — exactly
what happens right at a real dislocation core, since the core is by definition too disordered for
PTM to match anything. Any loop that actually encircles a line defect is topologically forced to
pass near its core, so a small local circuit can never see the signal — it has to go around the
defect at a safe distance, through the surrounding good lattice.

**Fix: a global elastic mapping**, built via a spanning tree over every vertex reachable through
*resolved* edges (step iii) — pick an arbitrary root, propagate an accumulated "ideal-lattice
position" (well-defined, since a tree has no cycles, so no ambiguity from path choice). Then for
**every** resolved edge in the whole mesh, not just tree edges, compare its own direct ideal
vector against the tree's accumulated position difference: tree edges trivially match (that's how
they were assigned); any other resolved edge with a mismatch above `min_burgers_norm` is a genuine
Burgers-vector signal — the tree path routes arbitrarily far around the disordered core using
whatever resolved edges exist anywhere in the good lattice, so it never needs to touch the
unresolved atoms a small ring would. Confirmed edges seed a flood-fill through bad-tet adjacency
(a real dislocation's whole tube counts as confirmed, not just the specific tets touching a signal
edge); any "bad" tet never reached is reclassified "good" in place — exactly the filtering the
user asked for.

**Verified**: quadrupole dislocation case now finds **24318 confirmed signal edges**, only
114/12297 originally-bad tetrahedra filtered as noise (the real dislocation content is preserved),
leaving **5264 interface triangles**. Re-ran the earlier mixed BCC/FCC/SC system (phase boundaries,
no real dislocations) to make sure the fix didn't lose the filtering behavior that motivated this
step in the first place: still **0 confirmed signal edges**, all bad tets reclassified good, 0
interface triangles — phase boundaries correctly still rejected as "not a dislocation" (this whole
pipeline is a *dislocation* extractor specifically, not a general defect classifier — matches
DXA's own stated scope, not a limitation introduced here).

## Mesh topology fixes, then a real reference-implementation comparison

**Two visual complaints, two different real causes.** User: interface mesh looked coarse and
"open in places" compared to OVITO's own DXA output. Investigated with two targeted fixes, both
in `src/delaunay/`:
1. `compute_dxa_burgers_circuits`'s flood-fill (grow a confirmed core edge to its whole bad-tet
   cluster) used strict *face* adjacency between bad tets. Where a dislocation's cross-section
   pinches to a thin/irregular sliver, consecutive core tets can share only an edge or vertex,
   not a full face — this was splitting one continuous tube into disconnected mesh islands.
   Switched to *vertex* adjacency for this specific purpose (kept strict face adjacency for
   `compute_interface_mesh`'s own triangle-finding, which genuinely needs it).
2. `compute_interface_mesh`'s triangles had no consistent winding (vertex order came straight
   from the bad tet's own arbitrary listing) — ParaView's backface culling/lighting can make a
   face look "missing" depending on view angle even though the mesh is topologically fine. Fixed
   by orienting every triangle so its right-hand-rule normal points away from the bad tet's own
   4th vertex (outward, into the good region), via a signed-volume check.
Verified via a connected-components analysis of the raw triangle mesh (parsing the `.vtu`
directly) that topology was identical before/after the winding fix (1 dominant 5060-triangle
component + 22 small 12-triangle satellites, 0 single-triangle edges either way) — confirming the
"open" look really was a rendering artifact, not a real gap.

**User then supplied ground truth**: `output_dxa_ovito.vtk` (OVITO's own DXA output on the same
BCC Ta dislocation quadrupole, `quadrupole_dislo.xyz`/`compute_dxa_real_case.msp`) — a single
closed component, 1648 triangles, 6408 Å² surface area, 0 open edges, 1 junction. My own mesh's
main component: 5060 triangles, **19529 Å²** — a **~3x** ratio on *both* metrics together, meaning
a genuinely wider dislocation-core tube, not just finer triangulation. Swept every tolerance that
could plausibly explain it (`vector_tolerance` 1e-6→1000: zero effect; `angle_tolerance` 5°→30°:
~3.5%; `rmsd_cutoff` 0.05→0.3: ~43%; `rcut` 3.6→8.0Å: ~37%) — none close the gap, even combined.

**Investigated two independent reference DXA implementations** (`/home/lafourcadep/CODES/ANALYSE/
lammps-disloc`, a 2014 CNA-triangulation method; `/home/lafourcadep/CODES/ANALYSE/DXA_SOURCE/
DXA1.3.6`, Stukowski's own 2010 CNA-based DXA predecessor — neither is literally the 2012
Delaunay/PTM paper this pipeline follows, but both independently converged on the same design
principle) via the Explore agent. Neither has anything resembling my per-edge "both endpoints'
fitted orientations must agree within an angle tolerance" gate. Both: classify each atom's
structure **independently** (no cross-atom check), propagate the ideal-lattice mapping
**one-sidedly** (stamp the matching atom's own orientation onto its own edges, no rejection at
assignment time), and **defer all defect detection to a later global circuit-closure test**. My
two-sided per-edge gate is far more sensitive to a dislocation's ordinary long-range elastic
lattice-rotation field than either reference design — exactly why it was overclassifying strained
atoms well outside the true core.

**Redesigned steps (iii)/(iv) to match**, in `src/ptm/compute_dxa_edge_vectors.cpp`:
- Step (iii): an edge is resolved using **only its lower-indexed endpoint's own frame** — the
  other endpoint's orientation is no longer consulted at all. `angle_tolerance` (default raised
  15°→40°) is now just a one-sided sanity bound (is the nearest ideal direction unambiguous?), not
  a cross-atom agreement requirement.
- Step (iv): dropped the edge-vector-closure test entirely. `DXATetClassification` is now simply
  "do all 4 vertices individually match `target_structure`" (`DXAEdgeVectors::
  vertex_matches_target`, a new per-vertex field) — matching both references' actual criterion.
  `vector_tolerance` slot removed (no longer meaningful).

**Result: 5324 → 3768 triangles** (~30% reduction) from these two changes alone, main component
area no longer separately re-measured but tracks the triangle-count drop closely given the
per-triangle size is roughly uniform. Still ~2.3x wider than OVITO's 1648 — the remaining gap
traces entirely to how many atoms PTM itself rejects (already-swept `rmsd_cutoff`/`rcut` plateau
well short), pointing at a *further* missing piece: `lammps-disloc`'s own
`optimizeAssignedLatticeVectors()` explicitly described "shrinking the region treated as
disordered as much as topologically possible" via majority-vote reconciliation *before* any
dislocation is identified — something with no counterpart anywhere in this pipeline.

## `ptm_shrink_disorder`: recovering strained-but-not-defective atoms

New operator, `src/ptm/ptm_shrink_disorder.cu`, run right after `compute_ptm` (per user's
request) and before `ptm_fields`. For each particle PTM rejected (`ptm_type == PTM_MATCH_NONE`):
tally neighbor `ptm_type` values (physical neighbors within `rcut`, same cutoff as `compute_ptm`
itself); if a clear majority (`min_neighbor_fraction`, default 0.5) already match some structure,
recover this particle to that same type. Repeated for `n_iterations` (default 3) so a
newly-recovered particle can help recover its own neighbors in a later pass — the recovered region
grows outward from confidently-matched bulk, same idea (not the same mechanism) as the reference
implementation's relaxation pass. A genuine dislocation core's own neighbors are themselves
disordered too, so there's structurally no majority to recover it with — this only removes
spurious individual rejections, not real defects.

**Technical note**: needed each neighbor's own *identity* (not just relative position, which is
all `compute_ptm.cu`'s own env-building ever needed) to look up its current `ptm_type` in the flat
buffer. Found `ComputePairBuffer2<false,true>` (flipping the 2nd template bool, `UseNeighbors`)
exposes `buf.nbh.get(j,cell,part)` per neighbor slot for free — `DefaultComputePairBufferAppendFunc`
already calls `tab.nbh.set(...)` unconditionally, it's a no-op only when `UseNeighbors=false`
(`compute_pair_buffer.h`). Confirmed via the Explore agent this is a real, already-exercised
code path (unlike the framework's *dynamic*-field neighbor-auto-packaging mechanism, which turned
out to have zero in-tree precedent for a `mk_generic_real` field like `ptm_type` — avoided that
riskier, untested route entirely).

**Verified**: only **281/176309 particles** (0.16%) recovered on the quadrupole case, plateauing
regardless of `n_iterations` (3+) or `min_neighbor_fraction` (≤0.5) — a genuine structural limit
of one-shell consensus, not a tuning gap. Despite the small recovery count, **interface triangles
dropped from 3768 → 2794** (a single recovered atom removes every "bad" tet touching it). Combined
with the redesign above: **5324 → 2794 total, ~48% reduction**, now ~1.7x OVITO's 1648 instead of
~3.2x. Sanity-checked on the mixed BCC/FCC/SC system (genuine phase boundaries, not noise): 1122
recovered there too, but interface triangles barely moved (5324→5254) — correctly, since that
"bad" region is genuinely non-BCC bulk, not strain noise, and majority-vote recovery has no reason
to touch it. No cross-phase contamination observed.

**Next steps, in order:**
1. The remaining ~1.7x gap is still open — likely needs either a deeper look at how OVITO's own
   PTM/structure-identification step differs (neighbor scheme, defaults), or a genuinely different
   shrink mechanism (e.g. iterating shrink+re-classify together, or a two-shell consensus).
2. Steps (viii)-(ix): sweep along the (now further-filtered) interface mesh to extract each
   dislocation's line geometry (a 1D curve through the tube of confirmed-bad tetrahedra), and
   detect junctions where multiple dislocations meet. `InterfaceMesh::edge_triangles`'s adjacency
   will need to grow into real traversal for this.
3. `DXABurgersCircuits::burgers_vector` gives a per-edge closure-failure vector already, but
   distinct dislocation lines nearby haven't been segmented/labeled yet — likely needed once line
   extraction exists, so each extracted line gets one clean Burgers vector rather than a cloud of
   per-edge residuals.

## `compute_cna`: OVITO's actual DXA classifier, tried as a `struct_field` swap-in — made things worse

OVITO's own DXA modifier documentation states it uses **CNA** (Common Neighbor Analysis), not PTM,
for local structure identification (`doc/manual/.../dislocation_analysis.rst`) — a correction of an
earlier wrong assumption in this investigation. Implemented from scratch in a new `src/cna/`
plugin (`compute_cna.cu`, `ComputeCNA` + `CNAFields` operators), following OVITO's own
adaptive-cutoff algorithm read directly from
`ovito/src/ovito/crystalanalysis/modifier/structureanalysis/StructureAnalysis.cpp`: per-atom
adaptive cutoff from the mean of the nearest 12 (FCC/HCP/ICO) or nearest-8-rescaled (BCC, 14
neighbors) bond lengths, scaled by `(1+sqrt(2))/2`; common-neighbor bond graph among the nn
candidates; tally the classic 4-2-1/4-2-2/5-5-5 (12-neighbor) / 4-4-4/6-6-6 (14-neighbor)
signatures; accept FCC/HCP/ICO/BCC on an exact match, no tolerance. Self-contained, no PTM
dependency. Same flat-buffer + `*_fields` materialization convention as `compute_ptm`/`ptm_fields`,
same `PTM_MATCH_*` numbering so it drops in as `compute_dxa_edge_vectors`' `struct_field`.

**Bug found and fixed**: initial implementation used union-find **component node-count** as the
third CNA signature index ("maxChainLength"). That's only correct for the *closed-ring* families
(444/555/666), where node count and edge count coincide. For the *open-chain* families (421: two
disjoint edges, 2 nodes/1 edge each; 422: one 2-edge/3-node path + 1 isolated node) they don't —
the correct metric is the largest **edge count** among components, not node count. This silently
misfiled every genuine FCC/HCP atom into the wrong bucket (421→"422", HCP's own 422→size-3 mismatch
too), giving **0/81751 matched** on a pure FCC lattice while BCC (444/666, both rings) worked fine
at 16000/16000 — the ring-only coincidence is exactly why BCC validated clean before this was
caught. Fixed by re-tallying edges per final union-find root in a second pass. Re-verified: BCC
16000/16000, FCC 32000/32000, HCP 4000/4000 (exercises both 421 and 422 together).

**Real-case comparison, corrected**: an earlier comparison here claimed CNA and PTM matched almost
identically (126976 vs 126977) and concluded classifier choice wasn't the gap's cause — that
comparison used `compute_ptm` at `rmsd_cutoff: 0.1` (`compare_cna_ptm_quadrupole.msp`'s own
setting), not the actual pipeline's `0.2` (`compute_dxa_real_case.msp`). At the pipeline's real
tolerance, PTM matches **127719/176309** (only 281 owned atoms rejected) vs CNA's **126976**
(1024 rejected) — a **3.6x** difference. The user supplied OVITO's own real CNA output on the
identical file (`data/regression_new/delaunay/ovitodata/output_cna_ovito.xyz`): **1026** non-BCC
atoms — matching this codebase's own (now bug-fixed) `compute_cna` almost exactly (1024 vs 1026).
So the original "3x-wider-than-OVITO" mystery was, in large part, simply PTM's continuous RMSD fit
being far more strain-tolerant than any real discrete crystal-structure classifier — `compute_ptm`
was finding a defect core roughly 3.6x too small. `struct_field: cna_type` is the fix; `struct_field
choice alone isn't the lever` (this file's own prior conclusion) was wrong, see below.

## `compute_atomistic_interface_mesh`: the real DXA mesh construction, in a new operator

Per-user request, and per the OVITO-documentation dead end above (`## Surface mesh...`, further
down this file if present, or see git history): the *interface mesh itself*, not just the
classifier, needed a from-scratch reimplementation, matching what
`DXA_SOURCE/DXA1.3.6/src/analysis/interfacemesh/InterfaceMesh.cpp` actually does — mesh vertices
are the non-crystalline atoms themselves (real positions, not tet-derived points), and facets come
from each *crystalline* atom's own fixed local BCC lattice template (6 quads of first/second-shell
neighbor slots), never an independent tessellation of every atom. Implemented as a new,
self-contained operator, `src/delaunay/compute_atomistic_interface_mesh.cpp`
(`compute_atomistic_interface_mesh`), deliberately kept separate from `compute_interface_mesh.cpp`
(both stay useful/reusable, per the user's own framing — the tet-boundary approach isn't obsoleted).
Produces the same `InterfaceMesh` struct, so `write_interface_mesh` works unchanged. BCC only for
now (FCC/HCP need the "8 Thompson tetrahedra" template instead of BCC's 6 quads — same technique,
tables already extracted from DXA1.3.6 during this investigation, just not wired in yet).

**Lattice templates**: DXA1.3.6's own hardcoded numeric tables (`src/lattice/LatticeTypeBCC.cpp`)
were extracted and transcribed verbatim (8 first-shell `<111>`-type + 6 second-shell `<100>`-type
directions, the 6 quads' index tables) rather than re-derived — safer than trying to regenerate
this geometry from scratch. Per-atom neighbor-to-canonical-slot resolution reuses
`compute_dxa_edge_vectors`' own "rotate the real bond into the atom's `ptm_orientation` frame, snap
to nearest ideal direction" technique, but computed **independently for both endpoints of every
edge** (that operator only resolves the lower-indexed endpoint, by design — this operator needs
every crystalline atom's own full 14-slot map, not half of them), reusing its already-deduplicated
edge list and `vertex_matches_target` flags as input.

**Bug-hunting story, in order**:
1. First real run (PTM-driven `vertex_matches_target`, same as the rest of the pipeline): only
   **91-93 triangles** — dramatically short of OVITO's 1648. Instrumented rather than guessed:
   confirmed the slot-resolution rate was fine (~96%) and the disordered population was genuinely
   just **186 unique atoms** (matching a single-atom-wide dislocation-core-length estimate) — i.e.
   not a bug in the resolution math, but too few "holes" for the per-atom quad mechanism to have
   anything to connect (it needs *simultaneously* non-crystalline neighbors in specific relative
   slots — a lone atom in an otherwise-crystalline neighborhood satisfies none of BCC's quad
   conditions). Cross-checked against the reference source itself:
   `DXAInterfaceMesh::createBCCMeshEdges()` is a **literal empty function** in DXA1.3.6 — BCC gets
   no independent edge pre-population the way FCC/HCP do, so even DXA1.3.6's own hole-closing pass
   (which only walks *existing* edges) has nothing to grab onto around a sparse BCC core either.
   This ruled out "just port more of closeFacetHoles" as the primary fix.
2. Added a bounded hole-closing pass anyway (walks the open-edge subgraph's clean simple loops —
   every loop vertex has exactly one required outgoing direction, derived from the opposite of
   whatever direction the existing bordering triangle already winds that edge — and fan-
   triangulates each loop in its own traversal order, which preserves the induced orientation
   automatically). Branch points (3+ open edges) or genuine dangling edges are left open, not
   forced. Went 71→93 triangles on the PTM-driven run — real, but small next to the scale of the
   gap, confirming (1) wasn't primarily a "few residual gaps" problem.
3. Root cause was the *classifier feeding it*, not the mesh code (see the corrected comparison
   above): switched `dxa_edge_vectors`' `struct_field` to `cna_type`. Non-crystalline population
   jumped from 186 to **1026 unique atoms** — an exact match with OVITO's own CNA output on the
   identical file. Re-ran the atomistic mesh: **2020 triangles**, only 34 open edges remaining
   (out of 924948 total tessellation edges) — down from ~3x off to **~22% over** OVITO's 1648.

**Current best pipeline** (`data/regression_new/delaunay/compute_dxa_real_case_cna_atomistic.msp`):
`compute_ptm` (still needed for `ptm_orientation` — `compute_cna` doesn't produce per-atom
orientation) → `compute_cna`/`cna_fields` (classification) → `ghost_update_opt` on both
`cna_type`+`ptm_orientation` → `compute_delaunay` → `compute_dxa_edge_vectors` with
`struct_field: cna_type` → `compute_atomistic_interface_mesh`.

**Next steps, in order:**
1. Close the remaining ~22% gap (2020 vs 1648): likely candidates are the 34 still-open edges
   (branch points this operator's simple hole-closer can't resolve — DXA1.3.6's full bounded
   backtracking search is the natural upgrade), or minor edge-resolution/`angle_tolerance`
   differences from OVITO's own exact pipeline.
2. FCC/HCP support (Thompson-tetrahedra templates, already extracted from DXA1.3.6, not yet wired
   into `compute_atomistic_interface_mesh`) if a non-BCC target_structure is ever needed here.
3. Steps (viii)-(ix) (dislocation line extraction, junction detection) — unchanged from before,
   now more promising to pursue against this operator's much-closer-to-OVITO mesh.

Test files: `data/regression_new/delaunay/compute_cna_test.msp` (pure-lattice CNA validation),
`compute_dxa_real_case_atomistic.msp` (atomistic mesh, PTM-driven classification — the 91-triangle
under-coverage case), `compute_dxa_real_case_cna_atomistic.msp` (atomistic mesh, CNA-driven — the
2020-triangle result), `ovitodata/output_cna_ovito.xyz` + `ovitodata/output_dxa_ovito.vtk`
(user-supplied OVITO ground truth: 1026 non-BCC atoms, 1648-triangle reference mesh).
