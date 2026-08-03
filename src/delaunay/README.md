# Delaunay / PTM / DXA pipeline

Goal: build towards a DXA-style dislocation extraction pipeline in exaStamp, following
Stukowski, Bulatov, Arsenlis, *"Automated identification and indexing of dislocations in
crystal interfaces"*, Model. Simul. Mater. Sci. Eng. 20 (2012) 085007
(local copy: `/home/lafourcadep/Bureau/DXA_ref.pdf`).

**Regression data layout (split 2026-08-03)**: `data/regression_new/delaunay/` now holds only the
raw Delaunay-tessellation test (`compute_delaunay.msp`, MPI/ghost tet-count correctness -- no DXA
operator involved at all); every actual DXA test (`.msp`, `.xyz` inputs, OVITO ground truth under
`ovitodata/`, outputs under `paraview/`) moved to the sibling `data/regression_new/dxa/` folder.

DXA's published pipeline has 9 steps. Status:

| Step | What it does | Status |
|---|---|---|
| (i) | Atomic structure identification: crystal type **and** local lattice orientation per atom | crystal type: **use `compute_cna`** (`src/cna/`), not `compute_ptm` -- PTM's continuous RMSD fit under-classifies real defects (see "`compute_cna`" section below). Orientation: **done** via PTM's `ptm_orientation` (`compute_cna` doesn't produce one) |
| (ii) | Space-filling Delaunay tessellation | **done** — this directory |
| (iii) | Assign an ideal lattice vector to each tessellation edge | **done** — `compute_dxa_edge_vectors`, `src/ptm/`, see below (use `struct_field: cna_type`) |
| (iv) | Classify each tetrahedron good/bad from edge compatibility | **done, matches OVITO's documented criterion** — `compute_dxa_tet_classification`, `src/ptm/`: bad if any of its 6 edges is unresolved (fixed from an earlier vertex-matching test that undershot the docs, though empirically identical at this system's `angle_tolerance=40`, see below) |
| (v) | Build the interface mesh (2D manifold separating good/bad regions) | **done, two alternative implementations, open question on which is right** — `compute_interface_mesh` (tet-classification-boundary, matches OVITO's own documented construction) and `compute_atomistic_interface_mesh` (DXA1.3.6's *older* atom-centric construction, empirically closer to OVITO's actual triangle count) — see "Remaining work" item 1 |
| (vi)-(vii) | Burgers circuit construction, real-vs-noise defect filtering | **corrected mid-session** — circuits must be built on the interface mesh's own edges (OVITO's documented algorithm), not a whole-crystal graph. `compute_dxa_mesh_burgers_circuits` (new, `src/delaunay/`) does this correctly on the atomistic mesh; `compute_dxa_burgers_circuits` (older, whole-crystal graph) is kept for the tet-boundary mesh's own bad-tet-reclassification use. See below |
| (vi)-(ix) | Burgers circuits + line/junction extraction | **rewritten to match OVITO's own real (non-public) source, `DislocationTracer.cpp`** — `compute_dxa_circuit_sweep` now does territorial exclusion in the seed search itself, lockstep incremental growth over increasing trial-circuit length, and incremental (not post-hoc) 2-arm merging. Validated against exact `.ca` ground truth from the user: **screw dipole exact 2/2** (no merging even needed); **quadrupole ~19-20 vs. true 11**, but 100% precision (every fragment maps onto a real line, zero false positives) — pure over-fragmentation, not noise. **One open bug**: circuits can drift across a real 3+-way junction into a neighboring dislocation's core during ordinary growth (not a merge-decision issue — reproduced with all merging disabled) — see item 3 below for the full diagnosis and a promising lead (OVITO's `CrystalPathFinder`, item 4). `compute_dxa_dislocation_lines`/`compute_dxa_mesh_dislocation_lines` (earlier line-tracing attempts) remain superseded; `compute_dxa_mesh_burgers_circuits`/`compute_dxa_burgers_circuits` remain independently valid for other uses. |

## Remaining work (checked here before resuming — see linked sections below for full detail)

0. **IN PROGRESS: full re-architecture to match OVITO's real, current DXA source exactly** (user
   request, checked directly against `ovito/src/ovito/crystalanalysis/modifier/dxa/*.{h,cpp}` --
   the actual 2025 source, not the old DXA1.3.6 predecessor or the public manual). This traces the
   drift bug and interface-mesh gap to something deeper than item 4's original CrystalPathFinder
   lead: OVITO's real elastic mapping doesn't use a PTM-style continuous orientation fit at all --
   `StructureAnalysis::determineLocalStructure` does a DISCRETE graph-topology match (CNA signature
   + neighbor-bond-graph isomorphism) against fixed reference tables, and
   `InterfaceMesh::createMesh`'s good/bad criterion is `ElasticMapping::isElasticMappingCompatible`
   (a genuine per-tetrahedron Burgers-circuit-closure + Frank-rotation test), not our current
   per-vertex/per-edge-count heuristics. Full plan (4 stages, user confirmed "full port, in order"):
   1. **DONE, validated** -- discrete neighbor-slot classifier + cluster graph, replacing PTM's role
      entirely for this pipeline. New files: `lattice_structure.h/.cpp` (verbatim BCC/FCC/HCP
      tables + point-group symmetry-permutation search, unit-tested standalone: FCC/BCC both find
      the correct 48-element Oh group, HCP the correct 12-element local group), `cluster_graph.h/
      .cpp` (Cluster/ClusterTransition/ClusterGraph, unit-tested standalone), `compute_dxa_lattice_
      correspondence.cpp` (per-atom backtracking permutation match --
      `StructureAnalysis::determineLocalStructure` port), `dxa_lattice_clusters_algo.cpp` +
      `compute_dxa_lattice_clusters.cpp` (BFS cluster growth + inter-cluster transitions --
      `buildClusters`/`connectClusters` port). **Verified on the real installed build**
      (`compute_dxa_lattice_clusters.msp`, same unrotated 16000-atom BCC Ta lattice `compute_ptm`
      was originally validated on): **16000/16000 owned particles matched, exactly 1 cluster, 0
      transitions** -- matches physical expectation exactly for a single perfect grain.
   2. **DONE, validated** -- `CrystalPathFinder` port. New files: `dxa_crystal_path.h` (header),
      `dxa_crystal_path_finder.cpp` (grid-independent `dxa_crystal_path_find`, unit-tested
      standalone on the same synthetic BCC lattice as stage 1: **756/756 direct-neighbor pairs and
      11/11 multi-hop, non-directly-bonded pairs exactly reproduce the real spatial vector** once
      transformed back through the cluster's own orientation fit), `compute_dxa_crystal_path_edge_
      vectors.cpp` (grid operator, `assignIdealVectorsToEdges` port). **Along the way, found and
      fixed a real integration gap**: `compute_dxa_lattice_correspondence` originally only
      classified owned cells (matching compute_ptm/compute_cna's own convention), leaving every
      ghost-copy atom unclassified -- since Delaunay tessellation vertices routinely land on ghost
      particles, this left ~7% of edges unresolved purely from ghost atoms never having a
      correspondence at all (measured: 104861/113012, 92.8%, on the plain BCC Ta system). Widened
      the classifier to cover the full grid including ghost cells (no PTM/CNA-style named-field
      `ghost_update_opt` sync needed, since `DXALatticeCorrespondence` isn't a named field anyway,
      and a ghost is its own real atom with its own real neighborhood -- classifying it directly is
      both simpler and correct for this project's single-MPI-rank scope, where every "ghost" is a
      periodic self-image with a fully populated local environment). Fixed: **113012/113012
      (100%)** edges resolved on the plain BCC Ta system. Documented caveat in the operator's own
      comment: a genuine multi-rank cross-rank ghost near the outer halo edge could still see a
      truncated neighbor count, same "ghost-fringe-trust" concern `compute_delaunay.cpp` already
      has for tets -- not yet exercised multi-rank. **Verified on the real quadrupole dislocation
      case** too (`compute_dxa_crystal_path_real_case.msp`): 152571/176309 matched (owned+ghost),
      **3 clusters, 0 transitions** (plausible: a real dislocation network can locally separate
      otherwise-good BCC regions from each other with no directly-bonded matched-atom path between
      them), **914276/924948 (98.85%) edges resolved** -- correctly less than 100% now that a real
      defect network exists, unlike the perfect-lattice case.
   3. **DONE, validated** -- `isElasticMappingCompatible` port (`dxa_elastic_mapping_compatible.cpp`
      + `compute_dxa_elastic_mapping_tet_classification.cpp`): a genuine per-tetrahedron
      Burgers-circuit-closure + Frank-rotation test on the new `DXACrystalPathEdgeVectors` data,
      replacing the old per-vertex/edge-resolution-count heuristics entirely. **On the plain,
      perfect BCC Ta system: 91200/91200 (100%) tetrahedra good** -- exactly right, zero defects.
      **On the real quadrupole dislocation case: 755770/768041 (98.4%) good, 12271 bad (1.6%)** -- a
      small, tight defect-core fraction, qualitatively much tighter than either old heuristic ever
      achieved (tet-boundary mesh: 5164 triangles off ~172559 tets, ~3%; atomistic mesh: 2020
      triangles vs OVITO's real 1648) -- though not yet directly comparable in the same units until
      stage 4 actually builds an interface mesh from this classification and gives a triangle count
      to compare against OVITO's 1648 reference.
   4. **DONE, validated -- and this is the actual fix for the open drift/over-fragmentation bug.**
      New file `compute_dxa_elastic_interface_mesh.cpp`: same tet-boundary triangle-extraction
      algorithm as `compute_interface_mesh.cpp`, fed by the new elastic-mapping classification
      (stage 3) instead of the old edge-resolution-count criterion; edge ideal vectors come directly
      from `DXACrystalPathEdgeVectors` (every interface-mesh edge IS a real tessellation edge here,
      unlike the atomistic mesh, so no re-derivation needed). Produces the *same* `InterfaceMesh`
      struct `compute_dxa_circuit_sweep` already consumes -- **the sweep itself needed zero code
      changes**, confirming the earlier analysis that OVITO's own `DislocationTracer` never
      re-checks elastic-mapping consistency during sweep; the fix has to come from (and only from)
      the interface mesh being built correctly in the first place.

      **Real quadrupole case, no additional tuning at all**: interface mesh has **5164 triangles**
      (vs. OVITO's 1648 -- more on this below) but the metric that actually matters, dislocation
      count after the sweep, is **18 physical dislocations vs. the true 11** -- beating this
      project's previous best (19-20, achieved only after extensive move-set tuning on the
      atomistic mesh) on the very first run of the new architecture, with clean stop reasons (0
      max-length, 0 self-closure, 78 junction, 0 open-edge, 0 exhausted). **Screw dipole case
      (regression check)**: still **exact 2/2**, lengths 110.0/110.4 Å, matching the previous best
      result exactly (4 open-mesh-edge stops, same known non-periodic-image-summed synthetic-field
      boundary artifact as before).

      **Honest open finding, not yet resolved**: the interface mesh's own triangle count (5164) is
      *higher* than both OVITO's reference (1648, ~3.1x) and this project's own current live
      edge-resolution-count pipeline (2794, ~1.7x) -- so by triangle-count alone the new, more
      principled classifier currently looks numerically worse, even though it produces a
      substantially better final dislocation count. Two candidate causes identified, neither yet
      confirmed: (a) `CA_LATTICE_VECTOR_EPSILON=1e-3`/`CA_TRANSITION_MATRIX_EPSILON=1e-4` (OVITO's
      own published values) may be tighter than this dataset's real MD relaxation noise tolerates
      -- back-of-envelope: 0.01 Å thermal noise at this system's lattice constant is already ~3x
      the 1e-3 lattice-unit epsilon: not yet swept/measured. (b) tried loosening
      `crystal_path_steps` 4->8: **zero effect** (edge-resolved count identical bit-for-bit),
      ruling out path-length as the lever -- the ~1.15% unresolved edges are dominated by genuinely
      unclassified-atom endpoints (`cluster==0` short-circuit before any path search even runs), not
      a too-short search radius. Not yet investigated further given the strong downstream result
      already achieved; worth an epsilon sweep next if the triangle-count gap itself becomes
      important (e.g. for a defect-mesh visualization matching OVITO's own more tightly, as opposed
      to just correct dislocation-line topology).

      **Follow-up investigation, same day: epsilon and ghost-margin both ruled out; classifier
      itself confirmed exactly correct.** Swept `CA_LATTICE_VECTOR_EPSILON`/
      `CA_TRANSITION_MATRIX_EPSILON` over a 9-point grid (1e-3 to 5e-2, ~50x range): **bit-identical
      result every time** (755770/768041 good, 5164 triangles) -- epsilon has zero effect, ruled
      out definitively. Added a diagnostic counting atoms that pass the aggregate CNA count but fail
      the exact bond-topology permutation match: **0**, on the real quadrupole case -- the discrete
      classifier's extra graph-isomorphism stage never rejects anything CNA's own aggregate count
      wouldn't already reject. Then directly compared owned-only match counts (compute_cna never
      touches ghost cells, so its own raw "matched/total" is misleadingly denominated over
      owned+ghost): **this classifier's owned-only match is 126976/128000, bit-identical to
      compute_cna's own owned-only match on the same file** -- full confirmation the classifier
      itself is exactly faithful, not the source of the gap. Tried widening the ghost halo margin
      (`rcut_max` 6.0->10.0 ang) in case ghost atoms near the halo's outer edge were getting
      spuriously rejected from a truncated neighbor search: raised owned+ghost matched count
      (152571->175285) but **the downstream tet-classification and triangle count were completely
      unchanged** (bit-identical 755770/768041, 5164) -- ruling out ghost margin too, since
      `compute_delaunay`'s own centroid-in-owned-cell trust rule already restricts which tets (and
      therefore which ghost vertices) matter, and none of the newly-recovered far-ghost atoms were
      referenced by any kept tet anyway.

      **Follow-up, next session: the "smoothed defect mesh" hypothesis above was WRONG, disproven
      directly.** User pointed out a real, runnable OVITO Pro Python interpreter exists locally
      (`/home/lafourcadep/CODES/VISU/ovito-pro-3.14.1-x86_64/bin/ovitos`) -- previously assumed
      unavailable (OVITO's own source-only checkout has an unbuilt `build/` dir, no compiled
      binary). This is a genuinely reusable capability going forward: OVITO ground truth no longer
      has to be user-supplied by hand, it can be generated directly for any test case via `ovitos` +
      `DislocationAnalysisModifier` + `export_file(..., "vtk/trimesh", key="dxa-interface-mesh")`
      (the defect mesh is a separate object, `key="dxa-defect-mesh"`).

      Re-generated the quadrupole's own interface mesh fresh via `ovitos`: **1648 triangles, 6400.9
      Å², identical line count to the user-supplied `output_dxa_ovito.vtk`** -- confirms that file
      IS the raw interface mesh (not the defect mesh, which for this same system exports as a
      degenerate 12-triangle mesh, clearly a different/oddly-behaved object, not what was being
      compared against all along). The "different metric" hypothesis is dead; the 3.1x triangle/area
      gap on the quadrupole is a real, apples-to-apples comparison.

      **But then generated the screw dipole's own OVITO interface mesh for the first time ever**
      (no such reference existed before this session) and got a real surprise: **OVITO's own screw-
      dipole interface mesh is 3004 triangles, 9479 Å²** -- our own screw-dipole mesh (3462
      triangles, 12617 Å²) is only **1.15x the triangles, 1.33x the area** of OVITO's real
      reference -- nowhere near the quadrupole's 3.1x gap. (A same-session self-consistency estimate
      using area/length ratios, made before this real ground truth existed, wrongly concluded the
      screw dipole was proportionally *worse* than the quadrupole -- retracted; that heuristic's
      implicit assumption, that BCC Ta's tube width-per-length should be similar between the two
      test systems, turns out false even in OVITO's own real output: OVITO's own screw-dipole
      area/length ratio is ~41 Å vs. its own quadrupole's ~10 Å, i.e. OVITO's real tube is *itself*
      proportionally much fatter on this specific (unrelaxed, synthetic-displacement-field) dataset
      than on the real MD-relaxed quadrupole.) **Conclusion: the ~3x-wider-region gap is concentrated
      in the quadrupole's junction-dense, real-relaxed-MD regime specifically, not a generic
      property of the elastic-mapping classifier** -- a much more localized, actionable lead than
      "epsilon" or "smoothed-mesh-mismatch" ever were. Next step if pursued: compare where in the
      quadrupole's own mesh the extra area concentrates (near real junctions vs. along ordinary line
      segments) to isolate the cause further.

      Fresh OVITO references saved for reuse: `ovitodata/output_dxa_screw_dipole_interface_ovito.vtk`
      (3004 triangles, the new ground truth), `ovitodata/output_dxa_screw_dipole_defect_ovito.vtk`,
      `ovitodata/output_dxa_quadrupole_defect_ovito.vtk` (the degenerate 12-triangle defect mesh,
      kept for reference even though it turned out not to be what `output_dxa_ovito.vtk` was).

      **Follow-up, same session: found and fixed the actual root cause -- a whole missing port
      step, `ElasticMapping::assignVerticesToClusters()`.** User asked to keep pushing for an exact
      match. Systematically re-tested and ruled out, each with a clean measurement: (a) a
      cluster-graph-transition mechanism (only 11/924948 edges affected, negligible); (b) the
      `DXA_MAX_NEIGHBORS` append cap silently dropping neighbors (0 drops measured); (c)
      `StructureAnalysis::formSuperClusters()` -- confirmed via `grep` that OVITO's real DXA
      pipeline (`DislocationAnalysisEngine.cpp`) never even calls it (only a *different* modifier,
      Elastic Strain, does) -- a real, clean dead end, not a bug; (d) classifier fidelity at the
      SET level (not just count): dumped both our own and OVITO's own per-atom rejected-atom id
      lists via `ovitos` and diffed them directly -- **our rejected set is a perfect subset of
      OVITO's, differing by exactly 2 atoms out of 1026** (OVITO rejects 2 extra borderline atoms
      we accept) -- as close to "exact" as classification gets, and far too small to explain a 3x
      gap; (e) tessellation density -- cross-checked our own Geogram-based tessellation's tets/edge-
      per-atom ratio against an independent SciPy/Qhull triangulation of the identical raw point
      cloud: 6.00/7.23 (ours) vs 6.05/7.05 (Qhull) tets,edges per atom -- normal, unremarkable,
      ruling out "our tessellation is unusually dense" too.

      **The real cause, found by re-reading `ElasticMapping::assignIdealVectorsToEdges` once more
      and noticing it calls `clusterOfVertex()`, not `structureAnalysis().atomCluster()`, for its
      own gate check** (`if(cluster1->id==0 || cluster2->id==0) continue;`) -- `clusterOfVertex()`
      returns a value from a SEPARATE, previously-unported method,
      `ElasticMapping::assignVerticesToClusters()`, which propagates a cluster id to **every**
      tessellation vertex (not just classified atoms) by flood-filling outward through ordinary
      tessellation-edge adjacency from already-clustered vertices -- entirely distinct from
      `StructureAnalysis::buildClusters`/`connectClusters` (both purely atom-classification-level,
      restricted to each atom's own *native* 14-neighbor list). Our port only ever had the
      atom-classification-level cluster (0 for any unclassified atom, permanently), and used THAT
      raw value for the edge-resolution gate -- meaning every tessellation edge touching *any* of
      the ~1024 unclassified atoms was rejected outright before `CrystalPathFinder` ever got a
      chance to route around it via its own reverse-neighbor-search mechanism (which was ported
      correctly and does still use the raw, un-propagated cluster internally, exactly matching
      `CrystalPathFinder::findPath`'s own `structureAnalysis().atomCluster()` call --
      only the *outer gate check* was using the wrong cluster source). New function
      `dxa_propagate_vertex_clusters()` (`dxa_crystal_path_finder.cpp`, BFS flood-fill over the full
      tessellation-edge graph -- OVITO's own "repeat until no change" scan converges to an
      equivalent result, just less efficiently) now feeds the gate check, the re-expression target,
      and the edge's own recorded cluster transition in `compute_dxa_crystal_path_edge_vectors.cpp`
      -- exactly mirroring which of the two cluster sources OVITO's own code uses at each specific
      point.

      **Result, quadrupole (no other changes)**: edge resolution 98.85% -> **99.98%** (924768/
      924948); tets good 98.4% -> **99.87%** (767017/768041, bad tets 12271 -> 1024); **interface
      mesh 5164 -> 1692 triangles vs. OVITO's 1648 -- a 2.7% difference, down from 3.13x.**
      Downstream circuit sweep (no code changes there either, confirming again the fix belongs
      entirely at the mesh-construction level): **8 physical dislocations vs. the true 11** (was
      18-21), total length 612.62 Å vs. OVITO's real coarsened 629.16 Å (2.6% off). Per-line
      breakdown against ground truth is very clean: our 3 shortest lines (6.56, 7.23, 13.25 Å)
      closely match OVITO's own 3 short junction-type lines (9.05, 9.4, 15.08 Å); 3 of our
      "classical" lines land right on 3 of OVITO's 8 (68.29≈68.24, 72.69/72.77≈72.24/73.68); and our
      remaining two long lines (152.93, 218.9 Å) sum to 371.83 Å, matching the sum of OVITO's
      remaining 5 classical lines (381.47 Å) within 2.5% -- **the sweep is finding the same 11
      physical dislocations, but the incremental two-arm merge logic is occasionally still
      over-merging 2-3 real, distinct dislocations that meet at one junction into a single long
      segment.** This is now a narrow, well-characterized remaining gap in `compute_dxa_circuit_
      sweep`'s own merge heuristic (likely `MERGE_GRACE_ROUNDS` needing retuning now that more edges
      resolve and arm growth timing has changed), not a classification, tessellation, or
      elastic-mapping problem -- those are now effectively solved.

      **Screw dipole regression, more nuanced**: interface mesh dropped 3462 -> 941 triangles
      (OVITO's own reference: 3004) -- now *under*, not over. But the actual dislocation output is
      completely unaffected (still exact 2/2, lengths ~110/110 Å, matching before and OVITO's own
      228.8 Å reference closely) -- the lost mesh area doesn't touch the real dislocation cores.
      Likely explanation, not fully confirmed: this specific synthetic test file has an
      already-documented confound (no periodic-image summing in its own construction, leaving a
      genuine spurious strained ribbon at the domain boundary, `239 cutoff edges` here) -- before
      this fix, part of that ribbon was included as "bad" purely via the missing-edge shortcut
      (not real physics); now that real edges resolve there and get a fair Frank/Burgers test (240
      genuine failures now, vs. 0 before), much of that region correctly comes back "good" instead
      of being auto-flagged bad. Not chased further this session -- the quadrupole (real, MD-relaxed
      ground truth) is the reliable signal, and it improved dramatically.

      **Follow-up, same session: found and fixed a second real bug, this time in
      `compute_dxa_circuit_sweep`'s own merge logic**, by re-reading OVITO's real
      `DislocationTracer::joinSegments()` + `DislocationNode::connectNodes()`/`formsJunctionWith()`
      in detail. OVITO's real 2-way-vs-3+-way junction decision is a proper ring-union structure
      (`junctionRing`, a circular linked list): while scanning a stopped circuit's *entire*
      boundary, it calls `connectNodes()` for *every* distinct adjacent circuit it touches,
      naturally building a ring whose size (`countJunctionArms()`) directly says "how many circuits
      meet here" -- `armCount>=3` is a real junction (kept separate), `armCount==2` merges. Our own
      port only ever recorded a *single* `blocking_node` per stopped node -- when a circuit's
      boundary directly touched 2+ *different* other circuits at once (exactly the local signature
      of a real 3+-way junction), whichever was encountered *last* while scanning silently
      overwrote the earlier one, discarding the direct evidence and falling back on a *global*
      incoming-count proxy that isn't equivalent. Fixed: `NodeState::blocking_node` ->
      `blocking_nodes` (collects every distinct touched node, not just the last), plus a new
      `resolved_distinct_blockers()` helper (resolves each through any prior merge chain and dedups)
      used everywhere a merge decision is made -- if a node's own resolved set has size != 1, that
      alone is now definitive, immediate, local evidence of a real 3+-way junction (no need to even
      wait on the incoming-count/grace-period checks).

      **Result on the quadrupole**: raw segments 12 -> 9, absorbed-via-merge 4 -> 1 -- a real,
      measured reduction in improper merges (matches the mechanism fix directly). Final dislocation
      count stayed at **8** (unchanged) and total length improved slightly (612.6 -> 615.5 Å, now
      2.2% off OVITO's 629.2 Å, down from 2.6%) -- the fix is verified correct and real, but the
      remaining 8-vs-11 gap has shifted: it's no longer primarily an incorrect-merge problem (that
      mechanism is now much more locally principled, matching OVITO's own), it looks more like a
      **raw seed-discovery gap** -- only 9 independent segments get seeded in the first place for a
      system with 11 real dislocations + 7 real junctions, before any merging even happens. Not yet
      investigated; the natural next place to look is `try_seed_from`'s own local trial-circuit
      search (why doesn't every real arm get its own independent seed before growth starts
      colliding with a neighbor's territory) rather than anything in the merge logic itself.

      **Follow-up, same session: diagnosed the seed-discovery gap precisely -- it's an
      order-dependent territorial race in `try_seed_from`'s sequential vertex scan, confirmed but
      not yet fixed.** Compared against OVITO's real `DislocationTracer::findPrimarySegments()`:
      structurally very similar (same "stop at first valid closing edge" BFS, same territorial
      exclusion checks) -- one real, but currently inert, faithfulness gap found: OVITO's search
      also verifies the accumulated Frank-rotation matrix agrees between the two BFS paths meeting
      at a candidate closing edge (`frankRotation.equals(neighborStruct->tm, ...)`), not just that
      the Burgers vector is nonzero; our own `try_seed_from` never checks this. Doesn't currently
      matter for either test system (both are single-cluster, so every transition is trivially the
      identity) but would matter for a genuine multi-grain system -- worth porting for full fidelity
      even though it's not the cause of the current gap. **The real, confirmed cause**: reversed the
      seed-scan order (`for(root=n_vertices-1; ...; root--)` instead of ascending) as a diagnostic,
      with zero other changes -- **result jumped from 8 to 10 physical dislocations** (vs. the true
      11), with all 3 short junction-type lines now found (previously only 2) and only one remaining
      compound-merged line (vs. two before). This conclusively confirms the gap is a genuine
      territorial race: whichever of two nearby real dislocation arms gets tried as a seed *first*
      (by raw vertex index) claims territory the other needs, and the second arm never gets its own
      independent seed at all -- it just silently disappears rather than erroring. **Not adopted as
      a fix**: reversing the scan order is itself just as arbitrary as the original ascending order
      -- it happens to do better on this one test case, but adopting it outright would be curve-
      fitting to the only ground truth available, not a real solution. A genuine fix would need a
      properly order-*independent* seed strategy (e.g. a "most locally-constrained region first"
      priority, or a truly simultaneous/parallel seeding pass) -- a real design effort, not
      attempted this session. Also re-tested `create_secondary_segment` (disabled since an earlier
      session, when it measurably made results worse) now that the real Frank-rotation-based
      elastic mapping exists: **still makes things worse** (22 dislocations, several exactly
      degenerate/zero-length) -- re-confirms the issue is specifically the hole-closing loop's own
      missing consistency check (not the general elastic-mapping fidelity, which is now good), left
      disabled.

      **Follow-up, same session: found and fixed a real, narrow classification bug behind the
      remaining gap, per user's "if the seed strategy is the same, why is our result different?"
      question.** Rather than accept the seed-ordering race as an unavoidable design limit, checked
      whether the *input* to the sweep (the interface mesh, and the classification feeding it) was
      truly identical to OVITO's, atom for atom, not just count-for-count. Clustered the 11
      dislocations' own endpoints from OVITO's real `.ca` ground truth into 7 junction positions
      (matching the known 6x3-way + 1x4-way topology exactly), then checked the 2 atoms where our
      classification still differed from OVITO's own (found earlier via the direct id-set diff, see
      item 0 above): **both sit within ~3.5 Å of the exact same real junction**, and straddle it
      almost symmetrically. Reproduced the BCC classification test independently in Python for these
      2 atoms and confirmed they correctly FAIL the BCC 14-neighbor sanity gate (matches OVITO) --
      but then found they also happen to **exactly satisfy the FCC 12-neighbor combinatorial test**
      (n421=12, a perfect topological match) purely from local strain distortion, even though the
      whole system is pure BCC. Root cause: `compute_dxa_lattice_correspondence` tested FCC/HCP
      *then* BCC unconditionally for every atom, regardless of what the pipeline actually wanted --
      unlike OVITO's real `DislocationAnalysisModifier`, whose `input_crystal_structure` is a single
      required choice that's the *only* structure ever tested per atom. Fixed: new required
      `target_structure` slot (`"BCC"`/`"FCC"`/`"HCP"`, default `"BCC"`), gating which family the
      functor even attempts -- exactly mirroring OVITO's own single-target design.

      **Result: classification now 126974/128000 -- an EXACT match to OVITO's own real
      DXA-internal count** (previously 126976, a 2-atom mismatch, now closed to zero). Interface
      mesh **1656 triangles vs. OVITO's 1648 -- 0.5% off**, down from 2.5-2.7% before this fix.
      Circuit sweep: **9 physical dislocations vs. the true 11** (up from 8), lengths [7.14, 7.19,
      12.63, 66.77, 67.34, 72.97, 73.79, 153.29, 153.58], total 614.71 Å vs. OVITO's 629.16 Å
      (2.3% off) -- the 3 short junction-type lines are now found essentially exactly (7.14/7.19/
      12.63 vs. OVITO's 9.05/9.4/15.08), and only 2 of OVITO's 8 classical lines remain compound-
      merged (down from a messier split before). Screw dipole regression: completely unaffected
      (still exact 2/2, mesh still 941 triangles) -- expected, since that test never had a spurious
      FCC-classified atom to begin with (single-cluster, no junction network). **This closes the
      classification-fidelity gap entirely** -- the residual 9-vs-11 count is now attributable
      solely to the still-open, still-unfixed order-dependent seed-race documented in the paragraph
      above, not to any remaining classification or mesh-construction discrepancy.

   5. **DONE, validated -- line coarsening + smoothing** (user request, matching OVITO's own
      `linePointInterval`/`lineSmoothingLevel` mechanism exactly). New files:
      `smooth_dxa_dislocation_lines.cpp` (`coarsen_dislocation_line`/`smooth_dislocation_line`, ported
      from `DislocationNetwork::coarsenDislocationLine()`/`smoothDislocationLine()`), plus a new
      `DXADislocationLines::core_size` field (parallel to `line_positions`, the sweeping circuit's own
      loop size when each point was recorded -- OVITO's own `DislocationSegment::coreSize`, used to
      weight the coarsening merge so points recorded near a junction/wide-circuit region get merged
      more aggressively than narrow, well-defined ones) populated by `compute_dxa_circuit_sweep`
      (`append_point`/`commit_new_segment`/`do_merge` all updated to track it alongside `line`).
      Coarsening: adaptive merge-group sizing (`target_point_interval`, default 2.5, OVITO's own
      default) with open-segment endpoints always pinned (so junction connectivity isn't disturbed)
      and a proper closed-loop "seam" point. Smoothing: 2D Taubin (SIGGRAPH 95) alternating
      lambda/mu relaxation, `target_smoothing_level` (default 1, OVITO's own default) iterations.
      One subtlety found tracing OVITO's own C++ scoping precisely: the closed-loop closing point
      deliberately reuses the two boundary half-interval passes' own accumulator (not the middle
      loop's last group) -- OVITO's middle loop declares its own shadowing local variables of the
      same name, easy to miss porting from the raw source without noticing the shadowing.

      **Verified real quadrupole case**: 3018 raw points -> 141 coarsened+smoothed points (~21x
      reduction), total length 953.96 -> 599.69 Å (a real 37% reduction from tortuosity, not points
      lost -- raw sweep moves genuinely zig-zag, especially near junctions). **Screw dipole
      regression check** (an already near-straight line, minimal real tortuosity to remove): 803 ->
      34 points, but length barely changes (220.39 -> 220.01 Å, 0.17%) -- confirms the algorithm
      preserves a genuinely straight line's own length and only shortens real zig-zag, not a
      systematic bias. Segment/dislocation count is unaffected either way (coarsening runs strictly
      after the sweep, touches point shape only, not topology) -- the run-to-run dislocation-count
      variance seen between these two test runs (19 vs 23 on the quadrupole) is the same pre-existing
      OMP-triangle-emission-order non-determinism already documented elsewhere in this file, not
      caused by this feature.

   **Not yet touched, still on the OLD (PTM+angle-snap) path**: `compute_dxa_edge_vectors`,
   `compute_dxa_tet_classification`, `compute_interface_mesh`, `compute_atomistic_interface_mesh`,
   `compute_dxa_circuit_sweep`, and everything below in this file describing them -- none of that is
   broken or changed yet, this new work is purely additive so far (new files only, nothing rewired).
   Read this item first; the rest of the file describes the pre-existing (still currently used)
   pipeline until stages 2-4 above land.

1. **Which interface mesh should Burgers circuits actually be built on?** OVITO's own documentation
   states the interface mesh is "those triangular Delaunay facets having a good tetrahedral element
   on one side and a bad element on the other" — i.e. `compute_interface_mesh`'s tet-boundary
   construction, not the atom-centric `compute_atomistic_interface_mesh` (built from DXA1.3.6's
   source, which turns out to implement the *older*, 2010 Stukowski-Albe algorithm, not the current
   2012 Stukowski-Bulatov-Arsenlis one the documentation describes and OVITO actually runs).
   *However*, empirically the tet-boundary mesh gives 5164 triangles (even with the now-corrected,
   documented edge-resolution tet criterion) vs. the atomistic mesh's 2020 — OVITO's real reference
   is 1648, so the atomistic mesh is closer despite being architecturally the "wrong" one per the
   docs. Current decision (user-confirmed): build circuits on the atomistic mesh anyway
   (`compute_dxa_mesh_burgers_circuits`, see below) since it's tighter/less noisy, revisit if this
   turns out to matter once line-sweeping is working. Separately, the atomistic mesh's own
   remaining ~22% gap vs. OVITO (2020 vs 1648) is still open — see the "not yet tried" list in its
   own section below (DXA1.3.6's real neighbor-sorting algorithm, or its actual half-edge
   `removeUnnecessaryFacets`/`duplicateSharedMeshNodes`/`fixMeshEdges` machinery).
2. **DXA steps (vi)-(ix): done and verified against the actual 2012 paper, cross-segment merging
   added, one remaining gap.** `compute_dxa_circuit_sweep` implements the real algorithm end-to-end:
   local bounded-BFS trial-circuit search (step vi), Burgers vector from the seed circuit's own
   edges (step vii), halfedge-based sweep with shrink-before-expand priority and explicit per-facet
   ownership (step vii), self-closure/junction stop detection (step ix) — see its own section below
   for the full history (two prior bugs found and fixed: a first attempt missing the move-
   priority/facet-ownership mechanism entirely, built against the OVITO manual's summary rather than
   the paper; a second seeding from a global spanning-tree signal that doesn't localize defects,
   giving ~zero Burgers vectors for most lines) — plus Union-Find merging of raw segments across
   clean two-way (non-branching) junctions into `DXADislocationLines::dislocation_id` groups, since
   one physical dislocation is often discovered as several independently-seeded segments (see the
   "Follow-up: merging segments into physical dislocations" subsection below). Result on the real
   quadrupole: ~40 raw segments, all with physically real Burgers vectors (0.58-0.88, vs. BCC Ta's
   real a/2⟨111⟩=0.866), merged down by only a couple in this particular (genuinely branchy) network.
   **One gap remains**: a real multi-way junction's shared node position/connectivity isn't
   reconstructed — arms are correctly kept as separate segments/ids, but
   `DXADislocationLines::junction_vertices` stays empty (only counted/classified, not assembled into
   a shared coordinate). Not needed for length statistics, would matter for skeleton-graph rendering
   or node-degree analysis downstream.

   **Follow-up (2026-08-03): attempted a fix for the visible symptom of this gap — junction endpoints
   not coinciding in space — REVERTED, wrong approach identified.** User noticed (visually, via
   Paraview) that lines meeting at a junction don't share an exact endpoint; measured a real 0.5-2 Å
   gap per junction on the quadrupole by comparing every segment's raw endpoint against its nearest
   other-segment endpoint in the `.ca` output. First attempt: a post-growth pass that clustered every
   unmerged Junction-stopped node with whichever other node its own `blocking_nodes` resolved to
   (through `resolved_distinct_blockers`, i.e. through the merge-resolution chain), then snapped every
   member's endpoint to the cluster's average. This measured as 0.0 Å gaps everywhere on one run — but
   the user immediately caught (visually) that it was wrong: it mixed up which points belong to which
   junction, adding spurious excess line length. **Root cause of the wrong fix**: `resolve_node`/
   `resolved_distinct_blockers` exist to answer "which segment-chain does this territory currently
   belong to" for MERGE bookkeeping, not "where is this node physically right now" — `facet_owner[T]`
   freezes whichever node claimed T at whatever point during ITS OWN growth that happened to be, which
   can be long before that node's own eventual final stopping position (it may have kept growing well
   past T before halting somewhere else). Resolving a blocker through the merge chain can land on a
   node whose CURRENT final endpoint is the far end of an already-grown, already-merged chain,
   potentially spatially unrelated to the actual junction location — averaging that in corrupts the
   geometry instead of reconciling it.

   **Second attempt, same day: fixed correctly, checked directly against OVITO's real
   `DislocationTracer::joinSegments`** (`ovito/src/ovito/crystalanalysis/modifier/dxa/
   DislocationTracer.cpp`, lines ~1105-1327, read at the user's request before retrying). OVITO's real
   mechanism: builds a `junctionRing` (circular linked list) from `Edge::circuit` -- which node
   currently owns the *live boundary edge* right there, right now -- not from historical interior-face
   ownership; for `armCount>=3` it computes the ring's average position and **extends** every arm with
   one brand new point reaching that shared center (`line.push_back(...)`), it never overwrites the
   real last recorded point. Reimplemented with the equivalent distinction already present in this
   operator's own two ownership layers (see this file's own header comment): use `edge_owner` (live,
   continuously updated, matches `Edge::circuit`), never `facet_owner`/the merge-resolution chain
   (permanent historical record, matches nothing OVITO uses for this). For a stopped node's own final
   loop boundary edge (a,b), whichever node currently owns the exact reverse edge (b,a) is
   definitionally still sitting right there (or frozen exactly where it stopped, since `do_merge` never
   touches `edge_owner`/`facet_owner`) -- cluster only nodes that are BOTH still-standalone (unmerged)
   Junction survivors and directly, currently share such a live reverse edge, then **append** (not
   overwrite) one new point per arm at the cluster's average position, duplicating that arm's own last
   `core_size` for the new point (matching OVITO's `coreSize.push_back(coreSize.back())`).

   **Verified carefully this time** (a bit-exact before/after A/B wasn't possible -- the interface mesh
   itself varies run-to-run even at `OMP_NUM_THREADS=1`, upstream of this operator entirely -- so
   verified the fixed run's own output directly instead, across 5 separate runs): every junction
   endpoint gap is **exactly 0.0 Å every time** (segment counts varied 8-11 run to run, gap was 0.0 Å
   regardless); the newly appended point's own jump distance is always in the same 2-8 Å range as every
   other consecutive-point jump in that same line (checked point-by-point) -- i.e. a small, physically
   reasonable last hop, not the wild spatially-unrelated jump the first (reverted) attempt produced;
   total raw sweep length stayed in the same 800-820 Å ballpark across all 5 runs, no blow-up. Screw
   dipole regression unaffected either way (0 junction stops in that test case, confirmed unchanged:
   exact 2/2, ~110/110 Å).
3. **Line count still doesn't match the quadrupole's real ground truth — over-fragmentation,
   substantially reduced but not fixed.** User provided OVITO's actual `.ca` (Crystal Analysis)
   output for both test cases (`ovitodata/output_dxa_no_coarsening.ca`,
   `output_dxa_screw_dipole_no_coarsening.ca`), giving *exact* ground truth (parsed directly,
   including decoding the `DISLOCATION_JUNCTIONS` circular-linked-list format from OVITO's own
   exporter source): **quadrupole = 11 dislocations** (8 classical a/2⟨111⟩ at |b|=0.866 exactly,
   lengths 85-107 Å; 3 junction-type a⟨100⟩ at |b|=1.0 exactly, lengths 11-19 Å; connected via 7
   junction nodes — six 3-way, one 4-way, **zero clean 2-way pass-throughs anywhere**);
   **screw dipole = exactly 2 dislocations** (|b|=0.866 exact, ~114.5 Å each, self-closing via the
   periodic wrap).

   User also obtained OVITO's own actual (non-public) DXA source (`DislocationTracer.cpp`) directly,
   which let this be diagnosed against the real algorithm instead of just the paper. Confirmed two
   structural gaps versus the first sweep implementation: (1) territorial exclusion (already-claimed
   mesh edges/facets) must be baked into the *seed search* itself, not just the sweep; (2) segments
   must grow in lockstep via an outer loop over increasing trial-circuit length, not be swept to full
   completion one at a time. `compute_dxa_circuit_sweep.cpp` was rewritten around both. This alone,
   with only a post-hoc (end-of-run) merge step, made the screw dipole *perfect* (2/2) but made the
   quadrupole *worse* (116, up from the original 41) — diagnosed as: once two adjacent redundant
   seeds mutually block each other, post-hoc merging correctly records them as one dislocation for
   length statistics, but neither keeps *growing through the gap* afterward, inviting more redundant
   reseeding nearby each round. Fixed by moving 2-arm merging *inside* the growth loop (run every
   round, splicing the absorbed segment's far end onto the survivor as a still-growing node) —
   matching what OVITO's own `joinSegments` does. Result: quadrupole down to **69** (repeatable
   improvement: 116 → 87 with post-hoc merge only → 69 with incremental merge-and-continue). Screw
   dipole stays exactly **2/2**.

   **Still 69 vs. 11 on the quadrupole (~6x too many)** — since the real topology has *zero* 2-way
   pass-throughs, the remaining fragments are almost certainly still redundant, independently-seeded
   duplicates of the same physical arcs that happen to cluster in groups of 3+ near each other
   (correctly refused a merge by the 3+-arm rule, since the code can't tell a genuine 3-way junction
   from 3 spurious neighbors). Burgers-vector character: 62/69 cluster near ⟨111⟩ (0.5-0.9), only 3
   show any hint of the axis-aligned junction character (magnitudes 0.8-1.2) vs. the true 8-classical
   + 3-junction split — still not resolved. Two documented simplifications kept from the first
   version (simpler 2-move set instead of OVITO's 5; no `createSecondarySegment`) are the leading
   suspects for the residual redundant seeding, not yet tried.

   Checked in with the user at this point (69 vs 11, real progress but not closing the gap) rather
   than continuing to iterate further without confirming direction — see
   [[feedback_dxa_collab_style]] for why.

   **Follow-up, same day: two more real bugs found and fixed, and a diagnostic that changed the
   plan.** First, verified precisely (not just guessed) that the 69 fragments were the right kind of
   problem: mapped every one of them against OVITO's exact ground-truth line positions and found
   **100% precision** — all 69 landed within 0.3-2 Å of one of the 11 real lines, zero false
   positives anywhere else. So the seed search and territory are correct; only fragment count is
   wrong. User then asked for a `write_dxa_ca_file` operator to load our own output directly into
   OVITO (built — see its own section below) and pushed for continued exact-fidelity work rather
   than a pragmatic geometric-proximity consolidation shortcut, backed by another concrete
   observation: OVITO's own uncoarsened `.ca` output has ~95 points for one 85 Å dislocation (i.e.
   its raw sweep runs ~95 elementary moves *without ever stopping*), while `line_coarsening` is a
   separate, unrelated point-density-reduction step applied *after* tracing (same 11 dislocations
   with or without it) — ruling out "the missing piece is a coarsening pass" and confirming the fix
   has to be in the sweep mechanics themselves.

   Found two concrete bugs by direct comparison with the exact move preconditions in
   `DislocationTracer.cpp`:
   - **`tryRemoveTwoCircuitEdges` ("remove-2") was being explicitly skipped, not just omitted.** The
     shrink-move code had `if(a==c) continue;` for the exact vertex pattern (a "spike" a→b→a in the
     loop) this move exists to collapse — a plausible source of real fragmentation on the
     quadrupole's irregular mesh (vs. the screw dipole's clean uniform tube, which never develops
     spikes and worked perfectly even before this fix). Implemented it, but naively let a size-4
     loop shrink to a degenerate 2-vertex pseudo-loop, which broke the screw dipole (2→4) via
     downstream modular-index assumptions expecting >=3 vertices — fixed by requiring the
     precondition to hold for size >=5 in this vector-based (not raw pointer-spliced) implementation.
   - **`trySweepTwoFacets` was entirely missing** — the move that slides the loop boundary sideways
     across two adjacent unclaimed triangles sharing a common far apex, needed when neither triangle
     individually offers a valid single-facet move. Implemented directly in vertex-list terms (find
     the two triangles owning consecutive loop edges (a,b),(b,c); if their two "third vertices" are
     the same vertex w and the resulting new edges aren't already claimed, replace b with w in
     place — same loop length, both facets claimed). `tryRemoveThreeCircuitEdges` ("remove-3") was
     checked and found structurally impossible to trigger here (it requires the loop to revisit a
     vertex it's already visited a few steps earlier, which the expand move's own self-intersection
     guard already prevents) — not a gap, just dead code if implemented.

   **Result: quadrupole 69 → ~19-20** (run-to-run variance from the known non-deterministic mesh
   triangle order), screw dipole **stays exact at 2/2, now needing zero merging at all** (both seeds
   grow to their full ~114 Å length directly, matching ground truth almost exactly without ever
   fragmenting). Re-ran the same spatial-mapping check against ground truth: most of the 11 real
   lines are now recovered as 1 clean piece; a handful are still split into 2-3, and at least 2-3
   fragments show total length *exceeding* their own ground-truth line's full length — a sign of at
   least one incorrect cross-dislocation merge somewhere, not just remaining under-fragmentation.
   Not yet root-caused. `write_dxa_ca_file` (new operator, see below) writes our own result straight
   into OVITO's own file format for the user to inspect visually going forward, rather than only via
   ad hoc point-cloud comparison scripts.

   **Follow-up, same day: root-caused the length-overshoot bug — real, fixed it, then found a
   second, deeper one that's still open.** User asked to pursue both the overshoot root-cause and
   full OVITO move-set fidelity in parallel ("DO both").
   - **Bug found and fixed: `facet_owner` staleness after a merge.** `facet_owner[T]` freezes
     whichever node first claimed triangle T. Once that node gets absorbed into a two-arm merge, the
     frozen id is stale — a genuine third arm arriving later at a real 3-way junction would record a
     `blocking_node` pointing at an already-retired node, so it was never counted against the
     already-committed 2-arm merge's own exclusivity check, letting two of the three real arms get
     incorrectly spliced together before the third (often slower-growing, e.g. a short a⟨100⟩
     junction segment) revealed the true 3-way topology. Fixed with a `resolve_node` indirection
     (each merged-away node points to the surviving node representing that chain, chased/compressed
     like a Union-Find) applied everywhere `blocking_node` affects a merge decision. Also added a
     grace period (a candidate 2-arm merge must look exclusively mutual for several consecutive
     rounds before committing) as an extra safety margin against timing races between differently-
     paced arms.
   - **A second, deeper, still-open bug found while isolating the first one.** Directly tested
     with ALL merging disabled (not even the final forced pass): a single, entirely unmerged, raw
     grown segment still had 81 of its own points tracing one real ground-truth line and 10 tracing
     a *different* one it shares a real junction with — proving the chimera problem isn't (only) a
     merge-decision bug at all: a circuit's own growth can drift, step by step, across a real 3+-way
     junction from encircling one dislocation's core into encircling a neighboring one's, with every
     individual shrink/expand/sweep-two-facets move staying locally valid throughout (no facet
     double-claim, no size blowup — `n_stop_maxlen` is 0). Nothing in this implementation
     re-validates that a move keeps the circuit's own local elastic mapping self-consistent — this is
     exactly the kind of drift OVITO's Frank-rotation check exists to catch, which this operator
     omitted as "not needed for a single-grain system" (true for the *original* purpose of that
     check — grain-boundary compatibility — but it turns out to also serve as a general "is this
     circuit still encircling the same real defect" guard, which single-grain systems need too).
     Implementing an equivalent would mean tracking per-vertex local lattice orientation drift during
     growth, not just a one-time compatibility test — substantially more work than anything else in
     this file, and not yet attempted. See `compute_dxa_circuit_sweep.cpp`'s own "KNOWN OPEN BUG"
     header comment for the full technical detail.
   - **`createSecondarySegment` implemented, then measured and disabled.** User asked for it anyway
     ("the other half of 'do both'") even after the drift bug was found. Implemented faithfully:
     every round, each dangling node's own boundary is scanned for edges whose opposite side is
     genuinely unclaimed (a "hole"); that hole's own perimeter is walked, and if it's a valid,
     nonzero-Burgers loop bordering >=2 distinct known segments, it's committed as a new segment.
     Measured result on the quadrupole: **made things worse, not better** (16-20 → 37 physical
     dislocations), and 8 of the 37 were exactly degenerate (a single point, zero length) — spurious
     tiny "holes" that pass the nonzero-Burgers-vector + touches-2-segments test without being real
     defects. This is the same root gap as the drift bug above: OVITO's own version guards this with
     its Frank-rotation consistency check on the hole's own loop, which this operator doesn't have.
     **Disabled** (the loop body is intact and documented in `compute_dxa_circuit_sweep.cpp`, gated
     off with a `false` condition) rather than shipped as a regression — re-enable once an equivalent
     consistency check exists to reject spurious holes.
4. **Parameter/algorithm consistency check against OVITO's real DXA workflow** (user request,
   checked by reading `DislocationAnalysisModifier`/`StructureAnalysis`/`ElasticMapping` directly).
   **Already correctly matched**: our `compute_cna.cu` already implements the same *adaptive*
   per-atom cutoff CNA as OVITO's real `StructureAnalysis::determineLocalStructure` (same formula
   constants, e.g. `(1+√2)/2`) — `rcut` is just the neighbor-search radius the adaptive cutoff is
   computed within, not a naive fixed classification cutoff. `max_circuit_length`/
   `circuit_stretchability` defaults (14/9) exactly match OVITO's own hardcoded defaults.
   **Real differences found**: (a) OVITO computes lattice orientation via its own per-*cluster*
   least-squares fit tied directly to CNA-identified bonds (`StructureAnalysis::identifyStructures`);
   we use a separate `compute_ptm` per-*atom* RMSD fit layered on top of CNA classification — a
   plausible contributor to the Burgers-vector magnitude undershoot seen on the screw dipole. (b)
   **OVITO assigns each tessellation edge its ideal lattice vector via `CrystalPathFinder`** — a
   graph walk connecting two atoms entirely *through the good crystal region* (never stepping
   through a defective atom), robust even when the two atoms aren't direct neighbors. Our
   `compute_dxa_edge_vectors` instead requires *both* endpoints to individually be good crystal and
   does a direct angle-tolerance snap using one endpoint's own orientation. This is a concrete,
   promising lead for the open growth-drift bug above (item 3), since `CrystalPathFinder` is
   specifically designed to stay robust near defects/junctions — exactly the failure regime found.
   Not yet traced through to whether/how it would change `compute_atomistic_interface_mesh.cpp`'s
   own edge-vector derivation (what the circuit sweep actually consumes) — worth investigating next.
5. **FCC/HCP support for `compute_atomistic_interface_mesh`** — BCC-only right now. Needs the "8
   Thompson tetrahedra" template (different from BCC's 6 quads) — tables already extracted from
   DXA1.3.6 during this investigation (see git history / session log), just not transcribed into
   the operator yet.

**Superseded** by item 0's full elastic-mapping re-architecture (stages 1-4), which measurably beats
this on every metric (0.5% interface-mesh error vs. this path's 2.5-3x, exact classification match).
Recommended pipeline now (see `data/regression_new/dxa/compute_dxa_elastic_sweep_real_case.msp`):
`compute_dxa_lattice_correspondence` (discrete classifier, needs `target_structure`) →
`compute_dxa_lattice_clusters` → `compute_delaunay` → `compute_dxa_crystal_path_edge_vectors` →
`compute_dxa_elastic_mapping_tet_classification` → `compute_dxa_elastic_interface_mesh` →
`compute_dxa_circuit_sweep` → `smooth_dxa_dislocation_lines`. The paragraph below (PTM+CNA+angle-snap
path) is kept for historical context only -- its own regression `.msp` files were removed in the
2026-08-03 `data/regression_new/delaunay` cleanup (218MB -> 17MB) since nothing exercises this path
anymore; the operators themselves (`compute_dxa_edge_vectors`, `compute_dxa_tet_classification`,
`compute_interface_mesh`, `compute_atomistic_interface_mesh`) are still in the codebase, just
unexercised by any current test.

`compute_ptm` (orientation only) → `compute_cna`/`cna_fields` (classification) →
`ghost_update_opt` on both `cna_type`+`ptm_orientation` → `compute_delaunay` →
`compute_dxa_edge_vectors` with `struct_field: cna_type` → `compute_atomistic_interface_mesh`.

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

**Verified** (test file removed in an earlier cleanup, same 16000-atom BCC Ta lattice as
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

**Superseded**: the "~1.7x gap" framing above predates the `compute_cna` classifier finding
further down this file (PTM was simply under-classifying defects, not a mesh-construction issue).
See "Remaining work" at the top of this file for the current, accurate state.

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
identical file (`data/regression_new/dxa/ovitodata/output_cna_ovito.xyz`): **1026** non-BCC
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

**Current best pipeline** (test file removed in the 2026-08-03 cleanup, superseded by the elastic-mapping pipeline, item 0 above):
`compute_ptm` (still needed for `ptm_orientation` — `compute_cna` doesn't produce per-atom
orientation) → `compute_cna`/`cna_fields` (classification) → `ghost_update_opt` on both
`cna_type`+`ptm_orientation` → `compute_delaunay` → `compute_dxa_edge_vectors` with
`struct_field: cna_type` → `compute_atomistic_interface_mesh`.

**Attempted closing the remaining 22% gap — one fix was wrong, one was neutral, the gap itself
is elsewhere:**
1. *Tried*: quad-level dedup, to catch a hypothesized case of two different crystalline atoms both
   emitting a full quad over the identical 4 hole vertices (via different diagonal splits, so the
   existing exact-triangle dedup misses it). **Actively made things worse** (2020→2000 triangles
   but 34→47 open edges) and was reverted. Root issue: that situation is usually *not* a duplicate —
   it's two genuinely distinct facet sheets (e.g. the top and bottom surface of a 1-plane-thin
   disordered layer) that happen to coincide at the same 4 vertex positions. Deleting one side is
   exactly the wrong move; DXA1.3.6's own `duplicateSharedMeshNodes` handles this case by *splitting*
   the shared nodes so both sheets keep their own facets, not by picking one to delete.
2. *Tried*: generalized the hole-closer from "clean simple loops only" to a bounded DFS/
   backtracking search that also covers branch points (3+ open edges at a vertex), matching
   DXA1.3.6's `constructFacetRecursive` in spirit (minus its Burgers-vector-zero validity check,
   which needs an ideal vector for hole-to-hole edges that doesn't exist here). **Measured zero
   additional loops closed** on the real quadrupole case — found exactly the same 10 loops the
   simpler unique-outgoing-edge walker already found. Also tried raising the max loop length from
   8 to 20 edges: still zero additional closures, so the remaining 34 opens aren't loop-length-
   limited either. Kept anyway (strictly more general, no measured downside), but it isn't the
   lever here.
3. **Conclusion**: the 22% gap (2020 vs 1648) is *not* dominated by hole-closing/dedup gaps at the
   scale these two fixes could reach — 34 open edges out of ~3000 edge-uses is too small a
   contributor on its own. The likely remaining candidates are more fundamental: (a) DXA1.3.6's own
   `orderBCCAtomNeighbors` uses bond-connectivity + orientation-consistency checks to sort a
   crystalline atom's neighbors into canonical slots, not plain nearest-ideal-direction snapping
   (this operator's `resolve_slot`) — could disagree on slot assignment for atoms right at the
   disordered boundary, where directions are least clean; (b) the real `removeUnnecessaryFacets`
   (quad-diagonal flipping + genuine redundant-facet-pair removal, a proper half-edge operation, not
   the vertex-set-level dedup tried and reverted above) and `duplicateSharedMeshNodes`/`fixMeshEdges`
   machinery is still fully unported. Both are bigger investments than what's been tried so far.

See "Remaining work" at the top of this file for what's next (line extraction can proceed against
this mesh as-is; the 22% gap and FCC/HCP support are tracked there too).

Test file removed in an earlier cleanup (was pure-lattice CNA validation),
`compute_dxa_real_case_atomistic.msp` (atomistic mesh, PTM-driven classification — the 91-triangle
under-coverage case), `compute_dxa_real_case_cna_atomistic.msp` (atomistic mesh, CNA-driven — the
2020-triangle result), `ovitodata/output_cna_ovito.xyz` + `ovitodata/output_dxa_ovito.vtk`
(user-supplied OVITO ground truth: 1026 non-BCC atoms, 1648-triangle reference mesh).

## `compute_dxa_dislocation_lines` + `write_dxa_dislocation_lines`: DXA steps (viii)-(ix), first cut — SUPERSEDED, architecturally wrong

**Correction (user, citing `ovito/doc/manual/reference/pipelines/modifiers/dislocation_analysis.rst`,
"Technical background")**: this whole approach is built on a misunderstanding of how DXA actually
constructs Burgers circuits and dislocation lines. The documentation is explicit: trial circuits
are closed sequences of **interface mesh edges** (not raw atom-to-atom hops through the disordered
atom population), found by enumerating circuits in order of increasing length until one has a
non-zero Burgers vector, and the dislocation line itself is generated by *sweeping* that circuit
along the mesh, taking its **center of mass at each step** — none of which this operator does. It
instead thins the raw hole-atom adjacency graph directly, which is exactly why it fragmented into
293 tiny noise-dominated pieces (see below): raw atom positions/connectivity are noisy in a way a
well-defined 2D surface isn't. Kept in the codebase for now (not deleted) since the general
topological-thinning code may still be reusable, but **do not use this for real results** — see
`compute_dxa_mesh_burgers_circuits` below for the corrected direction, and "Remaining work" at the
top of this file for what's still missing (the actual sweep/line-tracing).

Extracts 1D dislocation line geometry and junctions from the disordered ("hole") atom population.
New files: `src/delaunay/compute_dxa_dislocation_lines.cpp` (extraction),
`src/delaunay/write_dxa_dislocation_lines.cpp` (VTK output), `src/delaunay/include/exaStamp/
delaunay/dxa_dislocation_lines.h` (`DXADislocationLines` struct).

**Design choice: atom-graph-based, not a surface-mesh sweep.** DXA1.3.6's own reference approach
(`DXATracing.cpp`, `burgersSearchWalkEdge`) sweeps an elementary Burgers circuit stepwise around the
interface mesh's own tube (a half-edge walk with per-step elastic-mapping re-evaluation). Skipped in
favor of working directly on the hole-atom adjacency graph instead: since the atomistic interface
mesh's own vertices already *are* the disordered atoms (`compute_atomistic_interface_mesh.cpp`),
their adjacency already encodes almost all the topology a mesh-sweep would have to rediscover, at a
fraction of the implementation cost. Traded off against a less rigorous notion of "skeleton" (graph
thinning approximates a medial axis for tube-like topology, but isn't a rigorous one).

**Algorithm**: build the hole-hole adjacency graph (from `DXAEdgeVectors`, reusing its
`vertex_matches_target`/already-deduplicated edges) → **topological thinning**: repeatedly find any
atom with graph-degree > 2 whose removal would not disconnect its neighbors (i.e. not an
articulation point — Tarjan's algorithm, recomputed fresh after every single removal for
correctness, cheap at this data scale) and remove it, until none remain → trace the thinned graph
(walk from every degree≠2 vertex along its own edges to the next such vertex = an open line; any
edges left over form a pure degree-2 cycle = a closed loop with no junctions) → per-line Burgers
vector by averaging every `DXABurgersCircuits` "confirmed signal edge" whose midpoint's nearest
skeleton point belongs to that line.

**Why topological thinning, not fixed-depth erosion (the first version tried and discarded)**: the
hole-atom population is not a 1D chain — it's a genuine tube *surface* (atoms wrap around the
tube's own circumference, not just along its length). Measured both failure modes of a simpler
fixed-depth erosion (BFS distance from the crystalline-touching "interface" shell, keep atoms with
depth ≥ some threshold) on the real quadrupole case:
- `min_core_depth=1` (drop just the interface shell): eroded almost everything away — only **2 tiny
  4-point fragments** survived (length ~10 Å each) — the tube's local radius is barely more than 1
  atomic layer almost everywhere, so removing that one layer leaves nothing to trace a line through.
- `min_core_depth=0` (no erosion, use every hole atom): **1015 junctions out of 1026 atoms** — since
  every atom on a tube's own circumference has several neighbors *around* the ring in addition to
  along its length, almost everything looks like a "junction" (degree ≥ 3) even though it's just an
  ordinary point on a smooth tube surface.

Topological thinning fixes both failure modes at once, since it naturally adapts to locally-varying
tube radius instead of applying one uniform depth everywhere.

**Result on the real quadrupole case, and the open problem**: 1026 hole atoms → 419 skeleton atoms
→ **293 lines, 160 junctions**. The algorithm runs correctly and produces real structure (the
longest extracted lines, 27.5/22.2/22.0 Å, are plausible dislocation-segment lengths), but the
output is clearly **over-fragmented**: median line length is only 3.5 Å, and 169/293 lines are under
5 Å — way more pieces than the handful of lines a "quadrupole" (nominally ~4 dislocations) implies.
Likely cause: real atomic data has thermal/positional noise, and a strict graph-connectivity test
(articulation points, computed from an exact distance-cutoff-based edge graph) is sensitive to small
local irregularities that create spurious "bottlenecks" — genuinely disconnecting in the *discrete
graph* sense, even though they're not real physical branch points. **Not yet implemented**: a
post-processing simplification pass (merge junction nodes within a small radius of each other into
one; prune/absorb short dangling segments below a length threshold into their neighboring line) —
needed before this output is actually usable, tracked in "Remaining work" at the top of this file.

Also confirmed via a quick clustering check (not in the codebase, one-off analysis) that the real
disordered population is **one single connected network**, not 4 separate lines — a genuine
"quadrupole" defect network with real junctions where segments meet, not 4 independent objects. This
validates using a junction-aware graph approach over a simpler per-component (e.g. PCA-per-blob)
method that was considered and rejected before implementing thinning.

Test file removed in the 2026-08-03 cleanup (superseded for line extraction, see item 2 above).

## `compute_dxa_mesh_burgers_circuits`: DXA steps (vi)-(vii), corrected — circuits on the interface mesh itself

Stage 1 of the correction above (working in stages, checking in after each, per user's direction).
Implements the actual documented mechanism: a Burgers circuit is a closed sequence of **interface
mesh edges**, and its Burgers vector is the sum of their ideal lattice vectors. `compute_
dxa_burgers_circuits.cpp`'s existing spanning-tree/fundamental-cycle mechanism is exactly the right
*technique* for "enumerate circuits in order of increasing length" (a spanning tree's non-tree
edges each close exactly one fundamental cycle — the shortest one through that specific edge — with
no combinatorial search needed) — it was just applied to the wrong graph (the whole crystal's
resolved-edge graph, not the interface mesh's own edges). New operator `compute_dxa_mesh_burgers_
circuits.cpp` does the same computation restricted to the mesh.

**Prerequisite found along the way**: interface mesh edges connect two *non-crystalline* atoms,
which have no orientation/ideal-vector of their own — so where does an edge's ideal vector even
come from? DXA1.3.6's own source answers this: from the **generating crystalline atom's own
resolved template slots** (`latticeVectors[v1]-latticeVectors[v]` in its `InterfaceMesh.cpp`), not
from either mesh vertex. `compute_atomistic_interface_mesh.cpp` already computes exactly this
per-atom slot resolution when building triangles — it just wasn't exposing it. Extended
`InterfaceMesh` with a new `edge_ideal_vector` field (v0<v1 canonical direction, first-seen-wins if
two different generating atoms disagree slightly) and `compute_atomistic_interface_mesh.cpp` now
populates it directly from the same slot data already computed for triangle emission. Verified this
addition doesn't change the triangle count (still 2020 on the real case) before moving on.

**Result on the real quadrupole case**: **666 confirmed signal edges out of 2986** mesh edges that
have a resolved ideal vector (22.3%). Sanity-checked the actual Burgers vector magnitudes against
known crystallography: BCC Ta's real full dislocation is a/2⟨111⟩, magnitude 0.866 in
`bcc_ideal_raw`'s own units — the measured histogram peaks in the 0.7-0.8 bin (202/666 edges), with
a plausible tail out to ~1.7-1.8 (likely junction regions where a circuit inadvertently encloses
more than one dislocation's worth of signal). This is a real, physically-plausible signal — a sharp
contrast with the previous (wrong) approach's output.

**Answering "why not build circuits on the atomistic (thin) mesh, since it's already closer to
OVITO's triangle count than the tet-boundary one?"**: this operator does exactly that — it consumes
`compute_atomistic_interface_mesh`'s output specifically (not `compute_interface_mesh`'s
tet-boundary one), for exactly the reason above (the atomistic mesh's vertices are the actual
disordered atoms, giving a much tighter, less noisy surface to build circuits on than the coarser
tet-boundary mesh, which independently confirmed 5164 triangles with the exact-documented
edge-resolution tet criterion — see "Remaining work" item 1 at the top of this file, still open).

**Not yet done (next stage)**: the actual circuit **sweep** — advancing a confirmed signal circuit
step-by-step along the interface mesh (an "advancing front" on the triangulated surface) and
recording its center of mass at each step as the dislocation line's vertex, plus junction handling
where sweeps merge or the circuit needs to stretch past a kink. This is what actually produces
smooth 1D lines instead of just a set of flagged edges. `compute_dxa_dislocation_lines.cpp` (the
superseded operator above) should eventually be replaced by this, not extended.

Test file removed in the 2026-08-03 cleanup (operator still in the codebase, just unexercised).

## `compute_dxa_mesh_dislocation_lines`: DXA steps (viii)-(ix), stage 2 — seeded from confirmed circuits, still fragmented

Stage 2 of the correction (see `compute_dxa_mesh_burgers_circuits` section above for stage 1).
New operator: takes the confirmed signal edges from `compute_dxa_mesh_burgers_circuits`, builds a
"core" vertex set from their endpoints, pulls in the interface mesh's own full edge connectivity
among just those vertices (a signal edge is just one arbitrary closing edge of its own fundamental
cycle, not necessarily touching its geometric neighbors along the tube — surrounding mesh edges
restore that), then reuses `compute_dxa_dislocation_lines`' own validated topological-thinning +
tracing machinery on this much smaller, pre-validated graph instead of the raw noisy hole-atom
population. Output slot deliberately named `dxa_dislocation_lines` (same as the superseded
operator's own) so `write_dxa_dislocation_lines` auto-wires — don't run both operators together.

**This is explicitly an approximation of OVITO's documented advancing-front sweep** (see the
superseded section's own correction), not a literal port of it — ponytail-noted in the file's own
header comment. A real sweep would recompute and advance an actual circuit step by step, taking its
center of mass; this instead thins the validated core region once and uses skeleton atoms' own
positions directly. Revisit if results still don't look right after further tuning.

**Result on the real quadrupole case**: 552 core atoms (from 590 signal edges) → 329 skeleton atoms
→ **177 lines, 107 junctions**. A real, measured improvement over the superseded raw-hole-atom
approach (293 lines / 160 junctions): median line length **3.5→5.45 Å**, longest line **27.5→40.2
Å**. Still fragmented (82/177 lines under 5 Å) — confirms seeding from physically-validated signal
is meaningfully better than raw noisy connectivity, but doesn't fully solve fragmentation on its
own. Next candidates, not yet tried: a real advancing-front sweep implementation, or a
simplification pass (merge close junctions, prune/absorb short segments below a length threshold).

Test file removed in the 2026-08-03 cleanup (superseded for line extraction, see item 2 above).

## `compute_dxa_circuit_sweep`: the real advancing-front sweep — verified against the actual 2012 paper, working

Per user's explicit request to implement the actual sweep. First attempt (see git history) was
built against the OVITO *manual's* summary only, which turns out to omit the sweep's actual
mechanism entirely — it only says a hard length limit stops the circuit at a junction, with no
detail on *how*. That gap was diagnosed empirically (first attempt produced lines 700-1460 Å long
in a ~131 Å box, wandering the whole connected defect network) before checking the primary source:
Stukowski, Bulatov, Arsenlis 2012 (`/home/lafourcadep/Bureau/DXA_ref.pdf`, secs 2.4-2.6), which
turned out to specify a materially different and much more precise mechanism. Rewrote against that.

**The actual mechanism (verified via `pdftotext`-extracted, directly-quoted paper text)**:
- The interface mesh is a proper **halfedge structure**: each triangle (v0,v1,v2) owns 3 directed
  halfedges (v0→v1),(v1→v2),(v2→v0) — a manifold mesh has exactly one triangle owning any given
  directed pair.
- Seed circuit: genuine shortest cycle through a confirmed signal edge (BFS excluding the direct
  edge, bounded by `max_circuit_length`, OVITO's own default 14) — same as the first attempt.
- **Sweep, one elementary move at a time, with an explicit priority the first attempt didn't have**:
  "*Moves that reduce the length of the circuit are given precedence over moves that extend it.*"
  A **shrink** move: if 2 consecutive circuit halfedges (a→b),(b→c) are owned by the *same*
  unclaimed triangle {a,b,c}, replace them with the single edge (a→c) — cutting the corner off that
  triangle, claiming it. Only if *no* shrink is available does the circuit **expand**: absorb the
  third vertex of the (unclaimed) triangle owning one of its own halfedges.
- **Explicit per-facet ownership, which the first attempt didn't have at all**: "*Once a mesh facet
  has been traversed by an advancing circuit, it is marked as belonging to the current dislocation
  segment and no other circuit is allowed to sweep the same triangle again.*" A sweep genuinely
  halts — not via an arbitrary cap — when every candidate move's facet is already claimed (by
  itself: a closed loop; by another segment: a real junction) or doesn't exist (an open mesh edge).
  This is the piece that was actually missing, and it's why the first attempt's ad-hoc safety caps
  (path-length cap, "old ground" revisit detection) were symptom patches, not the fix — removed now
  that the real mechanism replaces them.
- Line vertex at each move = the circuit's own center of mass, literally as documented.

**Result on the real quadrupole case — dramatic, physically sensible improvement**: 43 lines
extracted, with **line lengths 1.7-82.2 Å (median 19.4 Å)** in the ~131 Å box — no more absurd
wandering. Sweep stop reasons are now real, diagnostic signals instead of an arbitrary cap: 53
max-length, **17 genuine self-closures** (the loop met its own earlier territory — a closed
dislocation loop or a very short segment), **16 genuine junction collisions** (met another
segment's claimed territory — an actual dislocation junction), 0 open-mesh-edge, 0
exhausted-no-move. `n_seeds_skipped` covers signal edges already consumed by an earlier line's
sweep.

**Still scoped out (ponytail, see file header comment)**: junction *connectivity* — splicing which
lines meet at a node into a proper multi-arm representation (matching the CA file format's own
circular-linked-list convention) — isn't built. The operator counts/classifies *why* each sweep
stopped, but doesn't yet record *which other segment* a junction collision was with, or assemble
that into `DXADislocationLines::junction_vertices` (left empty). This is the natural next piece if
junction connectivity is needed downstream.

`DXADislocationLines` gained a `line_positions` field (real `Vec3d`s, since this operator's line
vertices are synthetic swept-circuit centroids, not atom indices) — `write_dxa_dislocation_lines.cpp`
was extended to use it when present.

### Follow-up correction: seed discovery was borrowing the wrong signal, Burgers vectors were ~zero

The version above still seeded each sweep from `compute_dxa_mesh_burgers_circuits`' own confirmed
signal edges (a global, arbitrary-root spanning-tree residual) and used that residual directly as
the line's Burgers vector. Per user request, fixed to compute the Burgers vector properly: sum
`InterfaceMesh::edge_ideal_vector` directly around the operator's *own* seed loop — "the discrete
line integral over the initial forward circuit", literally as the paper specifies (§2.6), not a
value borrowed from a different computation.

**This surfaced a much bigger problem than a display value**: doing this revealed that **most seed
loops had ~zero Burgers vector** (37 of 41 lines) — only 4 showed a real signal. Root cause: a
global spanning-tree residual for edge (a,b) only proves *some* defect exists somewhere along the
(possibly very long) loop `root→tree-path→a→b→tree-path-reversed→root` — it says nothing about
whether the defect is anywhere *near* (a,b) itself. The genuinely local shortest cycle through (a,b)
usually doesn't enclose anything at all, because the real defect the spanning tree detected is
elsewhere along that long path. **Confirmed signal edges were the wrong tool for seeding a local
sweep** — a real, methodological bug, not a proxy-value inconvenience.

Fixed by implementing the paper's actual step (vi) directly, dropping the dependency on
`compute_dxa_mesh_burgers_circuits` entirely (this operator is now fully self-contained, needing
only the interface mesh): for each candidate mesh vertex, run one bounded-depth BFS (depth ≤
`max_circuit_length/2`) rooted at it; every non-tree edge found closes a genuinely *local*
fundamental cycle (bounded by construction, unlike the old global spanning tree); collect all such
candidates and keep the shortest one whose own Burgers vector exceeds `min_burgers_norm` —
equivalent to "circuits of increasing length until a non-zero one is found," computed in one BFS
pass instead of literally re-searching at each length. Iterate over every mesh vertex (skipping
ones already consumed by an earlier line) until the whole mesh is covered, matching step (viii).

**Result on the real quadrupole case — every single line is now real**: 39 lines, from 39 local
trial-circuit searches attempted (100% hit rate — expected here, since this atomistic mesh's
vertices *are* the disordered atoms themselves, so almost any starting point is genuinely near the
core network). **All 39 Burgers vector magnitudes cluster tightly between 0.58 and 0.88** — right
around BCC Ta's real a/2⟨111⟩ full-dislocation magnitude (0.866) — with individual components
consistent with the expected ⟨111⟩-type pattern. Line lengths mostly 0.7-101 Å (median 15.5 Å), but
~12 lines are suspiciously short (under 2 Å) — plausibly the "seed collides with existing territory
almost immediately" artifact flagged separately (not yet fixed, tracked as item 1 in "Remaining
work").

Test file removed in the 2026-08-03 cleanup (ran this same operator on the old, worse atomistic mesh -- superseded by `data/regression_new/dxa/compute_dxa_elastic_sweep_real_case.msp`).

### Follow-up: single screw dislocation dipole, a controlled test with no real junctions

To isolate whether the quadrupole's over-fragmentation (see "Remaining work" item 3) is a general
sweep/seeding bug or specific to that dataset's real junction network, generated a synthetic BCC Ta
sample with exactly two straight a/2⟨111⟩ screw dislocations (opposite sign, a periodic dipole — the
minimal way to embed a real dislocation under full 3D periodic boundary conditions). Built with `data/regression_new/dxa/gen_screw_dipole.py`: orthogonal simulation frame
x=[1,-1,0], y=[1,1,-2], z=[1,1,1] (line direction), lattice tiled exactly via rotate-and-crop of the
standard 2-atom BCC basis, then displaced with the exact isotropic elastic screw solution
`u_z = b/(2π) * (atan2(y-y0,x-x1) - atan2(y-y0,x-x2))` for two cores at the same y, separated along x
by half the box (a periodic-compatible dipole placement) — no relaxation run afterward
(`max_iteration: 0`), so atoms sit on the raw continuum-displaced positions.
File: `data/regression_new/dxa/screw_dislo_dipole.xyz` (128963 atoms, 140x146x114 Å box). Test
file removed in the 2026-08-03 cleanup (superseded by `data/regression_new/dxa/compute_dxa_elastic_sweep_screw_dipole.msp`).

**Result: the sweep/merge code found essentially the right answer.** 3 raw segments (not 43-like
fragmentation), 0 max-length stops: 2 real, long dislocations (146 Å and 184 Å) each correctly
extracted as a **self-closed loop** — exactly the expected topology for a straight line under PBC,
since sweeping along it eventually re-enters its own already-claimed territory after 1+ periodic
images. Their Burgers vectors point opposite directions (dominant component +0.68 vs -0.63),
consistent with a dipole, though undershooting the ideal a/2⟨111⟩=0.866 magnitude somewhat (no
relaxation was run, so the raw elastic core distorts the local circuit fit). The 3rd segment (1.15 Å,
open) is very likely a construction artifact, not a code bug: the dipole separation (70 Å) is
comparable to the distance from each core to the box edge (35 Å) rather than in the well-separated
far-field regime, so the naive two-term (non periodic-image-summed) displacement field doesn't
cancel exactly at the periodic boundary — visible as 234 open mesh edges (vs. 34 for the quadrupole)
and 749 failed local seed attempts (only 3 succeeded), i.e. a broad spurious non-crystalline ribbon
along the boundary that mostly (correctly) has zero Burgers vector.

**Conclusion**: on a case with no real junctions, the sweep/merge mechanism recovers almost exactly
the right topology. This points the quadrupole's 41-vs-12 over-fragmentation toward the
junction-handling / dense-network regime specifically (redundant reseeding along real branch
networks, as suspected), not a general defect in the sweep algorithm itself. If revisiting this test,
use a proper periodic-image-summed dipole field (or a much larger box relative to separation) to
remove the small 1.15 Å boundary artifact.

### Follow-up: merging segments into physical dislocations for length statistics

Per user request ("topologically a line can be an ensemble of segments... important to get the
final dislocation length statistics"): a single physical dislocation is very often discovered as
several separate raw segments, purely because they were seeded independently and happened to sweep
into each other's already-claimed territory (`StopReason::Junction`). That collision alone doesn't
mean a real 3+-arm branch — it only means "not the first segment to reach this facet". Reporting
each raw segment's length separately would badly fragment the true per-dislocation length.

Fixed with Union-Find: `sweep_move`/`sweep_direction` now also report `blocking_segment` — the id of
whichever *other* segment's claimed facet actually stopped a sweep (when `reason==Junction`). For
each segment B, `incoming_count[B]` counts how many other segments' sweeps were stopped by B. If a
segment A's sweep stopped at B and `incoming_count[B]==1`, there's no branching decision to make — A
and B are the same continuous line, `union(A,B)`. This chains transitively across longer clean
pass-through runs. A real multi-way junction (`incoming_count[B]>=2`) is left un-merged there on
purpose: each arm keeps its own id, since which of the >=2 incoming segments is "the real
continuation" isn't decidable from this signal alone. New `DXADislocationLines::dislocation_id`
field (one per raw segment/line) carries the merge grouping; `write_dxa_dislocation_lines.cpp` now
emits it as an `Int32` `dislocation_id` CellData array so ParaView can color/group merged segments.
New `n_dislocations` OUTPUT slot reports the post-merge count.

**Result on the real quadrupole case**: raw segment count varies run-to-run (41-43, since
`compute_atomistic_interface_mesh`'s OMP-parallel triangle emission order — and therefore this
operator's vertex-iteration seed order — isn't deterministic across runs); a representative run gave
43 raw segments, only 2 of which met the clean-pass-through merge criterion, giving 41 physical
dislocations, lengths 0.74-120.6 Å (median 9.6 Å, total network length 924.6 Å). Most junction stops
in this particular network land on a segment with `incoming_count >= 2` (a real multi-way hub, left
un-merged on purpose) rather than a clean 1-in pass-through — expected for a genuinely branchy
quadrupole dislocation network. So the merging logic is verified working and structurally correct,
but on this test case it only resolves a small minority of the raw segments; the shortest merged
dislocation is still ~0.74 Å, about the same order as before merging. The previously-flagged "~12
suspiciously short (<2 Å) lines" (item 1 in "Remaining work") are therefore only partly explained by
pass-through fragmentation — most are genuinely short arms terminating at a real multi-way junction,
not an artifact merging should eliminate.

## `write_dxa_ca_file`: writes our own result into OVITO's own `.ca` format

New operator, `src/delaunay/write_dxa_ca_file.cpp`. Writes a `DXADislocationLines` (only
`compute_dxa_circuit_sweep` populates the fields this needs: `line_positions`, `burgers_vector`)
into OVITO's own "Crystal Analysis" file format, reverse-engineered directly from OVITO's real
exporter/importer source (`CAExporter.cpp`/`CAImporter.cpp`, obtained by the user, not from public
docs) — so our result can be loaded straight into OVITO for visual side-by-side comparison against
its own DXA output, rather than only via ad hoc point-cloud comparison scripts.

Deliberately simplified relative to a real OVITO-written file:
- Only one bare `STRUCTURE_TYPE` stub is declared (enough for the importer's header parsing, not a
  faithful reproduction of OVITO's real per-structure Burgers vector family tables).
- `DISLOCATION_JUNCTIONS`: this operator doesn't reconstruct real multi-way junction connectivity
  (see `compute_dxa_circuit_sweep`'s own "ponytail" note), so every dislocation is written as a
  trivial self-referential 2-cycle — OVITO will render every line as an independent, unconnected
  segment, exactly what this operator actually knows, no more.
- No native multi-piece convention exists for `.ca` (unlike `write_dxa_dislocation_lines.cpp`'s VTK
  `.pvtu` pieces) — **fixed 2026-08-03** (was single-rank-only before, a real bug the user hit
  running with MPI>1: every rank raced to open/truncate the same output path independently). Now
  gathers every rank's own lines to rank 0 via `MPI_Gatherv` (serialized into a flat, self-delimiting
  `double` buffer: burgers.x/y/z, npoints, then npoints*4 doubles), only rank 0 writes the file,
  re-numbering every line 0..N-1 as it writes (the original per-rank index doesn't need to survive
  the gather). Verified with 1/2/4 MPI ranks on the real quadrupole case: loads cleanly in real OVITO
  Pro (`ovitos`) every time, segment count growing with rank count (15/18 vs single-rank's 11) exactly
  as expected — MPI domain decomposition splits the interface mesh, so a line crossing a rank boundary
  becomes 2 independent segments (one per rank's own local circuit sweep); this is the same, already-
  documented cross-rank fragmentation any per-rank-local DXA analysis has, not a gather bug.

Wired into `compute_dxa_elastic_sweep_real_case.msp`/`compute_dxa_elastic_sweep_screw_dipole.msp`,
writing to `ovitodata/our_quadrupole_result.ca` / `ovitodata/our_screw_dipole_result.ca` (same
directory as the user's own ground-truth `.ca` files, for easy side-by-side loading in OVITO).

**Follow-up, same day: "every dislocation shows up as 'Other' in OVITO" — real bug, found and
fixed, but not where suspected.** User noticed OVITO's own DXA classification always labeled our
dislocations "Other" and suspected `BURGERS_VECTOR_FAMILY`. That table WAS genuinely incomplete
(only the catch-all "Other" family was declared) and got fixed first — now a faithful copy of
OVITO's real 4-family BCC table (`Other`, `1/2<111>`, `<100>`, `<110>`, verbatim reference
vectors/colors cross-checked against a real OVITO-exported `.ca` file). **But this alone didn't fix
it** — measured directly via `ovitos` (`data.tables['disloc-lengths']`/`data.tables['disloc-counts']`,
the actual per-type length/count breakdown OVITO computes on import) that everything still classified
as "Other" after the family-table fix.

Root-caused by bisection against the real OVITO reference file rather than guessing: took the known-
good reference `.ca`, and one at a time (a) changed a dislocation's own Burgers vector to a
sign/permutation-symmetric variant, (b) forced `CLUSTER_ORIENTATION` to identity, (c) renumbered
`STRUCTURE_TYPE`'s own id, (d) spliced our own real dislocation records into the good header — none
of these broke classification. Only forcing `CLUSTER_SIZE` from its real value down to **0** broke
it completely (every dislocation instantly became "Other"). This operator had always written a
hardcoded `CLUSTER_SIZE 0` (it never had the real atom count available at all) — OVITO's importer
apparently treats a zero-size cluster as invalid for Burgers-vector-family matching, regardless of
how correct the family table itself is. Fixed by taking `DXALatticeClusters` as a new input: cluster
1's real `atom_count` (summed across MPI ranks via `MPI_Reduce`, since `compute_dxa_lattice_clusters`
is per-rank local) now goes into `CLUSTER_SIZE`, and its real least-squares `orientation` fit (rank
0's own local value — display/diagnostic only in OVITO, no cross-rank reconciliation needed) now
goes into `CLUSTER_ORIENTATION` instead of a hardcoded identity matrix.

**Verified via `ovitos`'s own classification tables, not just "loads without error"**: quadrupole
case now splits cleanly into `1/2<111>` and `<100>` with zero "Other", at 1 and 2 MPI ranks alike
(counts differ slightly per run from the same pre-existing cross-rank/mesh-order non-determinism
documented elsewhere in this file, but the "Other" bucket is empty every time). While in this code,
also fixed a smaller, unrelated gap noticed along the way: the per-point trailing value in the
`DISLOCATIONS` section (OVITO's own "core size") was hardcoded to a constant 0 for every point even
though `DXADislocationLines::core_size` already tracks the real value (same data
`smooth_dxa_dislocation_lines` consumes) — now written through properly.

**Follow-up, same day: header now declares all 5 of OVITO's own structure types, not just the one
used.** User pointed out that a real OVITO-written `.ca` always declares all 5 built-in structure
types (fcc/hcp/bcc/diamond/hex_diamond) even when a given DXA run only ever analyzes one of them.
This operator previously only declared the single "bcc" type it actually uses (as id 1, an arbitrary
choice). Now writes all 5 verbatim (same names/reference vectors/colors as OVITO's own file), with
`bcc` at OVITO's own real id (3) — `CLUSTER_STRUCTURE` updated to match. Purely a format-fidelity
change (already confirmed the numeric STRUCTURE_TYPE id doesn't affect classification correctness,
see the CLUSTER_SIZE finding above), verified via `ovitos`: `data.dislocations.crystal_structures`
now lists all 5 names correctly (plus OVITO's own built-in "Unidentified structure" at index 0), and
classification is still exactly correct (0 "Other").

**Follow-up, same session: built.** New operator, `src/delaunay/dxa_lattice_cluster_fields.cpp`,
copying `DXALatticeClusters::atom_cluster` (0 = unresolved) into a named per-particle grid field —
identical `exanb::compute_cell_particles` pointwise-copy pattern as `cna_fields`/`ptm_fields`
(`src/cna/compute_cna.cu`), just reading a `std::vector<int32_t>` member off `DXALatticeClusters`
instead of a standalone `CudaMMVector<double>` slot. Verified on the real quadrupole case via
`write_xyz` (`dxa_lattice_cluster_fields: { cluster_field: dxa_cluster_id }` right after
`compute_dxa_lattice_clusters`, then `write_xyz: { fields: [id, type, dxa_cluster_id] }`): the field's
own value distribution matched the operator's own printed stats exactly (126974 atoms at cluster id
1, 1026 at 0/unresolved, out of 128000 owned particles). Wired into `compute_dxa_elastic_sweep_
real_case.msp` right after `compute_dxa_lattice_clusters` (cheap, pointwise, no extra output file by
default -- only materializes the field, doesn't write anything on its own).

## `write_ovito_interface_mesh`: writes our own interface mesh in OVITO's own legacy VTK format

New operator, `src/delaunay/write_ovito_interface_mesh.cpp` (2026-08-03, user request). Unlike
`write_interface_mesh.cpp`'s own XML `.pvtu`/per-rank-piece convention, this instead byte-for-byte
matches the *legacy* VTK ASCII format OVITO Pro itself writes for its own interface-mesh export
(cross-checked directly against a real OVITO-exported reference file): `# vtk DataFile Version 3.0`
header, `DATASET UNSTRUCTURED_GRID`, `POINTS n double`, `CELLS n 4n` (`"3 v0 v1 v2"` rows,
`VTK_TRIANGLE`), `CELL_TYPES n` (all `5`), then both `CELL_DATA`/`POINT_DATA` carrying one
`SCALARS cap unsigned_char` field. Legacy VTK has no multi-piece convention at all (it's inherently
one self-contained file), so — unlike `write_interface_mesh.cpp` — this operator does a real
`MPI_Gatherv` of every rank's own (deduplicated, referenced-only) vertices and triangle connectivity
to rank 0, rebasing each rank's own local vertex indices by that rank's running point-count offset
before writing; only rank 0 touches the filesystem.

**Honest caveat on `cap`**: a real OVITO file can have `cap==1` triangles (confirmed on the real
screw-dipole reference: 506/3004 capped) — synthetic triangles OVITO adds to close the surface at a
non-periodic domain boundary. This operator's own `InterfaceMesh` never synthesizes such triangles
(a boundary edge is simply left open, see `InterfaceMesh::edge_triangles`'s own doc comment), so `cap`
is always written as `0` here. This matches OVITO's file *format* exactly (loads identically into
OVITO/ParaView), not its capping *semantics* — synthesizing real cap triangles would be a separate,
substantially bigger feature, not attempted.

**Verified** with the real quadrupole case, 1/2/4 MPI ranks, loaded via real `ovitos`
(`import_file(path, input_format="vtk/legacy/mesh")`, `data.triangle_meshes['mesh']`): single-rank
gives `vertex_count=816` (exact match to OVITO's own reference), `face_count=1660` (OVITO's own:
1648, the same small residual gap already documented elsewhere in this file); 2 ranks gives
824/1631, 4 ranks gives 837/1624 — every case loads cleanly with vertex/face counts matching the
file's own declared header exactly, confirming the cross-rank index rebasing is correct.
