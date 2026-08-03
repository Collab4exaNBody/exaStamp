/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <exaStamp/delaunay/lattice_structure.h>
#include <exaStamp/delaunay/cluster_graph.h>
#include <vector>
#include <cstdint>

namespace exaStamp
{
  using namespace exanb;

  // Output of compute_dxa_lattice_correspondence (StructureAnalysis::determineLocalStructure port,
  // see that operator's file header): per atom, a DISCRETE neighbor-to-canonical-slot mapping,
  // found by matching the atom's own physical bond graph against the fixed reference tables
  // (lattice_structure.h) -- this replaces PTM's role in the DXA elastic mapping. Flat arrays,
  // indexed the same way as compute_ptm/compute_cna's own buffers
  // (grid->cell_particle_offset_data()); only OWNED particles are touched by the classifier itself,
  // but the array is sized to cover ghosts too so compute_dxa_lattice_clusters (which needs
  // neighbor identity across rank/domain-decomposition boundaries) can read a ghost's own
  // classification once its field is synchronized (ghost_update_opt), same convention as ptm_type.
  struct DXALatticeCorrespondence
  {
    std::vector<uint8_t> structure_type; // per atom: LatticeStructureType (0 = unmatched)

    // atom*DXA_MAX_NEIGHBORS + slot -> flat particle index of the physical neighbor occupying
    // canonical slot `slot` of this atom's OWN bond-topology match (see lattice_structure.h) --
    // NOT yet resolved to a globally-consistent absolute frame across atoms (that discrete
    // labeling ambiguity, up to the structure's own point-group symmetry, is what
    // compute_dxa_lattice_clusters' BFS cluster growth resolves). -1 where structure_type==0 or
    // the slot exceeds this atom's own num_neighbors.
    std::vector<int64_t> neighbor_atom;

    inline int64_t neighbor( size_t atom, int slot ) const { return neighbor_atom[ atom * DXA_MAX_NEIGHBORS + slot ]; }
  };

  // Output of compute_dxa_lattice_clusters (StructureAnalysis::buildClusters +
  // connectClusters port): resolves DXALatticeCorrespondence's per-atom labeling ambiguity into
  // globally-consistent clusters (contiguous grains sharing one discrete orientation labeling) and
  // the transitions between adjacent clusters (grain/defect boundaries) -- this is what
  // CrystalPathFinder (compute_dxa_crystal_path_edge_vectors, next stage) actually walks.
  struct DXALatticeClusters
  {
    ClusterGraph graph;

    std::vector<int32_t> atom_cluster;              // per atom: cluster id (0 = unresolved)
    std::vector<int32_t> atom_symmetry_permutation;  // per atom: index into dxa_lattice_structure(cluster.structure_type).permutations

    // Number of valid entries in this atom's own DXALatticeCorrespondence::neighbor_atom row. For
    // a classified atom this is just dxa_lattice_structure(its structure_type).num_neighbors. For
    // an UNCLASSIFIED atom (structure_type==LATTICE_OTHER) it starts at 0 and grows as
    // connectClusters appends classified neighbors that reference it -- exactly OVITO's own
    // mechanism (StructureAnalysis::connectClusters' "add this atom to the neighbor's own list of
    // neighbors" step) for letting CrystalPathFinder's reverse-neighbor-search branch step INTO a
    // disordered atom as a path intermediate, even though it has no native lattice correspondence
    // of its own.
    std::vector<int32_t> atom_neighbor_count;

    // Diagnostic: number of classified neighbors that couldn't be appended to an unclassified
    // atom's own row because it was already at capacity (DXA_MAX_NEIGHBORS) -- if nonzero, some
    // real physical neighbor relationships are being silently dropped, which CrystalPathFinder's
    // reverse-neighbor-search branch would otherwise have been able to use.
    long n_appends_dropped_diagnostic = 0;

    // The physical neighbor occupying raw stored slot `slot` NEVER changes after
    // DXALatticeCorrespondence assigns it -- use DXALatticeCorrespondence::neighbor(atom,slot)
    // directly for that. What atom_symmetry_permutation resolves is only the INTERPRETATION of
    // that same raw slot's ideal direction (see below), exactly OVITO's own split between
    // getNeighbor() [fixed at classification time] and neighborLatticeVector() [reinterpreted per
    // cluster] -- this is the mechanism that makes a whole cluster's slot labeling mutually
    // consistent despite each atom having found its own, independently arbitrary, matching
    // permutation during classification.

    // this atom's own ideal lattice vector for its RAW stored slot `slot` (same slot index as
    // DXALatticeCorrespondence::neighbor(atom,slot)), in its cluster's own shared local frame --
    // fixed reference table entry, no per-atom rotation involved (see lattice_structure.h file
    // header: this is the whole point of the discrete approach).
    inline Vec3d ideal_vector( size_t atom, int slot, const LatticeStructure& s ) const
    {
      return s.lattice_vectors[ s.permutations[ atom_symmetry_permutation[atom] ].permutation[slot] ];
    }
  };

  // Core algorithm behind compute_dxa_lattice_clusters -- StructureAnalysis::buildClusters() +
  // connectClusters() port. Deliberately grid/operator-independent (only needs each atom's already
  // -flattened real position) so it can be unit-tested standalone; the operator itself just
  // flattens grid positions and calls this. May append entries to lc.neighbor_atom (see
  // DXALatticeClusters::atom_neighbor_count's own comment). pos.size() must equal
  // lc.structure_type.size().
  void dxa_build_lattice_clusters( DXALatticeCorrespondence& lc, const std::vector<Vec3d>& pos, DXALatticeClusters& result );
}
