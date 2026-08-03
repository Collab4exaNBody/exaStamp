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

#include <onika/math/basic_types.h>
#include <vector>
#include <array>
#include <cstdint>

namespace exaStamp
{
  using namespace exanb;

  // Fixed, per-structure-type tables driving OVITO's real DXA elastic-mapping machinery
  // (ovito/src/ovito/crystalanalysis/modifier/dxa/StructureAnalysis.{h,cpp},
  // initializeListOfStructures() / determineLocalStructure()): a DISCRETE graph-topology match
  // (CNA bond-signature + neighbor-bond-graph isomorphism against these fixed reference tables),
  // not PTM's continuous least-squares/RMSD orientation fit. See src/delaunay/README.md for why
  // this replaces compute_ptm's role in the DXA pipeline specifically (PTM's continuous fit turned
  // out to have a materially different, wider strain tolerance than the discrete topology test
  // OVITO's real DXA source actually uses for the elastic mapping).
  //
  // OVITO splits "CoordinationStructure" (bond graph + CNA signature) from "LatticeStructure"
  // (adds a second neighbor shell, needed only for diamond's dual-shell CNA scheme) -- collapsed
  // into one struct here since diamond is out of scope (not tested, no test data), so for
  // FCC/HCP/BCC the two are identical 1:1 tables anyway.
  enum LatticeStructureType : int
  {
    LATTICE_OTHER = 0,
    LATTICE_FCC,
    LATTICE_HCP,
    LATTICE_BCC,
    NUM_LATTICE_TYPES
  };

  static constexpr int DXA_MAX_NEIGHBORS = 14; // BCC's own count, the largest of FCC(12)/HCP(12)/BCC(14)

  // One symmetry element of the structure's point group, found by brute-force permutation search
  // (see lattice_structure.cpp) over ways to relabel the fixed lattice_vectors table onto itself
  // via a single rigid rotation/reflection. permutation[slot] = which lattice_vectors[] entry
  // occupies canonical neighbor slot `slot` under this symmetry element (permutation[*][slot=0..)
  // is the identity mapping for the first entry, always present first).
  //
  // OVITO's own SymmetryPermutation also carries `product`/`inverseProduct` tables (composition of
  // symmetry elements), used only to re-align a cluster's orientation with a user-supplied
  // "preferred crystal orientation" list (an OVITO GUI feature). This pipeline has no such input
  // (_preferredCrystalOrientations stays empty, exactly like OVITO's own code path when it's
  // unset), so that whole step is a structural no-op here -- ponytail: dropped, add back only if
  // a preferred-orientation input is ever wired in.
  struct SymmetryPermutation
  {
    Mat3d transformation = onika::math::make_identity_matrix();
    std::vector<int> permutation;
  };

  struct LatticeStructure
  {
    LatticeStructureType type = LATTICE_OTHER;
    int num_neighbors = 0;
    std::vector<Vec3d> lattice_vectors; // fixed ideal-lattice-space neighbor directions (reference/unrotated frame, primitive-lattice-constant units)

    // per-slot-pair bond adjacency (symmetric, diagonal false), num_neighbors*num_neighbors flat.
    std::vector<uint8_t> bonded;
    inline bool is_bonded(int a, int b) const { return bonded[static_cast<size_t>(a) * num_neighbors + b] != 0; }

    // per-slot CNA family index: BCC 0=6-6-6 (1st shell, slots 0-7), 1=4-4-4 (2nd shell, slots
    // 8-13); FCC always 0 (4-2-1); HCP 0=4-2-1 (out-of-plane), 1=4-2-2 (in-plane, z==0).
    std::array<int, DXA_MAX_NEIGHBORS> cna_signature{};

    // two other neighbor slots, both bonded to `slot` and forming a non-coplanar triplet with it
    // -- used to fit the misorientation matrix between two adjacent atoms' own slot assignments
    // (buildClusters/connectClusters).
    std::array<std::array<int, 2>, DXA_MAX_NEIGHBORS> common_neighbors{};

    // the structure's point-group symmetry elements (identity first) -- resolves the discrete
    // labeling ambiguity left after a per-atom bond-topology match (many permutations can satisfy
    // the same local graph match; buildClusters picks the one consistent with already-visited
    // neighbors in the same cluster).
    std::vector<SymmetryPermutation> permutations;
  };

  // Lazily builds and caches the fixed structural tables on first call (thread-safe, std::call_once
  // internally). type must be LATTICE_FCC, LATTICE_HCP or LATTICE_BCC.
  const LatticeStructure& dxa_lattice_structure(LatticeStructureType type);

  // Sorts the [from,to) slice of a permutation into DESCENDING order using only the values
  // actually present in it (a bitmap over 0..max-1) -- makes it the lexicographically LAST
  // permutation of that suffix, so the next std::next_permutation() call jumps straight past every
  // remaining permutation sharing the current (rejected) prefix. Shared by lattice_structure.cpp's
  // own symmetry search and any per-atom neighbor-to-canonical-slot backtracking match that needs
  // the same permutation-search-with-pruning technique (ported from OVITO's bitmapSort()).
  inline void dxa_bitmap_sort_desc(std::vector<int>& v, int from, int to, int max)
  {
    unsigned long long bits = 0;
    for(int i=from; i<to; i++) { bits |= (1ull << v[i]); }
    int pos = from;
    for(int i=max-1; i>=0; i--) { if( bits & (1ull << i) ) { v[pos++] = i; } }
  }
}
