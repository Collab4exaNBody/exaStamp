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

#include <exaStamp/delaunay/dxa_lattice_correspondence.h>

#include <algorithm>
#include <cmath>
#include <deque>

// DXA elastic-mapping step, resolving compute_dxa_lattice_correspondence's per-atom labeling
// ambiguity into globally-consistent clusters -- port of StructureAnalysis::buildClusters() +
// connectClusters() (ovito/src/ovito/crystalanalysis/modifier/dxa/StructureAnalysis.cpp). See
// compute_dxa_lattice_correspondence.cpp's own file header for why this two-stage design exists.
//
// buildClusters: BFS from each not-yet-visited classified atom, growing a cluster by testing, for
// every unassigned same-structure-type neighbor, whether ITS OWN (independently, arbitrarily
// chosen) neighbor-to-slot labeling is related to the current atom's by one of the structure's
// known point-group symmetry elements -- fit via three non-coplanar reference atoms shared between
// the two, exactly as OVITO does it. If a match is found, the neighbor joins the cluster with that
// symmetry element recorded (atom_symmetry_permutation); if not, growth simply doesn't cross that
// bond (a real grain/defect boundary). Also accumulates a least-squares cluster orientation
// (display/diagnostic only, not consumed by the elastic mapping itself).
//
// connectClusters: for every bond crossing from one already-assigned cluster into a DIFFERENT one,
// fits the same three-non-coplanar-atom misorientation matrix (this time between two different
// clusters' own frames, not within one) and records it as a ClusterTransition -- this is what lets
// CrystalPathFinder (next stage) route a path across a low-angle sub-boundary or between grains.
// Also implements OVITO's "extend an unclassified atom's own neighbor list" step: when a classified
// atom borders a genuinely unclassified one (structure_type==LATTICE_OTHER), that unclassified
// atom's own DXALatticeCorrespondence::neighbor_atom row gets the classified atom appended to it --
// with no native lattice correspondence of its own, this is the only way a later path search can
// step INTO it as an intermediate (CrystalPathFinder's own "reverse neighbor search" branch).
//
// Deliberately grid/operator-independent (see dxa_lattice_correspondence.h) so it's unit-testable
// standalone -- see /tmp scratch test used to validate this port before wiring it into a real
// simulation pipeline.
namespace exaStamp
{
  using namespace exanb;

  namespace
  {
    inline bool mat3_close( const Mat3d& a, const Mat3d& b, double eps )
    {
      return std::abs(a.m11-b.m11)<=eps && std::abs(a.m12-b.m12)<=eps && std::abs(a.m13-b.m13)<=eps
          && std::abs(a.m21-b.m21)<=eps && std::abs(a.m22-b.m22)<=eps && std::abs(a.m23-b.m23)<=eps
          && std::abs(a.m31-b.m31)<=eps && std::abs(a.m32-b.m32)<=eps && std::abs(a.m33-b.m33)<=eps;
    }

    // Linear scan of an atom's own (possibly appended-beyond-native-count) neighbor row for a
    // specific physical neighbor -- mirrors StructureAnalysis::findNeighbor().
    inline int find_neighbor_slot( const DXALatticeCorrespondence& lc, const DXALatticeClusters& lcl, int64_t atom, int64_t target )
    {
      const int count = lcl.atom_neighbor_count[atom];
      for(int j=0;j<count;j++) { if( lc.neighbor(atom,j) == target ) { return j; } }
      return -1;
    }
  }

  void dxa_build_lattice_clusters( DXALatticeCorrespondence& lc, const std::vector<Vec3d>& pos, const std::vector<uint64_t>& global_id, DXALatticeClusters& result )
  {
    const size_t n_particles = lc.structure_type.size();

    result.graph = ClusterGraph{};
    result.atom_cluster.assign( n_particles, 0 );
    result.atom_symmetry_permutation.assign( n_particles, 0 );
    result.atom_neighbor_count.resize( n_particles );
    for(size_t a=0;a<n_particles;a++)
    {
      const auto st = static_cast<LatticeStructureType>( lc.structure_type[a] );
      result.atom_neighbor_count[a] = ( st == LATTICE_OTHER ) ? 0 : dxa_lattice_structure(st).num_neighbors;
    }

    // --- buildClusters ---
    // Seed order is sorted by global_id, NOT local array index -- see this function's own
    // declaration comment (dxa_lattice_correspondence.h) for why: the local index is an arbitrary
    // accident of this rank's own particle layout, and basing the seed (hence the arbitrary
    // reference orientation each new cluster grows from) on it let independent computations of the
    // same physical region disagree near a real defect.
    std::vector<size_t> seed_order( n_particles );
    for(size_t a=0;a<n_particles;a++) { seed_order[a] = a; }
    std::sort( seed_order.begin(), seed_order.end(), [&]( size_t a, size_t b ) { return global_id[a] < global_id[b]; } );

    for(size_t seed_idx : seed_order)
    {
      const size_t seed = seed_idx;
      if( result.atom_cluster[seed] != 0 ) { continue; }
      const auto st = static_cast<LatticeStructureType>( lc.structure_type[seed] );
      if( st == LATTICE_OTHER ) { continue; }

      const int cluster_id = result.graph.create_cluster( static_cast<int>(st) );
      result.atom_cluster[seed] = cluster_id;
      result.atom_symmetry_permutation[seed] = 0;
      result.graph.cluster(cluster_id).atom_count = 1;

      Mat3d orientationV = onika::math::make_zero_matrix();
      Mat3d orientationW = onika::math::make_zero_matrix();

      std::deque<int64_t> to_visit(1, static_cast<int64_t>(seed));
      while( !to_visit.empty() )
      {
        const int64_t cur = to_visit.front(); to_visit.pop_front();
        const LatticeStructure& ls = dxa_lattice_structure( static_cast<LatticeStructureType>( lc.structure_type[cur] ) );
        const auto& perm = ls.permutations[ result.atom_symmetry_permutation[cur] ].permutation;

        for(int slot=0; slot<ls.num_neighbors; slot++)
        {
          const int64_t nbr = lc.neighbor( cur, slot );
          if( nbr < 0 ) { continue; }

          const Vec3d latticeVector = ls.lattice_vectors[ perm[slot] ];
          const Vec3d spatialVector = pos[nbr] - pos[cur];
          orientationV = orientationV + tensor( latticeVector, latticeVector );
          orientationW = orientationW + tensor( spatialVector, latticeVector );

          if( result.atom_cluster[nbr] != 0 ) { continue; }
          if( lc.structure_type[nbr] != lc.structure_type[cur] ) { continue; }

          // Three non-coplanar reference atoms shared between cur and nbr: two common-neighbor
          // slots plus cur itself -- exactly OVITO's own buildClusters() fit.
          Mat3d tm1{}, tm2{};
          {
            Vec3d col1[3], col2[3];
            bool proper = true;
            for(int i=0;i<3 && proper;i++)
            {
              int64_t atomIndex;
              if( i != 2 )
              {
                const int cslot = ls.common_neighbors[slot][i];
                atomIndex = lc.neighbor( cur, cslot );
                col1[i] = ls.lattice_vectors[ perm[cslot] ] - ls.lattice_vectors[ perm[slot] ];
              }
              else
              {
                atomIndex = cur;
                col1[i] = Vec3d{0,0,0} - ls.lattice_vectors[ perm[slot] ];
              }
              if( atomIndex < 0 ) { proper = false; break; }
              const int j = find_neighbor_slot( lc, result, nbr, atomIndex );
              if( j < 0 ) { proper = false; break; }
              col2[i] = ls.lattice_vectors[j]; // nbr's own RAW slot order -- it has no symmetry permutation assigned yet
            }
            if( !proper ) { continue; }
            tm1 = onika::math::make_mat3d( col1[0], col1[1], col1[2] );
            tm2 = onika::math::make_mat3d( col2[0], col2[1], col2[2] );
          }

          if( std::abs(determinant(tm1)) <= 1e-9 ) { continue; }
          const Mat3d tm2inv_check = inverse(tm2);
          if( !std::isfinite(tm2inv_check.m11) ) { continue; }
          const Mat3d transition = tm1 * tm2inv_check;

          for(size_t si=0; si<ls.permutations.size(); si++)
          {
            if( mat3_close( transition, ls.permutations[si].transformation, 1e-4 ) )
            {
              result.atom_cluster[nbr] = cluster_id;
              result.atom_symmetry_permutation[nbr] = static_cast<int>(si);
              result.graph.cluster(cluster_id).atom_count++;
              to_visit.push_back( nbr );
              break;
            }
          }
        }
      }

      const Mat3d orientationVinv = inverse(orientationV);
      result.graph.cluster(cluster_id).orientation = orientationW * orientationVinv;
    }

    // --- connectClusters, part 1: extend unclassified neighbors' own rows ---
    long n_appends_dropped = 0; // DIAGNOSTIC: classified neighbors that couldn't be appended because the unclassified atom's own row was already full (capacity DXA_MAX_NEIGHBORS)
    for(size_t a=0; a<n_particles; a++)
    {
      const int cluster1_id = result.atom_cluster[a];
      if( cluster1_id == 0 ) { continue; }

      const LatticeStructure& ls = dxa_lattice_structure( static_cast<LatticeStructureType>( lc.structure_type[a] ) );

      for(int slot=0; slot<ls.num_neighbors; slot++)
      {
        const int64_t nbr = lc.neighbor( a, slot );
        if( nbr < 0 ) { continue; }
        if( result.atom_cluster[nbr] != 0 ) { continue; } // only unclassified neighbors get extended

        const int count = result.atom_neighbor_count[nbr];
        if( count < DXA_MAX_NEIGHBORS )
        {
          lc.neighbor_atom[ nbr*DXA_MAX_NEIGHBORS + count ] = static_cast<int64_t>(a);
          result.atom_neighbor_count[nbr] = count + 1;
        }
        else
        {
          ++n_appends_dropped;
        }
      }
    }
    result.n_appends_dropped_diagnostic = n_appends_dropped;

    // --- connectClusters, part 2: fit transitions between distinct clusters ---
    for(size_t a=0; a<n_particles; a++)
    {
      const int cluster1_id = result.atom_cluster[a];
      if( cluster1_id == 0 ) { continue; }
      const LatticeStructure& ls = dxa_lattice_structure( static_cast<LatticeStructureType>( lc.structure_type[a] ) );
      const auto& perm = ls.permutations[ result.atom_symmetry_permutation[a] ].permutation;

      for(int slot=0; slot<ls.num_neighbors; slot++)
      {
        const int64_t nbr = lc.neighbor( a, slot );
        if( nbr < 0 ) { continue; }
        const int cluster2_id = result.atom_cluster[nbr];
        if( cluster2_id == 0 || cluster2_id == cluster1_id ) { continue; }

        // Skip if a transition already exists between the two clusters (just bump its area).
        bool found_existing = false;
        for( int ti : result.graph.cluster(cluster1_id).transition_indices )
        {
          if( result.graph.transitions[ti].cluster2 == cluster2_id )
          {
            result.graph.transitions[ti].area++;
            result.graph.transitions[ result.graph.transitions[ti].reverse_index ].area++;
            found_existing = true;
            break;
          }
        }
        if( found_existing ) { continue; }

        const LatticeStructure& ls2 = dxa_lattice_structure( static_cast<LatticeStructureType>( lc.structure_type[nbr] ) );
        const auto& perm2 = ls2.permutations[ result.atom_symmetry_permutation[nbr] ].permutation;

        Mat3d tm1{}, tm2{};
        {
          Vec3d col1[3], col2[3];
          bool proper = true;
          for(int i=0;i<3 && proper;i++)
          {
            int64_t atomIndex;
            if( i != 2 )
            {
              const int cslot = ls.common_neighbors[slot][i];
              atomIndex = lc.neighbor( a, cslot );
              col1[i] = ls.lattice_vectors[ perm[cslot] ] - ls.lattice_vectors[ perm[slot] ];
            }
            else
            {
              atomIndex = a;
              col1[i] = Vec3d{0,0,0} - ls.lattice_vectors[ perm[slot] ];
            }
            if( atomIndex < 0 ) { proper = false; break; }
            const int j = find_neighbor_slot( lc, result, nbr, atomIndex );
            if( j < 0 ) { proper = false; break; }
            col2[i] = ls2.lattice_vectors[ perm2[j] ]; // nbr's OWN cluster frame this time (already assigned)
          }
          if( !proper ) { continue; }
          tm1 = onika::math::make_mat3d( col1[0], col1[1], col1[2] );
          tm2 = onika::math::make_mat3d( col2[0], col2[1], col2[2] );
        }

        if( std::abs(determinant(tm1)) <= 1e-9 ) { continue; }
        const Mat3d tm1inv = inverse(tm1);
        if( !std::isfinite(tm1inv.m11) ) { continue; }
        const Mat3d transition = tm2 * tm1inv;

        const Mat3d check = transpose(transition) * transition;
        if( mat3_close( check, onika::math::make_identity_matrix(), 1e-4 ) )
        {
          const int ti = result.graph.create_transition( cluster1_id, cluster2_id, transition, 1 );
          result.graph.transitions[ti].area++;
          result.graph.transitions[ result.graph.transitions[ti].reverse_index ].area++;
        }
      }
    }
  }
}
