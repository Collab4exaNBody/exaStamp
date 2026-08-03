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

#include <exaStamp/delaunay/lattice_structure.h>

#include <algorithm>
#include <cmath>
#include <mutex>
#include <numeric>

// Table data (neighbor vectors, bond cutoffs, CNA signatures) and the symmetry-permutation search
// are ported verbatim from OVITO's real DXA source,
// ovito/src/ovito/crystalanalysis/modifier/dxa/StructureAnalysis.cpp,
// StructureAnalysis::initializeListOfStructures() -- see lattice_structure.h for why this exists.
namespace exaStamp
{
  using namespace exanb;

  namespace
  {
    inline bool vec3_close( const Vec3d& a, const Vec3d& b, double eps = 1e-6 )
    {
      return std::abs(a.x-b.x) <= eps && std::abs(a.y-b.y) <= eps && std::abs(a.z-b.z) <= eps;
    }

    inline bool mat3_close( const Mat3d& a, const Mat3d& b, double eps )
    {
      return std::abs(a.m11-b.m11)<=eps && std::abs(a.m12-b.m12)<=eps && std::abs(a.m13-b.m13)<=eps
          && std::abs(a.m21-b.m21)<=eps && std::abs(a.m22-b.m22)<=eps && std::abs(a.m23-b.m23)<=eps
          && std::abs(a.m31-b.m31)<=eps && std::abs(a.m32-b.m32)<=eps && std::abs(a.m33-b.m33)<=eps;
    }

    // Orthogonal iff M^T * M == Identity -- exactly OVITO's Matrix3::isOrthogonalMatrix() test
    // (accepts both proper rotations and reflections, same as the point-group symmetry elements
    // this search is looking for).
    inline bool mat3_is_orthogonal( const Mat3d& m, double eps = 1e-6 )
    {
      const Mat3d mtm = transpose(m) * m;
      return mat3_close( mtm, onika::math::make_identity_matrix(), eps );
    }

    Mat3d mat3_from_columns( const Vec3d& c0, const Vec3d& c1, const Vec3d& c2 )
    {
      return onika::math::make_mat3d( c0, c1, c2 );
    }

    // Finds two neighbor slots bonded to `slot`, both from `candidates`, that together with `slot`
    // form a non-coplanar (non-degenerate) triplet of directions. Ported from
    // initializeListOfStructures()'s "find two non-coplanar common neighbors for every neighbor
    // bond" loop.
    void find_common_neighbors( LatticeStructure& s )
    {
      for(int slot=0; slot<s.num_neighbors; slot++)
      {
        bool found = false;
        for(int i1=0; i1<s.num_neighbors && !found; i1++)
        {
          if( !s.is_bonded(slot,i1) ) continue;
          for(int i2=i1+1; i2<s.num_neighbors; i2++)
          {
            if( !s.is_bonded(slot,i2) ) continue;
            const Mat3d tm = mat3_from_columns( s.lattice_vectors[slot], s.lattice_vectors[i1], s.lattice_vectors[i2] );
            if( std::abs(determinant(tm)) > 1e-6 )
            {
              s.common_neighbors[slot] = { i1, i2 };
              found = true;
              break;
            }
          }
        }
      }
    }

    // Brute-force search (with pruning) for the structure's point-group symmetry elements: every
    // permutation `p` of {0..N-1} such that a single rigid rotation/reflection maps
    // lattice_vectors[i] -> lattice_vectors[p[i]] for every i simultaneously. Ported from
    // initializeListOfStructures()'s "Generate symmetry information" loop.
    void find_symmetry_permutations( LatticeStructure& s )
    {
      const int n = s.num_neighbors;

      // Three non-coplanar reference directions -- the rotation is fully determined by where
      // these three land.
      int nidx[3]; int found = 0;
      for(int i=0; i<n && found<3; i++)
      {
        // reuse whatever's already chosen plus the candidate, testing incrementally exactly like
        // the reference implementation (cross-product check at 2, determinant check at 3).
        if( found == 1 )
        {
          const Vec3d a = s.lattice_vectors[nidx[0]];
          const Vec3d b = s.lattice_vectors[i];
          if( norm2(cross(a,b)) <= 1e-6 ) { continue; }
        }
        else if( found == 2 )
        {
          const Mat3d tm = mat3_from_columns( s.lattice_vectors[nidx[0]], s.lattice_vectors[nidx[1]], s.lattice_vectors[i] );
          if( std::abs(determinant(tm)) <= 1e-6 ) { continue; }
        }
        nidx[found++] = i;
      }

      const Mat3d tm1 = mat3_from_columns( s.lattice_vectors[nidx[0]], s.lattice_vectors[nidx[1]], s.lattice_vectors[nidx[2]] );
      const Mat3d tm1inv = inverse(tm1);

      std::vector<int> permutation(n), last_permutation(n, -1);
      std::iota( permutation.begin(), permutation.end(), 0 );

      Mat3d transformation = onika::math::make_identity_matrix();
      do
      {
        int changed_from = 0;
        while( changed_from < n && permutation[changed_from] == last_permutation[changed_from] ) { ++changed_from; }
        last_permutation = permutation;

        if( changed_from <= nidx[2] )
        {
          const Mat3d tm2 = mat3_from_columns( s.lattice_vectors[permutation[nidx[0]]], s.lattice_vectors[permutation[nidx[1]]], s.lattice_vectors[permutation[nidx[2]]] );
          transformation = tm2 * tm1inv;
          if( !mat3_is_orthogonal(transformation) )
          {
            dxa_bitmap_sort_desc( permutation, nidx[2]+1, n, n );
            continue;
          }
          changed_from = 0;
        }

        int sort_from = nidx[2];
        int invalid_from = changed_from;
        for(; invalid_from<n; invalid_from++)
        {
          const Vec3d v = transformation * s.lattice_vectors[invalid_from];
          if( !vec3_close( v, s.lattice_vectors[permutation[invalid_from]] ) ) { break; }
        }

        if( invalid_from == n )
        {
          SymmetryPermutation sp;
          sp.transformation = transformation;
          sp.permutation.assign( permutation.begin(), permutation.begin()+n );
          s.permutations.push_back( std::move(sp) );
        }
        else
        {
          sort_from = invalid_from;
        }
        dxa_bitmap_sort_desc( permutation, sort_from+1, n, n );
      }
      while( std::next_permutation( permutation.begin(), permutation.end() ) );
    }

    LatticeStructure build_fcc()
    {
      LatticeStructure s; s.type = LATTICE_FCC; s.num_neighbors = 12;
      s.lattice_vectors = {
        { 0.5,  0.5,  0.0}, { 0.0,  0.5,  0.5}, { 0.5,  0.0,  0.5},
        {-0.5, -0.5,  0.0}, { 0.0, -0.5, -0.5}, {-0.5,  0.0, -0.5},
        {-0.5,  0.5,  0.0}, { 0.0, -0.5,  0.5}, {-0.5,  0.0,  0.5},
        { 0.5, -0.5,  0.0}, { 0.0,  0.5, -0.5}, { 0.5,  0.0, -0.5}
      };
      s.bonded.assign( 12*12, 0 );
      for(int a=0;a<12;a++) for(int b=a+1;b<12;b++)
      {
        const bool bonded = norm( s.lattice_vectors[a] - s.lattice_vectors[b] ) < (std::sqrt(0.5)+1.0)*0.5;
        s.bonded[a*12+b] = s.bonded[b*12+a] = bonded ? 1 : 0;
      }
      for(int i=0;i<12;i++) { s.cna_signature[i] = 0; }
      find_common_neighbors(s);
      find_symmetry_permutations(s);
      return s;
    }

    LatticeStructure build_hcp()
    {
      LatticeStructure s; s.type = LATTICE_HCP; s.num_neighbors = 12;
      const double s2 = std::sqrt(2.0), s6 = std::sqrt(6.0), s3 = std::sqrt(3.0);
      s.lattice_vectors = {
        { s2/4.0, -s6/4.0,  0.0}, {-s2/2.0,     0.0,  0.0}, {-s2/4.0,  s6/12.0, -s3/3.0},
        { s2/4.0,  s6/12.0, -s3/3.0}, { 0.0, -s6/6.0, -s3/3.0}, {-s2/4.0,  s6/4.0,  0.0},
        { s2/4.0,  s6/4.0,  0.0}, { s2/2.0,     0.0,  0.0}, {-s2/4.0, -s6/4.0,  0.0},
        { 0.0, -s6/6.0,  s3/3.0}, { s2/4.0,  s6/12.0,  s3/3.0}, {-s2/4.0,  s6/12.0,  s3/3.0}
      };
      s.bonded.assign( 12*12, 0 );
      for(int a=0;a<12;a++) for(int b=a+1;b<12;b++)
      {
        const bool bonded = norm( s.lattice_vectors[a] - s.lattice_vectors[b] ) < (std::sqrt(0.5)+1.0)*0.5;
        s.bonded[a*12+b] = s.bonded[b*12+a] = bonded ? 1 : 0;
      }
      for(int i=0;i<12;i++) { s.cna_signature[i] = ( s.lattice_vectors[i].z == 0.0 ) ? 1 : 0; }
      find_common_neighbors(s);
      find_symmetry_permutations(s);
      return s;
    }

    LatticeStructure build_bcc()
    {
      LatticeStructure s; s.type = LATTICE_BCC; s.num_neighbors = 14;
      s.lattice_vectors = {
        { 0.5,  0.5,  0.5}, {-0.5,  0.5,  0.5}, { 0.5,  0.5, -0.5}, {-0.5, -0.5,  0.5},
        { 0.5, -0.5,  0.5}, {-0.5,  0.5, -0.5}, {-0.5, -0.5, -0.5}, { 0.5, -0.5, -0.5},
        { 1.0,  0.0,  0.0}, {-1.0,  0.0,  0.0}, { 0.0,  1.0,  0.0}, { 0.0, -1.0,  0.0},
        { 0.0,  0.0,  1.0}, { 0.0,  0.0, -1.0}
      };
      s.bonded.assign( 14*14, 0 );
      for(int a=0;a<14;a++) for(int b=a+1;b<14;b++)
      {
        const bool bonded = norm( s.lattice_vectors[a] - s.lattice_vectors[b] ) < (1.0+std::sqrt(2.0))*0.5;
        s.bonded[a*14+b] = s.bonded[b*14+a] = bonded ? 1 : 0;
      }
      for(int i=0;i<14;i++) { s.cna_signature[i] = (i<8) ? 0 : 1; }
      find_common_neighbors(s);
      find_symmetry_permutations(s);
      return s;
    }
  }

  const LatticeStructure& dxa_lattice_structure( LatticeStructureType type )
  {
    static LatticeStructure tables[NUM_LATTICE_TYPES];
    static std::once_flag once;
    std::call_once( once, [&]()
    {
      tables[LATTICE_FCC] = build_fcc();
      tables[LATTICE_HCP] = build_hcp();
      tables[LATTICE_BCC] = build_bcc();
    });
    return tables[type];
  }
}
