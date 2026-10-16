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

#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/cuda/cuda.h>
#include <exanb/core/domain.h>

namespace exaStamp
{
  using onika::math::Vec3d;
  using onika::math::Mat3d;

  // Periodic images of the atoms of an intramolecular term (bond, angle, torsion, improper, pair).
  //
  // Forces of a term are applied to the owned copy of each atom (bonded forces are never folded back from ghosts),
  // so a term references the owned copy whenever it exists. When a molecule crosses a periodic boundary, these
  // copies may be one box length apart : the image shift of each atom, in box lengths along each grid axis, is
  // computed when term lists are rebuilt, so that atom k of the term is at r_k + n_k * L, with n_0 = 0.
  // Positions evolve continuously between two rebuilds, so the shifts stay valid until the next one.
  //
  // Encoding : 3 bits per axis (shift + 4, shift in [-4,3]), 9 bits per atom for atoms 1 to 3 -> 27 bits.
  using IntramolecularImage = uint32_t;

  struct IntramolecularImageCodec
  {
    static constexpr int BITS = 3;
    static constexpr int OFFSET = 4;
    static constexpr int MIN_SHIFT = -OFFSET;
    static constexpr int MAX_SHIFT = (1<<BITS) - 1 - OFFSET;
    static constexpr int MAX_ATOMS = 4;
  };

  // displacement (grid space) to add to the stored position of atom k of a term
  ONIKA_HOST_DEVICE_FUNC inline Vec3d intramolecular_image_shift( IntramolecularImage img, int k, const Vec3d& size_box )
  {
    using C = IntramolecularImageCodec;
    if( k == 0 || img == 0 ) return Vec3d{0.,0.,0.};
    const unsigned int bits = img >> ( (k-1) * 3 * C::BITS );
    constexpr unsigned int mask = (1u<<C::BITS) - 1;
    return Vec3d{ double( int( bits                & mask) - C::OFFSET ) * size_box.x
                , double( int((bits >>   C::BITS ) & mask) - C::OFFSET ) * size_box.y
                , double( int((bits >> 2*C::BITS ) & mask) - C::OFFSET ) * size_box.z };
  }

  // vector from atom a to atom b of a term (grid space), with stored positions ra, rb
  ONIKA_HOST_DEVICE_FUNC inline Vec3d intramolecular_delta( const Vec3d& ra, int a, const Vec3d& rb, int b, IntramolecularImage img, const Vec3d& size_box )
  {
    if( img == 0 ) return rb - ra;
    return ( rb + intramolecular_image_shift(img,b,size_box) ) - ( ra + intramolecular_image_shift(img,a,size_box) );
  }

  // host side : computes image shifts for the atoms of a term when term lists are rebuilt
  struct IntramolecularImageBuilder
  {
    Vec3d m_size_box = {0.,0.,0.};
    Mat3d m_xform = {};
    bool m_xform_is_identity = true;
    bool m_periodic[3] = {false,false,false};

    inline IntramolecularImageBuilder() = default;
    inline IntramolecularImageBuilder( const exanb::Domain& domain )
      : m_size_box( domain.extent() - domain.origin() )
      , m_xform( domain.xform() )
      , m_xform_is_identity( domain.xform_is_identity() )
      , m_periodic{ domain.periodic_boundary_x() , domain.periodic_boundary_y() , domain.periodic_boundary_z() }
    {}

    // shift (in box lengths) bringing d = r_k - r_anchor to its shortest physical length
    inline void closest_image( const Vec3d& d, int n[3] ) const
    {
      const double dd[3] = { d.x , d.y , d.z };
      const double L[3] = { m_size_box.x , m_size_box.y , m_size_box.z };
      for(int i=0;i<3;i++) n[i] = m_periodic[i] ? - int( std::round( dd[i] / L[i] ) ) : 0;
      if( m_xform_is_identity ) return;
      // a sheared cell : the shortest vector may differ by one more cell along each periodic axis
      int best[3] = { n[0], n[1], n[2] };
      double best_d2 = -1.0;
      for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++) for(int k=-1;k<=1;k++)
      {
        const int c[3] = { n[0] + i , n[1] + j , n[2] + k };
        if( (i!=0 && !m_periodic[0]) || (j!=0 && !m_periodic[1]) || (k!=0 && !m_periodic[2]) ) continue;
        const Vec3d v = m_xform * Vec3d{ d.x + c[0]*L[0] , d.y + c[1]*L[1] , d.z + c[2]*L[2] };
        const double d2 = dot(v,v);
        if( best_d2 < 0.0 || d2 < best_d2 ) { best_d2 = d2; best[0]=c[0]; best[1]=c[1]; best[2]=c[2]; }
      }
      n[0]=best[0]; n[1]=best[1]; n[2]=best[2];
    }

    // closest periodic image of the grid space vector d, any number of box lengths away
    inline Vec3d closest_delta( const Vec3d& d ) const
    {
      int n[3];
      closest_image( d , n );
      return Vec3d{ d.x + n[0]*m_size_box.x , d.y + n[1]*m_size_box.y , d.z + n[2]*m_size_box.z };
    }

    // physical length of a grid space vector
    inline double physical_norm( const Vec3d& d ) const
    {
      return m_xform_is_identity ? norm(d) : norm( m_xform * d );
    }

    // the closest image is unique for physical distances below half of the smallest width of the cell
    // (distance between opposite faces) along periodic directions
    inline double max_unique_distance() const
    {
      const Vec3d a = m_xform * Vec3d{ m_size_box.x , 0. , 0. };
      const Vec3d b = m_xform * Vec3d{ 0. , m_size_box.y , 0. };
      const Vec3d c = m_xform * Vec3d{ 0. , 0. , m_size_box.z };
      const double volume = std::abs( dot( a , cross(b,c) ) );
      double w = std::numeric_limits<double>::max();
      if( m_periodic[0] ) w = std::min( w , volume / norm(cross(b,c)) );
      if( m_periodic[1] ) w = std::min( w , volume / norm(cross(c,a)) );
      if( m_periodic[2] ) w = std::min( w , volume / norm(cross(a,b)) );
      return 0.5 * w;
    }

    // r : stored positions of the n_atoms atoms of the term (grid space)
    // anchor[k] : atom (k' < k) atom k is chemically bonded to (or close to), atom k is placed at its closest image from it.
    // returns false if a shift does not fit in the encoding
    inline bool build( const Vec3d* r, const int* anchor, int n_atoms, IntramolecularImage& img ) const
    {
      using C = IntramolecularImageCodec;
      int shift[C::MAX_ATOMS][3] = { {0,0,0} };
      img = 0;
      for(int k=1;k<n_atoms;k++)
      {
        const int a = anchor[k];
        const Vec3d ra = r[a] + Vec3d{ shift[a][0]*m_size_box.x , shift[a][1]*m_size_box.y , shift[a][2]*m_size_box.z };
        // r[k] + n*L is the image of atom k closest to the (shifted) anchor : n is the shift of atom k
        int n[3];
        closest_image( r[k] - ra , n );
        for(int i=0;i<3;i++)
        {
          shift[k][i] = n[i];
          if( shift[k][i] < C::MIN_SHIFT || shift[k][i] > C::MAX_SHIFT ) return false;
          img |= IntramolecularImage( shift[k][i] + C::OFFSET ) << ( ( (k-1) * 3 + i ) * C::BITS );
        }
      }
      // all shifts zero : encode as 0 so that kernels can skip the shift computation
      bool all_zero = true;
      for(int k=1;k<n_atoms;k++) for(int i=0;i<3;i++) if( shift[k][i] != 0 ) all_zero = false;
      if( all_zero ) img = 0;
      return true;
    }
  };

}
