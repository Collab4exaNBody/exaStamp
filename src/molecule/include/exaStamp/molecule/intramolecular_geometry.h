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

#include <cmath>
#include <onika/math/basic_types.h>
#include <onika/cuda/cuda.h>

namespace exaStamp
{

  // Bend angle theta between r1 = ra - rb and r2 = rc - rb, and its gradients with respect to r1 and r2.
  // theta is computed with atan2, which is accurate on the whole [0,pi] range (acos is not near 0 and pi).
  // The gradients are proportional to 1/sin(theta). sin(theta) is bounded below by bend_min_sin so that
  // (nearly) linear angles give large but finite forces instead of NaNs.
  struct BendGeometry
  {
    double theta = 0.0;
    Vec3d dtheta_dr1 = {0.,0.,0.};
    Vec3d dtheta_dr2 = {0.,0.,0.};
  };

  static constexpr double bend_min_sin = 1.0e-3;

  ONIKA_HOST_DEVICE_FUNC
  inline BendGeometry bend_geometry( const Vec3d& r1, const Vec3d& r2 )
  {
    const double n1sq = dot(r1,r1);
    const double n2sq = dot(r2,r2);
    const double n1n2 = sqrt( n1sq * n2sq );
    const double r1_dot_r2 = dot(r1,r2);
    const double r1_cross_r2 = norm( cross(r1,r2) );
    const double cos_theta = r1_dot_r2 / n1n2;
    const double raw_sin_theta = r1_cross_r2 / n1n2;
    const double sin_theta = ( raw_sin_theta > bend_min_sin ) ? raw_sin_theta : bend_min_sin;
    // dtheta/dr = -1/sin(theta) * dcos(theta)/dr
    return { atan2( r1_cross_r2 , r1_dot_r2 )
           , ( r1 * (cos_theta/n1sq) - r2 / n1n2 ) / sin_theta
           , ( r2 * (cos_theta/n2sq) - r1 / n1n2 ) / sin_theta };
  }

  // Dihedral angle phi of the chain a-b-c-d, from r10 = ra - rb, r12 = rc - rb, r23 = rd - rc.
  // phi is in [-pi,pi], computed with atan2 (exact sign, including phi = pi for a planar trans chain).
  // dphi_dra and dphi_drd are the gradients of phi with respect to the positions of the end atoms a and d.
  // If a-b-c or b-c-d is (numerically) collinear, phi is undefined : valid is false and gradients are zero.
  struct TorsionGeometry
  {
    double phi = 0.0;
    Vec3d dphi_dra = {0.,0.,0.};
    Vec3d dphi_drd = {0.,0.,0.};
    bool valid = false;
  };

  static constexpr double torsion_min_sin2 = 1.0e-20;

  ONIKA_HOST_DEVICE_FUNC
  inline TorsionGeometry torsion_geometry( const Vec3d& r10, const Vec3d& r12, const Vec3d& r23 )
  {
    const Vec3d V = cross(r10,r12);
    const Vec3d W = cross(r23,r12);
    const double V2 = dot(V,V);
    const double W2 = dot(W,W);
    const double r12sq = dot(r12,r12);
    const double norm_r12 = sqrt( r12sq );
    const double phi = atan2( dot(cross(V,W),r12) / norm_r12 , dot(V,W) );
    if( V2 <= torsion_min_sin2 * dot(r10,r10) * r12sq || W2 <= torsion_min_sin2 * dot(r23,r23) * r12sq )
    {
      return { phi , {0.,0.,0.} , {0.,0.,0.} , false };
    }
    return { phi , V * ( norm_r12 / V2 ) , W * ( - norm_r12 / W2 ) , true };
  }

}
