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
#include <ptm_constants.h>
#include <cmath>

// Quaternion/crystal-symmetry math shared by compute_grain_bond_misorientation.cu (GPU, per-bond
// disorientation) and grain_segmentation_algo.cpp (host, cluster-level running quaternion-sum
// accumulation during Node-Pair-Sampling merges). Ported VERBATIM from the vendored PTM library's
// own src/ptm/lib/ptm_quat.cpp (PM Larsen, MIT license -- the exact routine OVITO's own
// GrainSegmentationEngine calls), duplicated here rather than linked so the GPU pass can be
// genuinely CudaCompatible=true (the vendored ptm:: functions are plain host C++, not annotated for
// device compilation). See compute_grain_bond_misorientation.cu's own header comment for full scope
// notes (cubic 24-element / hcp "conventional" 12-element groups only, matching OVITO's own
// GrainSegmentation support -- ICO/diamond/graphene not supported there either).
namespace exaStamp
{
  ONIKA_HOST_DEVICE_FUNC inline void gb_quat_rot( const double* r, const double* a, double* b )
  {
    b[0] = r[0]*a[0] - r[1]*a[1] - r[2]*a[2] - r[3]*a[3];
    b[1] = r[0]*a[1] + r[1]*a[0] + r[2]*a[3] - r[3]*a[2];
    b[2] = r[0]*a[2] - r[1]*a[3] + r[2]*a[0] + r[3]*a[1];
    b[3] = r[0]*a[3] + r[1]*a[2] - r[2]*a[1] + r[3]*a[0];
  }

  ONIKA_HOST_DEVICE_FUNC inline void gb_generator_table( int structure_type, int & n_gen, double gen[24][4] )
  {
    if( structure_type == PTM_MATCH_FCC || structure_type == PTM_MATCH_BCC || structure_type == PTM_MATCH_SC )
    {
      n_gen = 24;
      const double g[24][4] = {
        {          1,          0,          0,          0 }, {  M_SQRT1_2,  M_SQRT1_2,          0,          0 },
        {  M_SQRT1_2,          0,  M_SQRT1_2,          0 }, {  M_SQRT1_2,          0,          0,  M_SQRT1_2 },
        {  M_SQRT1_2,          0,          0, -M_SQRT1_2 }, {  M_SQRT1_2,          0, -M_SQRT1_2,          0 },
        {  M_SQRT1_2, -M_SQRT1_2,          0,          0 }, {        0.5,        0.5,        0.5,        0.5 },
        {        0.5,        0.5,        0.5,       -0.5 }, {        0.5,        0.5,       -0.5,        0.5 },
        {        0.5,        0.5,       -0.5,       -0.5 }, {        0.5,       -0.5,        0.5,        0.5 },
        {        0.5,       -0.5,        0.5,       -0.5 }, {        0.5,       -0.5,       -0.5,        0.5 },
        {        0.5,       -0.5,       -0.5,       -0.5 }, {          0,          1,          0,          0 },
        {          0,  M_SQRT1_2,  M_SQRT1_2,          0 }, {          0,  M_SQRT1_2,          0,  M_SQRT1_2 },
        {          0,  M_SQRT1_2,          0, -M_SQRT1_2 }, {          0,  M_SQRT1_2, -M_SQRT1_2,          0 },
        {          0,          0,          1,          0 }, {          0,          0,  M_SQRT1_2,  M_SQRT1_2 },
        {          0,          0,  M_SQRT1_2, -M_SQRT1_2 }, {          0,          0,          0,          1 },
      };
      for(int i=0;i<24;i++) { gen[i][0]=g[i][0]; gen[i][1]=g[i][1]; gen[i][2]=g[i][2]; gen[i][3]=g[i][3]; }
    }
    else if( structure_type == PTM_MATCH_HCP )
    {
      n_gen = 12;
      const double SQRT3_2 = 0.86602540378443864676;
      const double g[12][4] = {
        {          1,          0,          0,          0 }, {    SQRT3_2,          0,          0,        0.5 },
        {    SQRT3_2,          0,          0,       -0.5 }, {        0.5,          0,          0,    SQRT3_2 },
        {        0.5,          0,          0,   -SQRT3_2 }, {          0,          1,          0,          0 },
        {          0,    SQRT3_2,        0.5,          0 }, {          0,    SQRT3_2,       -0.5,          0 },
        {          0,        0.5,    SQRT3_2,          0 }, {          0,        0.5,   -SQRT3_2,          0 },
        {          0,          0,          1,          0 }, {          0,          0,          0,          1 },
      };
      for(int i=0;i<12;i++) { gen[i][0]=g[i][0]; gen[i][1]=g[i][1]; gen[i][2]=g[i][2]; gen[i][3]=g[i][3]; }
    }
    else { n_gen = 0; }
  }

  // Rotates q into its symmetry-equivalent closest to the identity; returns the winning generator
  // index (needed by gb_map_quaternion_onto_target below), or -1 if unsupported structure_type.
  ONIKA_HOST_DEVICE_FUNC inline int gb_rotate_into_fundamental_zone( int structure_type, double* q )
  {
    double gen[24][4]; int n_gen = 0;
    gb_generator_table( structure_type, n_gen, gen );
    if( n_gen == 0 ) { return -1; }
    double best = 0.0; int bi = 0;
    for(int i=0;i<n_gen;i++)
    {
      const double t = fabs( q[0]*gen[i][0] - q[1]*gen[i][1] - q[2]*gen[i][2] - q[3]*gen[i][3] );
      if( t > best ) { best = t; bi = i; }
    }
    double f[4]; gb_quat_rot( q, gen[bi], f );
    q[0]=f[0]; q[1]=f[1]; q[2]=f[2]; q[3]=f[3];
    if( q[0] < 0.0 ) { q[0]=-q[0]; q[1]=-q[1]; q[2]=-q[2]; q[3]=-q[3]; }
    return bi;
  }

  // Returns the disorientation angle in DEGREES, or -1.0 if `structure_type` isn't one of the two
  // symmetry families supported here.
  ONIKA_HOST_DEVICE_FUNC inline double gb_disorientation_deg( int structure_type, const double* q0, const double* q1 )
  {
    double qinv[4] = { q0[0], -q0[1], -q0[2], -q0[3] };
    double qrot[4]; gb_quat_rot( qinv, q1, qrot );
    if( gb_rotate_into_fundamental_zone( structure_type, qrot ) < 0 ) { return -1.0; }
    double t = qrot[0];
    if( t > 1.0 ) t = 1.0; else if( t < -1.0 ) t = -1.0;
    return acos( 2.0*t*t - 1.0 ) * (180.0/M_PI);
  }

  // Ported verbatim from ptm_quat.cpp's own rotation_matrix_to_quaternion (Shepperd's method) --
  // needed because ptm_fields.cpp materializes PTM's orientation as a rotation-tensor grid field
  // (Mat3d), not the raw quaternion (which only exists in compute_ptm's own flat, MPI-ghost-unaware
  // output slot -- unusable here since neighbor bonds can cross into ghost territory).
  ONIKA_HOST_DEVICE_FUNC inline void gb_matrix_to_quat( const Mat3d& R, double* q )
  {
    const double r11=R.m11,r12=R.m12,r13=R.m13, r21=R.m21,r22=R.m22,r23=R.m23, r31=R.m31,r32=R.m32,r33=R.m33;
    q[0] = sqrt( fmax(0., (1.+r11+r22+r33)/4.) );
    q[1] = sqrt( fmax(0., (1.+r11-r22-r33)/4.) );
    q[2] = sqrt( fmax(0., (1.-r11+r22-r33)/4.) );
    q[3] = sqrt( fmax(0., (1.-r11-r22+r33)/4.) );
    int bi=0; double best=q[0];
    for(int i=1;i<4;i++) { if(q[i]>best) { best=q[i]; bi=i; } }
    if(bi==0)      { q[1]*=(r32-r23>=0.?1.:-1.); q[2]*=(r13-r31>=0.?1.:-1.); q[3]*=(r21-r12>=0.?1.:-1.); }
    else if(bi==1) { q[0]*=(r32-r23>=0.?1.:-1.); q[2]*=(r21+r12>=0.?1.:-1.); q[3]*=(r13+r31>=0.?1.:-1.); }
    else if(bi==2) { q[0]*=(r13-r31>=0.?1.:-1.); q[1]*=(r21+r12>=0.?1.:-1.); q[3]*=(r32+r23>=0.?1.:-1.); }
    else           { q[0]*=(r21-r12>=0.?1.:-1.); q[1]*=(r31+r13>=0.?1.:-1.); q[2]*=(r32+r23>=0.?1.:-1.); }
    const double n = sqrt( q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3] );
    q[0]/=n; q[1]/=n; q[2]/=n; q[3]/=n;
  }

  // OVITO's own cluster-level merge helper (GrainSegmentationEngine1::calculate_disorientation):
  // maps q's own symmetry-equivalent orientation closest to qtarget directly onto q (mutated in
  // place), returns the disorientation (degrees) between the ORIGINAL q and qtarget. Caller is
  // expected to then accumulate the mapped q (weighted by its own pre-mapping norm) into a running
  // cluster orientation sum -- see grain_segmentation_algo.cpp's own merge loop.
  inline double gb_map_quaternion_onto_target( int structure_type, const double* qtarget, double* q )
  {
    double qtemp[4]; double qtarget_inv[4] = { qtarget[0], -qtarget[1], -qtarget[2], -qtarget[3] };
    gb_quat_rot( qtarget_inv, q, qtemp );
    const int bi = gb_rotate_into_fundamental_zone( structure_type, qtemp );
    if( bi < 0 ) { return -1.0; }
    double gen[24][4]; int n_gen=0; gb_generator_table( structure_type, n_gen, gen );
    double mapped[4]; gb_quat_rot( q, gen[bi], mapped );
    if( mapped[0] < 0.0 ) { mapped[0]=-mapped[0]; mapped[1]=-mapped[1]; mapped[2]=-mapped[2]; mapped[3]=-mapped[3]; }
    // disorientation between original q and qtarget (same formula as gb_disorientation_deg, reusing
    // qtemp which is already q rotated into qtarget's own frame -- avoid a second quat_rot call)
    double t = qtemp[0];
    if( t > 1.0 ) t = 1.0; else if( t < -1.0 ) t = -1.0;
    const double disorientation_deg = acos( 2.0*t*t - 1.0 ) * (180.0/M_PI);
    q[0]=mapped[0]; q[1]=mapped[1]; q[2]=mapped[2]; q[3]=mapped[3];
    return disorientation_deg;
  }
}
