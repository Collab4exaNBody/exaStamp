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
#include <onika/cuda/cuda.h>

namespace exaStamp
{

  ONIKA_HOST_DEVICE_FUNC
  inline Vec3d periodic_r_delta(const Vec3d& r1, const Vec3d& r2, const Vec3d& size_box, double half_min_size_box)
  {
    Vec3d rij = r2 - r1;

    if(rij.x>  half_min_size_box) rij.x -= size_box.x;
    if(rij.x< -half_min_size_box) rij.x += size_box.x;
    if(rij.y>  half_min_size_box) rij.y -= size_box.y;
    if(rij.y< -half_min_size_box) rij.y += size_box.y;
    if(rij.z>  half_min_size_box) rij.z -= size_box.z;
    if(rij.z< -half_min_size_box) rij.z += size_box.z;
    assert(norm(rij)<half_min_size_box);
    assert( (rij != Vec3d{0,0,0}) );
    
    return rij;
  }

  ONIKA_HOST_DEVICE_FUNC
  inline Vec3d periodic_r_delta_loop(const Vec3d& r1, const Vec3d& r2, const Vec3d& size_box, double half_min_size_box)
  {
    Vec3d rij = r2 - r1;
    while(rij.x>  half_min_size_box) rij.x -= size_box.x;
    while(rij.x< -half_min_size_box) rij.x += size_box.x;
    while(rij.y>  half_min_size_box) rij.y -= size_box.y;
    while(rij.y< -half_min_size_box) rij.y += size_box.y;
    while(rij.z>  half_min_size_box) rij.z -= size_box.z;
    while(rij.z< -half_min_size_box) rij.z += size_box.z;
    assert(norm(rij)<half_min_size_box);
    assert( (rij != Vec3d{0,0,0}) );
    return rij;
  }

}
