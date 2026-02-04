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

#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types.h>
#include <exaStamp/particle_species/particle_specie.h>

namespace exaStamp
{

  struct KineticEnergyFunctor
  {
    const ParticleSpecie* __restrict__ m_species = nullptr;
    const unsigned int m_type = 0;    
    ONIKA_HOST_DEVICE_FUNC inline double operator () (double vx, double vy, double vz, unsigned int type) const
    {
      return m_species[type].m_mass * (vx*vx+vy*vy+vz*vz);
    }
    ONIKA_HOST_DEVICE_FUNC inline double operator () (double vx, double vy, double vz) const
    {
      return this->operator() (vx,vy,vz,m_type);
    }
  };

  struct KineticEnergyTensorFunctor
  {
    const ParticleSpecie* __restrict__ m_species = nullptr;
    const unsigned int m_type = 0;    
    ONIKA_HOST_DEVICE_FUNC inline Mat3d operator () (double vx, double vy, double vz, unsigned int type) const
    {
      return Mat3d{vx*vx,vx*vy,vx*vz,vy*vx,vy*vy,vy*vz,vz*vx,vz*vy,vz*vz} * m_species[type].m_mass;
    }
    ONIKA_HOST_DEVICE_FUNC inline Mat3d operator () (double vx, double vy, double vz) const
    {
      return this->operator() (vx,vy,vz,m_type);
    }
  };

  struct MomentumFunctor
  {
    const ParticleSpecie* __restrict__ m_species = nullptr;
    const unsigned int m_type = 0;    
    ONIKA_HOST_DEVICE_FUNC inline Vec3d operator () (double vx, double vy, double vz, unsigned int type) const
    {
      return Vec3d{vx,vy,vz} * m_species[type].m_mass;
    }
    ONIKA_HOST_DEVICE_FUNC inline Vec3d operator () (double vx, double vy, double vz) const
    {
      return this->operator() (vx,vy,vz,m_type);
    }
  };

  struct MassFunctor
  {
    const ParticleSpecie* __restrict__ m_species = nullptr;
    const unsigned int m_type = 0;    
    ONIKA_HOST_DEVICE_FUNC inline double operator () (unsigned int type) const
    {
      return m_species[type].m_mass;
    }
    ONIKA_HOST_DEVICE_FUNC inline double operator () () const
    {
      return this->operator() (m_type);
    }
  };

  struct ChargeFunctor
  {
    const ParticleSpecie* __restrict__ m_species = nullptr;  
    ONIKA_HOST_DEVICE_FUNC inline double operator () (unsigned int type) const
    {
      return m_species[type].m_mass;
    }
  };

}
