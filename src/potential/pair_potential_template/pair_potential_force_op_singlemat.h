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

#include "pair_potential_template.h"
#include <exanb/compute/compute_pair_traits.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <cmath>
#include <cstdlib>

#ifndef PRIV_NAMESPACE_NAME
#define PRIV_NAMESPACE_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_details)
#endif

// default definitions doing nothing
#ifndef DEBUG_ADDITIONAL_FIELDS
#define DEBUG_ADDITIONAL_FIELDS /**/
#endif

#ifndef DEBUG_ADDITIONAL_PARAMETERS
#define DEBUG_ADDITIONAL_PARAMETERS /**/
#endif

namespace exaStamp
{
  using namespace exanb;

  namespace PRIV_NAMESPACE_NAME
  {
    // helper functor, with templated call operator (with/without weights)
    struct ForceOp
    {
      onika::cuda::ro_shallow_copy_t< USTAMP_POTENTIAL_PARAMS > p;
      PairPotentialMinimalParameters pair_params;
      double ecut;

      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (Vec3d dr,double d2,double& fx,double& fy,double& fz,Mat3d& vir,DEBUG_ADDITIONAL_PARAMETERS CellParticlesT,size_t,size_t, double weight ) const
      {
        // double ep = 0.0;
        const double r = sqrt(d2);
        double e=0.0, de=0.0;
        
#       if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
        USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de , weight );
#       else
        USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de );
        e *= weight;
        de *= weight;
#       endif
        e -= ecut * weight;
        de /= r;
             
        const Vec3d fe = { de * dr.x , de * dr.y , de * dr.z };
        fx += fe.x;
        fy += fe.y;
        fz += fe.z;
        // ep += .5 * e;
        vir += tensor(fe,dr) * -0.5;
      }

      // ComputeBuffer less computation without virial
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (Vec3d dr,double d2,double& fx,double& fy,double& fz,DEBUG_ADDITIONAL_PARAMETERS CellParticlesT,size_t,size_t, double weight ) const
      {
        // double ep = 0.0;
        const double r = sqrt(d2);
        double e=0.0, de=0.0;

#       if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
        USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de , weight );
#       else
        USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de );
        e *= weight;
        de *= weight;
#       endif
        e -= ecut * weight;
        de /= r;

        fx += de * dr.x;
        fy += de * dr.y;
        fz += de * dr.z;
        // ep += .5 * e;
      }

      // ComputeBuffer less computation with virial
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (Vec3d dr,double d2,double& ep,double& fx,double& fy,double& fz,Mat3d& vir,DEBUG_ADDITIONAL_PARAMETERS CellParticlesT,size_t,size_t, double weight ) const
      {
        const double r = sqrt(d2);
        double e=0.0, de=0.0;

#       if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
        USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de , weight );
#       else
        USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de );
        e *= weight;
        de *= weight;
#       endif
        e -= ecut * weight;
        de /= r;
        
        const Vec3d fe = { de * dr.x , de * dr.y , de * dr.z };
        fx += fe.x;
        fy += fe.y;
        fz += fe.z;
        ep += .5 * e;
        vir += tensor(fe,dr) * -0.5;
      }

      // ComputeBuffer less computation without virial
      template<class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (Vec3d dr,double d2,double& ep,double& fx,double& fy,double& fz,DEBUG_ADDITIONAL_PARAMETERS CellParticlesT,size_t,size_t, double weight ) const
      {
        const double r = sqrt(d2);
        double e=0.0, de=0.0;
        
#       if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
        USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de , weight );
#       else
        USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de );
        e *= weight;
        de *= weight;
#       endif
        e -= ecut * weight;
        de /= r;
            
        fx += de * dr.x;
        fy += de * dr.y;
        fz += de * dr.z;
        ep += .5 * e;
      }

      // without virial computation
      template<bool UseWeights, bool UseNeighbors, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBuffer2<UseWeights,UseNeighbors>& tab,
        double& ep,
        double& ax,
        double& ay,
        double& az,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT
        ) const
      {
        using Mat3d = ::exanb::FakeMat3d;
        FakeMat3d virial; // unused virial
#       include "force_op_impl2.hxx"
      }

      // without virial computation
      template<bool UseWeights, bool UseNeighbors, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBuffer2<UseWeights,UseNeighbors>& tab,
        double& ep,
        double& ax,
        double& ay,
        double& az,
        Mat3d& virial,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT
        ) const
      {
#       include "force_op_impl2.hxx"
      }

      // without virial computation
      template<bool UseWeights, bool UseNeighbors, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBuffer2<UseWeights,UseNeighbors>& tab,
        double& ax,
        double& ay,
        double& az,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT
        ) const
      {
        [[maybe_unused]] double ep = 0.0;
        using Mat3d = ::exanb::FakeMat3d;
        FakeMat3d virial; // unused virial
#       include "force_op_impl2.hxx"
      }

      // without virial computation
      template<bool UseWeights, bool UseNeighbors, class CellParticlesT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (
        size_t n,
        ComputePairBuffer2<UseWeights,UseNeighbors>& tab,
        double& ax,
        double& ay,
        double& az,
        Mat3d& virial,
        DEBUG_ADDITIONAL_PARAMETERS
        CellParticlesT
        ) const
      {
        [[maybe_unused]] double ep = 0.0;
#       include "force_op_impl2.hxx"
      }

    };

  }

}

namespace exanb
{
	
  template<> struct ComputePairTraits<exaStamp::PRIV_NAMESPACE_NAME::ForceOp>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible = true;
    static inline constexpr bool BufferLessCompatible = true;
#   ifdef USTAMP_POTENTIAL_ENABLE_CUDA
    static inline constexpr bool CudaCompatible = true;
#   else
    static inline constexpr bool CudaCompatible = false;
#   endif
  };

}

