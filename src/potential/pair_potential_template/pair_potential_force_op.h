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
#include <exanb/compute/compute_pair_function.h>
#include <exaStamp/potential_factory/pair_potential.h>

#undef _CLASS_NAME
#undef _CLASS_BASE

#ifdef USTAMP_POTENTIAL_WITH_VIRIAL
#define _CLASS_NAME USTAMP_POTENTIAL_FORCE_VIRIAL_OP_NAME
#define _CLASS_BASE PairPotentialComputeVirialOperator
#else
#define _CLASS_NAME USTAMP_POTENTIAL_FORCE_OP_NAME
#define _CLASS_BASE PairPotentialComputeOperator
#endif

namespace exaStamp
{
  using namespace exanb;

  class _CLASS_NAME final : public _CLASS_BASE
  {
  public:

    inline _CLASS_NAME(const USTAMP_POTENTIAL_PARAMS& p) : m_potential_params(p) {}
    ~_CLASS_NAME() override final = default;

    inline void set_parameters( const PairPotentialParameters& params) override final
    {
      m_pair_params = params;
    }

    inline const std::string& name() const override final
    {
      static const std::string pot_name(USTAMP_POTENTIAL_STRING);
      return pot_name;
    }

    uint64_t signature() const override final;
    void set_rcut( double rcut ) override final ;
    
    void operator() (
      ComputePairBuffer2<false,false>& tab,
      double& ep,
      double& ax, double& ay,double& az
#     ifdef USTAMP_POTENTIAL_WITH_VIRIAL
      , Mat3d& virial
#     endif
      ) const noexcept override final ;

    void operator() (
      ComputePairBuffer2<true,false>& tab,
      double& ep,
      double& ax, double& ay,double& az
#     ifdef USTAMP_POTENTIAL_WITH_VIRIAL
      , Mat3d& virial
#     endif
      ) const noexcept override final ;

  private:
    const USTAMP_POTENTIAL_PARAMS m_potential_params;
    PairPotentialMinimalParameters m_pair_params;
    double m_ecut = 0.;
    double m_rcut = 0.;
  };

}

