
#pragma once

#include <onika/memory/allocator.h>
#include <onika/cuda/cuda.h>
#include <vector>
#include <utility>
#include <exanb/compute/compute_pair_buffer.h>

#include <exanb/fields.h>
#include <exanb/core/config.h>
#include <exanb/core/declare_field.h>

XSTAMP_DECLARE_FIELD(double , dEmb ,"Embedding derivative");

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  struct PhiRhoCutoff
  {
    double m_phi_cutoff;
    double m_rho_cutoff;
  };

  struct EamParticleEmbField
  {
    onika::memory::CudaMMVector< double > m_emb;
  };

  template<class EamMultimatPotentialParmT>
  struct EamMultimatPotentialScratch
  {
    onika::memory::CudaMMVector< PhiRhoCutoff > m_phi_rho_cutoff;
    onika::memory::CudaMMVector< uint8_t > m_pair_enabled;
    onika::memory::CudaMMVector< EamMultimatPotentialParmT > m_ro_potentials;
  };

  struct EAMSpecyPairInfo
  {
    int32_t m_type_a = 0;
    int32_t m_type_b = 0;
    double m_charge_a = 0.0;
    double m_charge_b = 0.0;
  };

  template<class EamParametersT>
  struct EamMultimatParameters
  {
    EamParametersT m_parameters;
    EAMSpecyPairInfo m_specy_pair;
    std::string m_type_a;
    std::string m_type_b;
  };

  template<class EamParametersT>
  struct EamMultimatParametersRO
  {
    EamParametersT m_parameters;
    EAMSpecyPairInfo m_specy_pair;
  };

}

