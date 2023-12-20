#pragma once

#include <exanb/fields.h>
#include <exanb/core/basic_types.h>

namespace exaStamp
{
  template<class PartArraysT, bool = PartArraysT::template HasField<field::_virial>::value >
  struct VirialFieldHelper
  {
    static inline void add_contribution(PartArraysT& array, size_t i, const exanb::Mat3d& vir)
    {
      array[field::virial][i] += vir;
    }
  };
  template<class PartArraysT>
  struct VirialFieldHelper<PartArraysT,false>
  {
    static inline void add_contribution(PartArraysT&, size_t, const exanb::Mat3d&) {}
  };

  template<class PartArraysT> 
  static inline void virial_add_contribution(PartArraysT& array, size_t i, const exanb::Mat3d& vir)
  {
    VirialFieldHelper<PartArraysT>::add_contribution(array,i,vir);
  }
}

