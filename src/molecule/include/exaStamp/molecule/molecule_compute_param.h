#pragma once

#include <exaStamp/molecule/molecule_species.h>
#include <exanb/core/basic_types_def.h>
#include <onika/oarray.h>
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>

#include <onika/memory/allocator.h>
#include <onika/cuda/ro_shallow_copy.h>

namespace exaStamp
{
  struct MoleculeGenericFuncParam
  {
    double p0 = 0.0;
    double p1 = 0.0;
    double p2 = 0.0;
    float x0 = 0.0;
    float coeff = 0.0;
    ONIKA_HOST_DEVICE_FUNC inline bool is_null() const { return p0==0.0 && p1==0.0 && p2==0.0 && x0==0.0f && coeff==0.0f; }
    inline bool operator < (const MoleculeGenericFuncParam& m) const
    {
      if( p0 < m.p0 ) return true;
      else if( p0 > m.p0 ) return false;

      if( p1 < m.p1 ) return true;
      else if( p1 > m.p1 ) return false;

      if( p2 < m.p2 ) return true;
      else if( p2 > m.p2 ) return false;

      if( x0 < m.x0 ) return true;
      else if( x0 > m.x0 ) return false;

      if( coeff < m.coeff ) return true;
      else if( coeff > m.coeff ) return false;
      
      return false;
    }
  };

  inline std::ostream& operator << ( std::ostream& out , const MoleculeGenericFuncParam& m )
  {
    return out<<"("<<m.p0<<","<<m.p1<<","<<m.p2<<","<<m.x0<<","<<m.coeff<<")";
  }
  
  struct MoleculeComputeParameterSet
  {
    onika::memory::CudaMMVector<MoleculeGenericFuncParam> m_func_params;
    onika::memory::CudaMMVector<LJExp6RFParms> m_pair_params;
    onika::memory::CudaMMVector<double> m_energy_correction; // per atome type energy correction
    onika::memory::CudaMMVector<Mat3d> m_virial_correction; // per atom type virial correction
    std::map< onika::oarray_t<int,4> , int > m_intramol_param_map;
  };

  struct MoleculeComputeParameterSetRO
  {
    const MoleculeGenericFuncParam * __restrict__ m_func_params = nullptr;
    const LJExp6RFParms * __restrict__ m_pair_params = nullptr;
    
    MoleculeComputeParameterSetRO() = default;
    MoleculeComputeParameterSetRO(const MoleculeComputeParameterSetRO&) = default;
    MoleculeComputeParameterSetRO(MoleculeComputeParameterSetRO&&) = default;
    MoleculeComputeParameterSetRO& operator = (const MoleculeComputeParameterSetRO&) = default;
    MoleculeComputeParameterSetRO& operator = (MoleculeComputeParameterSetRO&&) = default;

    inline MoleculeComputeParameterSetRO( const MoleculeComputeParameterSet& ml )
      : m_func_params( ml.m_func_params.data() )
      , m_pair_params( ml.m_pair_params.data() )
      {}
  };

  struct IntramolecularParameterIndexLists
  {
    onika::memory::CudaMMVector<int> m_bond_param_idx;
    onika::memory::CudaMMVector<int> m_angle_param_idx;
    onika::memory::CudaMMVector<int> m_torsion_param_idx;
    onika::memory::CudaMMVector<int> m_improper_param_idx;
  };

  struct IntramolecularParameterIndexListsRO
  {
    const int * const __restrict__ m_bond_param_idx = nullptr;
    const int * const __restrict__ m_angle_param_idx = nullptr;
    const int * const __restrict__ m_torsion_param_idx = nullptr;
    const int * const __restrict__ m_improper_param_idx = nullptr;
    
    IntramolecularParameterIndexListsRO() = default;
    IntramolecularParameterIndexListsRO(const IntramolecularParameterIndexListsRO&) = default;
    IntramolecularParameterIndexListsRO(IntramolecularParameterIndexListsRO&&) = default;
    IntramolecularParameterIndexListsRO& operator = (const IntramolecularParameterIndexListsRO&) = default;
    IntramolecularParameterIndexListsRO& operator = (IntramolecularParameterIndexListsRO&&) = default;

    inline IntramolecularParameterIndexListsRO( const IntramolecularParameterIndexLists& ipil )
      : m_bond_param_idx( ipil.m_bond_param_idx.data() )
      , m_angle_param_idx( ipil.m_angle_param_idx.data() )
      , m_torsion_param_idx( ipil.m_torsion_param_idx.data() )
      , m_improper_param_idx( ipil.m_improper_param_idx.data() )
      {}
  };

}

