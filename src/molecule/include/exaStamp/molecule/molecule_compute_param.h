#pragma once

#include <exaStamp/molecule/molecule_species.h>
#include <onika/math/basic_types_def.h>
#include <exanb/core/mt_concurrent_map.h>

#include <onika/oarray.h>
#include <exaStamp/potential/ljexp6rf/ljexp6rf.h>

#include <onika/memory/allocator.h>
#include <onika/cuda/ro_shallow_copy.h>

namespace exaStamp
{
  static inline uint64_t bond_key(int t0, int t1)
  {
    assert( t0 >= 0 && t0 < (1<<15) );
    assert( t1 >= 0 && t1 < (1<<15) );
    if( t0 > t1 ) std::swap( t0 , t1 );
    static constexpr int t2 = (1<<16) - 1;
    static constexpr int t3 = (1<<16) - 1;
    return ( uint64_t(t0)<<48 ) | ( uint64_t(t1)<<32 ) | ( uint64_t(t2)<<16 ) | uint64_t(t3);
  }

  static inline uint64_t angle_key(int t0, int t1, int t2)
  {
    assert( t0 >= 0 && t0 < (1<<15) );
    assert( t1 >= 0 && t1 < (1<<15) );
    assert( t2 >= 0 && t2 < (1<<15) );
    if( t0 > t2 ) std::swap( t0 , t2 );
    static constexpr int t3 = (1<<16) - 1;
    return ( uint64_t(t0)<<48 ) | ( uint64_t(t1)<<32 ) | ( uint64_t(t2)<<16 ) | uint64_t(t3);
  }

  static inline uint64_t torsion_key(int t0, int t1, int t2, int t3)
  {
    assert( t0 >= 0 && t0 < (1<<15) );
    assert( t1 >= 0 && t1 < (1<<15) );
    assert( t2 >= 0 && t2 < (1<<15) );
    assert( t3 >= 0 && t3 < (1<<15) );
    if( t0 > t3 ) { std::swap( t0 , t3 ); std::swap( t1 , t2 ); }
    else if (t0 == t3 && t1 > t2 ) {std::swap( t0 , t3 ); std::swap( t1 , t2 ); }
    return ( uint64_t(t0)<<48 ) | ( uint64_t(t1)<<32 ) | ( uint64_t(t2)<<16 ) | uint64_t(t3);
  }

  static inline uint64_t improper_key(int t0, int t1, int t2, int t3)
  {
    assert( t0 >= 0 && t0 < (1<<15) );
    assert( t1 >= 0 && t1 < (1<<15) );
    assert( t2 >= 0 && t2 < (1<<15) );
    assert( t3 >= 0 && t3 < (1<<15) );
    uint64_t key     = (uint64_t(t0)<<48) | (uint64_t(t1)<<32) | (uint64_t(t2)<<16) | uint64_t(t3);
    uint64_t alt_key = (uint64_t(t0)<<48) | (uint64_t(t2)<<32) | (uint64_t(t3)<<16) | uint64_t(t1);
    if( alt_key < key ) key = alt_key;
    alt_key          = (uint64_t(t0)<<48) | (uint64_t(t3)<<32) | (uint64_t(t1)<<16) | uint64_t(t2);
    if( alt_key < key ) key = alt_key;
    alt_key          = (uint64_t(t0)<<48) | (uint64_t(t2)<<32) | (uint64_t(t1)<<16) | uint64_t(t3);
    if( alt_key < key ) key = alt_key;
    alt_key          = (uint64_t(t0)<<48) | (uint64_t(t1)<<32) | (uint64_t(t3)<<16) | uint64_t(t2);
    if( alt_key < key ) key = alt_key;
    alt_key          = (uint64_t(t0)<<48) | (uint64_t(t3)<<32) | (uint64_t(t2)<<16) | uint64_t(t1);
    if( alt_key < key ) key = alt_key;
    return ((1ull<<48)-1) ^ key; // so that improper keys are different from torsion keys even if types are the same
  }

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

    inline bool operator == (const MoleculeGenericFuncParam& m) const
    {
      return p0==m.p0 && p1==m.p1 && p2==m.p2 && x0==m.x0 && coeff==m.coeff;
    }
  };

  inline std::ostream& operator << ( std::ostream& out , const MoleculeGenericFuncParam& m )
  {
    return out<<"("<<m.p0<<","<<m.p1<<","<<m.p2<<","<<m.x0<<","<<m.coeff<<")";
  }
  
  using ChemicalPairPotMap = exanb::MultiThreadedConcurrentMap< std::unordered_map< onika::oarray_t<uint64_t,2>, int > >;
  
  struct IntramolecularPairParams
  {
    LJExp6RFParms m_param;
    float m_pair_weight = 1.0f;
    float m_rf_weight = 1.0f;
    inline bool operator < (const IntramolecularPairParams& ipp) const
    {
           if( m_param < ipp.m_param ) return true;
      else if( ipp.m_param < m_param ) return false;
      else if( m_pair_weight < ipp.m_pair_weight ) return true;
      else if( m_pair_weight > ipp.m_pair_weight ) return false;
      else if( m_rf_weight < ipp.m_rf_weight ) return true;
      else if( m_rf_weight > ipp.m_rf_weight ) return false;
      else return false;
    }
    inline bool operator == (const IntramolecularPairParams& ipp) const
    {
      return m_param == ipp.m_param && m_pair_weight == ipp.m_pair_weight && m_rf_weight == ipp.m_rf_weight;
    }
  };
  
  struct MoleculeComputeParameterSet
  {
    onika::memory::CudaMMVector<MoleculeGenericFuncParam> m_func_params;
    onika::memory::CudaMMVector<IntramolecularPairParams> m_pair_params;
    onika::memory::CudaMMVector<double> m_energy_correction; // per atome type energy correction
    onika::memory::CudaMMVector<Mat3d> m_virial_correction; // per atom type virial correction
    std::map< uint64_t , int > m_intramol_param_map;
  };

  struct MoleculeComputeParameterSetRO
  {
    const MoleculeGenericFuncParam * __restrict__ m_func_params = nullptr;
    const IntramolecularPairParams * __restrict__ m_pair_params = nullptr;
    
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
    onika::memory::CudaMMVector<int> m_pair_param_idx;
  };

  struct IntramolecularParameterIndexListsRO
  {
    const int * const __restrict__ m_bond_param_idx = nullptr;
    const int * const __restrict__ m_angle_param_idx = nullptr;
    const int * const __restrict__ m_torsion_param_idx = nullptr;
    const int * const __restrict__ m_improper_param_idx = nullptr;
    const int * const __restrict__ m_pair_param_idx = nullptr;
    
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
      , m_pair_param_idx( ipil.m_pair_param_idx.data() )
      {}
  };

}

