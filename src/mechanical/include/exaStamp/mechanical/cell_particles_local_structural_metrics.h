#pragma once

#include <vector>
#include <exanb/core/grid_particle_field_accessor.h>

namespace exaStamp
{

  struct CellParticleLocalStructuralMetrics
  {
    // The size of vectors below equals the number of particles in cell
    // All variables are per-atom quantities

    // Structural metrics used as descriptors in ML models (Bispectrum for now)
    std::vector< std::vector<double> > bispectrum; // Bispectrum
    std::vector< std::vector<double> > steinhart;  // Steinhart

    // Crystal structure assigned using logistic regression (ML model)
    std::vector< int > crystal_structure;
    
  };

  using GridParticleLocalStructuralMetrics = std::vector< CellParticleLocalStructuralMetrics >;

  /*
   * Allows a particle field, stored as a flat 1D array, to be used as a particle field through the ExternalCellParticleFieldAccessor envelop
   */
  template<class FieldIdT> struct GridLocalStructuralMetricsAccessor
  {
    using value_type = typename FieldIdT::value_type;
    using reference_t = const value_type &;
    const CellParticleLocalStructuralMetrics * __restrict__ mech_data = nullptr;
    template<class CellsT>
    ONIKA_HOST_DEVICE_FUNC
    inline reference_t operator () ( size_t cell_i, size_t p_i , CellsT ) const
    {
      return get_from_mechanical_data( mech_data, cell_i, p_i, FieldIdT{} );
    }
  };

}

#define STRUCTURAL_DECLARE_FIELD(name,member) \
XNB_DECLARE_FIELD( std::remove_reference_t<decltype(exaStamp::CellParticleLocalStructuralMetrics{}.member[0])> , name , #name ); \
namespace exaStamp { \
  static inline const auto& get_from_mechanical_data( const CellParticleLocalStructuralMetrics * __restrict__ mech_data, size_t cell_i, size_t p_i, onika::soatl::FieldId<field::_##name> ) \
  { \
    return mech_data[cell_i].member[p_i]; \
  } \
  static inline auto structural_field( const CellParticleLocalStructuralMetrics * __restrict__ mech_data, onika::soatl::FieldId<field::_##name> ) \
  { \
    using FieldIdT = onika::soatl::FieldId<field::_##name>; \
    using ExternalFieldT = exanb::ExternalCellParticleFieldAccessor< GridLocalStructuralMetricsAccessor<FieldIdT> , FieldIdT >; \
    return ExternalFieldT{ .m_func = { mech_data } }; \
  } \
}

STRUCTURAL_DECLARE_FIELD(bispectrum,bispectrum)
STRUCTURAL_DECLARE_FIELD(steinhart,steinhart)
STRUCTURAL_DECLARE_FIELD(crystal,crystal_structure)


