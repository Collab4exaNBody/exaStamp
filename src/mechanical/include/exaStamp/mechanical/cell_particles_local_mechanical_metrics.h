#pragma once

#include <cstdint>
#include <vector>
#include <cstdlib>

#include <onika/math/basic_types.h>
#include <exanb/core/grid_particle_field_accessor.h>

namespace exaStamp
{

  struct CellParticleLocalMechanicalMetrics
  {
    // The size of vectors below equals the number of particles in cell
    // All variables are per-atom quantities
    
    // Continuum mechanics local measures (static)
    std::vector< Mat3d > F;           // Deformation gradient tensor
    std::vector< Mat3d > E;           // Green-Lagrange strain tensor
    std::vector< Mat3d > R;           // Pure rotation tensor
    std::vector< Mat3d > U;           // Pure stretch tensor
    std::vector< Vec3d > mu;          // Microrotation vector
    std::vector< Vec3d > s;           // Slip vector
    std::vector< Vec3d > l;           // Vector parallel to Burgers vector in the glide plane
    std::vector< Vec3d > m;           // Vector perpendicular to Burgers vector in the glide plane
    std::vector< Vec3d > n;           // Glide plane vector

    // Continuum mechanics local measures (dynamic)
    std::vector< Mat3d > L;           // Velocity gradient tensor
    std::vector< Vec3d > phi;         // Vorticity vector

    // Continuum mechanics local measures (dislocation analysis)
    std::vector< Mat3d> vector_gradient_tensor; // Spatial gradient of a vector quantity (in the present form it is the spatial gradient of the microrotation vector, to identify dislocations)
    std::vector< double > dislo_indic;          // Velocity gradient tensor
    std::vector< double > vis;                  // How close to a screw dislocation am I
    std::vector< double > coin;                 // How close to an edge dislocation am I
    std::vector< Vec3d > dislo_line;            // Vector parallel to dislocation line in dislocation glide plane
    std::vector< Vec3d > dislo_line_ortho;      // Vector perpendicular to dislocation line in dislocation glide plane
    
  };

  using GridParticleLocalMechanicalMetrics = std::vector< CellParticleLocalMechanicalMetrics >;

  /*
   * Allows a particle field, stored as a flat 1D array, to be used as a particle field through the ExternalCellParticleFieldAccessor envelop
   */
  template<class FieldIdT> struct GridLocalMechanicalMetricsAccessor
  {
    using value_type = typename FieldIdT::value_type;
    using reference_t = const value_type &;
    const CellParticleLocalMechanicalMetrics * __restrict__ mech_data = nullptr;
    template<class CellsT>
    ONIKA_HOST_DEVICE_FUNC
    inline reference_t operator () ( size_t cell_i, size_t p_i , CellsT ) const
    {
      return get_from_mechanical_data( mech_data, cell_i, p_i, FieldIdT{} );
    }
  };

}

#define MECHANICAL_DECLARE_FIELD(name,member) \
XSTAMP_DECLARE_FIELD( std::remove_reference_t<decltype(exaStamp::CellParticleLocalMechanicalMetrics{}.member[0])> , name , #name ); \
namespace exaStamp { \
  static inline const auto& get_from_mechanical_data( const CellParticleLocalMechanicalMetrics * __restrict__ mech_data, size_t cell_i, size_t p_i, onika::soatl::FieldId<field::_##name> ) \
  { \
    return mech_data[cell_i].member[p_i]; \
  } \
  static inline auto mechanical_field( const CellParticleLocalMechanicalMetrics * __restrict__ mech_data, onika::soatl::FieldId<field::_##name> ) \
  { \
    using FieldIdT = onika::soatl::FieldId<field::_##name>; \
    using ExternalFieldT = exanb::ExternalCellParticleFieldAccessor< GridLocalMechanicalMetricsAccessor<FieldIdT> , FieldIdT >; \
    return ExternalFieldT{ .m_func = { mech_data } }; \
  } \
}

MECHANICAL_DECLARE_FIELD(defgrad,F)
MECHANICAL_DECLARE_FIELD(greenlag,E)
MECHANICAL_DECLARE_FIELD(rot,R)
MECHANICAL_DECLARE_FIELD(stretch,U)
MECHANICAL_DECLARE_FIELD(microrot,mu)
MECHANICAL_DECLARE_FIELD(slip,s)
MECHANICAL_DECLARE_FIELD(burgerpar,l)
MECHANICAL_DECLARE_FIELD(burgerortho,m)
MECHANICAL_DECLARE_FIELD(glide,n)
MECHANICAL_DECLARE_FIELD(velgrad,L)
MECHANICAL_DECLARE_FIELD(vort,phi)
MECHANICAL_DECLARE_FIELD(vecgrad,vector_gradient_tensor)
MECHANICAL_DECLARE_FIELD(dislo,dislo_indic)
MECHANICAL_DECLARE_FIELD(vis,vis)
MECHANICAL_DECLARE_FIELD(coin,coin)
MECHANICAL_DECLARE_FIELD(dislol,dislo_line)
MECHANICAL_DECLARE_FIELD(dislolo,dislo_line_ortho)

