#pragma once

#include <cstdint>
#include <vector>
#include <cstdlib>

#include <exanb/core/basic_types.h>

namespace exaStamp
{

  using namespace exanb;

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

  template<class FieldIdT> struct IsMechanicalExternalField : public std::false_type {};

  namespace mechanical
  {
    template<class FieldIdT> struct ExternalField;
  }

# define MECHANICAL_DECLARE_FIELD(field,member) \
  namespace mechanical { \
  struct _##field {}; \
  template<> struct ExternalField< _##field > { \
    using value_type = std::remove_reference_t< decltype( CellParticleLocalMechanicalMetrics{}. member [0] ) >; \
    static const char* short_name() { return #field; } \
    static const char* name() { return #field; } \
    static inline value_type get_from_mechanical_data(const CellParticleLocalMechanicalMetrics * __restrict__ mech_data, size_t cell_i, size_t p_j) { return mech_data[cell_i]. member [p_j] ; } \
  }; \
  static constexpr ExternalField< _##field > field = {}; \
  } \
  template<> struct IsMechanicalExternalField< mechanical::ExternalField< mechanical::_##field > > : public std::true_type {};

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

  template<class FieldIdT> static inline constexpr bool is_mechanical_field_v = IsMechanicalExternalField<FieldIdT>::value ;

  template<class CellsT>
  struct GridParticleMechanicalAccessor
  {
    CellsT m_cells;
    const CellParticleLocalMechanicalMetrics * __restrict__ mech_data = nullptr;
    template<class FieldIdT>
    ONIKA_HOST_DEVICE_FUNC inline typename FieldIdT::value_type get(size_t cell_i, size_t p_j, const FieldIdT& f ) const
    {
      if constexpr ( is_mechanical_field_v<FieldIdT> )
      {
        return FieldIdT::get_from_mechanical_data( mech_data , cell_i, p_j );
      }
      if constexpr ( ! is_mechanical_field_v<FieldIdT> )
      {
        return m_cells[cell_i][f][p_j];
      }
      return {};
    }
  };

}

