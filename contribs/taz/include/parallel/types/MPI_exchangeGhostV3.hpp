/// @file 
/// @brief

#ifndef __MPI_EXCHANGE_GHOST_V3_HPP_INCLUDED
#define __MPI_EXCHANGE_GHOST_V3_HPP_INCLUDED


#include "utils/vec3/vec3.hpp"


///@brief Mpi type to exchange vec3
struct ExchangeGhostV3 {
  /// @brief Default constructor
  ExchangeGhostV3() {}
  /// @brief Destructor
  ~ExchangeGhostV3() {}

  /// @brief Constructor
  /// @param [in] id_ Particle ID
  /// @param [in] v3_ 3d data
  /// @param [in] cell_ Particle cell
  ExchangeGhostV3(uint id_, vec3<double> v3_, const vec3<int>& cell_) 
    : id(id_), v3(v3_), cell(cell_) {}

  uint id; ///< Particle ID
  vec3<double> v3; ///< Vec3 data
  vec3<int> cell; ///< Particle cell

};

#endif // __MPI_EXCHANGE_GHOST_VEL_HPP_INCLUDED
