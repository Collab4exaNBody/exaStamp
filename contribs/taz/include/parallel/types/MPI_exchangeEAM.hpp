/// @file 
/// @brief Definition of EAM related exchange types

#ifndef __MPI_EXCHANGE_EAM_HPP_INCLUDED
#define __MPI_EXCHANGE_EAM_HPP_INCLUDED


#include "utils/vec3/vec3.hpp"


/// @brief Type to exchange EAM data
struct ExchangeEAM {

  /// @brief Default constructor
  ExchangeEAM() {}
  /// @brief Destructor (nothing to do)
  ~ExchangeEAM() {}

  /// @brief Constructor
  /// @param [in] id_ Particle ID
  /// @param [in] emb_ EAM data
  /// @param [in] cell_ Particle cell
  ExchangeEAM(uint id_, double emb_, const vec3<int>& cell_) 
    : id(id_), emb(emb_), cell(cell_) {}

  uint id; ///< Particle ID
  double emb; ///< EAM data
  vec3<int> cell; ///< Particle cell

};

#endif // __MPI_EXCHANGE_EAM_HPP_INCLUDED
