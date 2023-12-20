/*
 * MPI_ghost.hpp
 *
 *  Created on: Oct 14, 2016
 *      Author: giarda
 */
/// @file
/// @brief Definition of types to exchange ghost versions of the particles

#ifndef MPI_GHOST_HPP_
#define MPI_GHOST_HPP_


#include "utils/vec3/vec3.hpp"


/// @brief Type to communicate objects of the ghost
/// @tparam P Class of the objects to communicate
template <class P> class MPI__Ghost {

public:

  /// @brief Shortcut for the class of the objects to communicate (not used)
  typedef P base_t;

  /// @brief Default constructor
  MPI__Ghost() {}

  /// @brief Destructor (nothing to do)
  ~MPI__Ghost() {}

  /// @brief Constructor from a cell and an object
  /// @tparam P Class of the objects to communicate
  /// @param [in] cell_ Coordinates of the recipient cell
  /// @param [in] base_ Object to communicate
  MPI__Ghost(const vec3<int>& cell_, const P& base_)
    : cell(cell_), base(base_) {}

  /// @brief Copy constructor
  /// @tparam P Class of the objects to communicate
	/// @param [in] p MPI_Ghost<p> to copy
  MPI__Ghost(const MPI__Ghost<P>& p)
    : cell(p.cell), base(p.base) {}

  /// @brief Assignment operator
  /// @tparam P Class of the objects to communicate
	/// @param [in] p MPI_Ghost<p> to copy
  MPI__Ghost<P>& operator = (const MPI__Ghost<P>& p) {
    cell = p.cell;
    base = p.base;
    return *this;
  }

  vec3<int> cell; ///< Coordinates of the recipient cell

  P base; ///< Object to communicate

};

#endif /* MPI_GHOST_HPP_ */
