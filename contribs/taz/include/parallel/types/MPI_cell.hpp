/// @file 
/// @brief Definition of cell environment MPI types

#ifndef __MPI_CELL_HPP_INCLUDED
#define __MPI_CELL_HPP_INCLUDED


#include "parallel/types/MPI_ghost.hpp"

#include "utils/neighbor.hpp"


/// @brief Type to communicate the owner node of ghost cells
typedef MPI__Ghost<int> MPI__GhostCellOwner;

/// @brief Type to communicate the owner nodes of the neighbors of a cell
struct MPI__CellEnv {

	/// @brief Default constructor
  MPI__CellEnv() : gid(), neighborOwner() {}

  /// @brief Destructor (nothing to do)
  ~MPI__CellEnv() {}

  /// @brief Copy constructor
	/// @param [in] cellEnv Cell environment to copy
  MPI__CellEnv(const MPI__CellEnv& cellEnv) 
    : gid(cellEnv.gid), neighborOwner() {

    for (uint i=0; i<Neighbor::num_neighbors; ++i)
      neighborOwner[i] = cellEnv.neighborOwner[i];
  }

  /// @brief Assignment operator
  /// @param [in] cellEnv Cell environment to copy
  MPI__CellEnv& operator = (const MPI__CellEnv& cellEnv) {

    gid = cellEnv.gid;

    for (uint i=0; i<Neighbor::num_neighbors; ++i)
      neighborOwner[i] = cellEnv.neighborOwner[i];

    return *this;

  }

  uint64_t gid; ///< Global ID of the cell
  uint32_t neighborOwner[Neighbor::num_neighbors] = {0}; ///< Neighboring cell owner node for each direction

};

#endif // __MPI_CELL_HPP_INCLUDED
