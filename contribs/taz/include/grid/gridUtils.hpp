/// @file
/// @brief Definition of the traversal manager

#ifndef __GRID_UTILS_HPP_INCLUDED
#define __GRID_UTILS_HPP_INCLUDED


#include <map>
#include <vector>

#include "utils/array/array.hpp"
#include "utils/array/extArray.hpp"
#include "utils/vec3/vec3.hpp"


class AnyGridInfo;
class RectilinearGridInfo;


/// @brief Enumeration of base traversals defining place where force computation is called
///
/// Associated with level to get traversals with cells
enum Traversal {
	ALL, ///< All cells
	INSIDE, ///< Inside cells
	EDGE ///< Edges cells
};


/// @brief Tool to manage traversals
class TraversalManager {

public:

	/// @brief Enumeration of traversals
  enum CellTraversal {
    EMPTY,         ///< No cell
    ALL,           ///< All cells (used for debug),
    REAL,          ///< Real cells
    GHOST,         ///< Ghost cells
    ALL_BUT_ONE,   ///< All cells minus one layer of ghost
    INSIDE_ONE,    ///< Real cells minus one layer
    EDGE_ONE,      ///< One layer of real cells on the edge
    INSIDE_GTHICK, ///< Real cells minus ghost thick layer
    EDGE_GTHICK,   ///< Ghost thick layer of real cells on grid edge
    EDGE_MAKE_NBR, ///< Ghost minus last layer plus ghost thick layer of real cells
    EDGE_GTHICK_PLUS_ONE  ///< Ghost thick layer of real cells plus one layer of ghost
  };

  TraversalManager();
  ~TraversalManager();

  void makeTraversals(AnyGridInfo*         info, const std::vector< vec3<int> >& coords);
  void makeTraversals(RectilinearGridInfo* info, const std::vector< vec3<int> >& coords);
  // void makeTraversals(VoronoiGridInfo*     info, const std::vector< vec3<int> >& coords) {}

  Array<uint>&  getTraversal(CellTraversal traversal);
  CellTraversal getCellTraversal(::Traversal T, uint level=0);

  /// @brief Print the number of cells in some traversals (not used)
  void print() {

    std::cout<< "EMPTY      " << " size : " << getTraversal(EMPTY      ).size() << std::endl;
    std::cout<< "ALL        " << " size : " << getTraversal(ALL        ).size() << std::endl;
    std::cout<< "REAL       " << " size : " << getTraversal(REAL       ).size() << std::endl;
    std::cout<< "GHOST      " << " size : " << getTraversal(GHOST      ).size() << std::endl;
    std::cout<< "ALL_BUT_ONE" << " size : " << getTraversal(ALL_BUT_ONE).size() << std::endl;
    std::cout<< "INSIDE_ONE " << " size : " << getTraversal(INSIDE_ONE ).size() << std::endl;
    std::cout<< "EDGE_ONE   " << " size : " << getTraversal(EDGE_ONE   ).size() << std::endl;

  }

protected:

  void addTraversal(CellTraversal traversal, ExtArray<uint>* cells);
  void addTraversal(CellTraversal traversal, std::vector<uint>& cells);
  void addTraversal(CellTraversal t1, CellTraversal t2);

  std::map< CellTraversal, Array<uint>* > traversals; ///< Map linking traversals to their cells
  std::map< CellTraversal, bool >        destroyMap; ///< Map indexing the traversals that must be destroyed (not a copy)

};

#endif // __GRID_UTILS_HPP_INCLUDED
