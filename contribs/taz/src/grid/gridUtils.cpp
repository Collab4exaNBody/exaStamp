/// @file
/// @brief Implementation of the traversal manager


#include <algorithm>
#include <vector>

#include "globals.hpp"
#include "referenceMap.hpp"

#include "grid/gridUtils.hpp"
#include "grid/anyGridInfo.hpp"
#include "grid/rectilinearGridInfo.hpp"
#include "grid/zone.hpp"


/// @brief Default constructor
///
///
TraversalManager::TraversalManager() 
  : traversals(), destroyMap() {
}


/// @brief Destructor
///
/// Delete arrays of cells
TraversalManager::~TraversalManager() {

  for (auto& elem : traversals) {

    if (destroyMap[elem.first]) {
      delete elem.second;
    }

    elem.second = nullptr;

  }
  
  traversals.clear();
  destroyMap.clear();

}


/// @brief Get the cells in specified traversal
/// @param [in] traversal Traversal
/// @return Cells of the traversal
Array<uint>& TraversalManager::getTraversal(CellTraversal traversal) {
  auto it = traversals.find(traversal);
  if (it!=traversals.end()) return *(it->second);
  else                      return *traversals[EMPTY];
}


/// @brief Convert a base traversal into a usable traversal (CellTraversal)
/// @param [in] T Base traversal
/// @param [in] level Level depending on ghost thickness used for the conversion, default=0
/// @return Usable traversal (with cells)
TraversalManager::CellTraversal TraversalManager::getCellTraversal(::Traversal T, uint level) {

  // NOTE : if ghostThickness > 3, this will have to change ...
  const uint& ghostThickness = Global::reference.getGhostThickness();
  level = auxMin(level, ghostThickness-1);

  CellTraversal t;

  switch (T) {

  case Traversal::ALL:
    if      (level==0)                t=CellTraversal::REAL;
    else if (level==ghostThickness-1) t=CellTraversal::ALL_BUT_ONE;
    //    else if (level==1)                t=CellTraversal::ALL_BUT_TWO;
    else                              t=CellTraversal::EMPTY;
    break;

  case Traversal::EDGE: 
    if      (level==0)                t=CellTraversal::EDGE_GTHICK;
    else if (level==ghostThickness-1) t=CellTraversal::EDGE_MAKE_NBR;
    else if (level==1)                t=CellTraversal::EDGE_GTHICK_PLUS_ONE;
    else                              t=CellTraversal::EMPTY;
    break;

  case Traversal::INSIDE:
    t=CellTraversal::INSIDE_GTHICK;
    break;

  default:
    t=CellTraversal::EMPTY;
    break;
  }

  return t;

}


/// @brief Add a traversal and his cells (stored in an ExtArray) to the map
/// @param [in] traversal Traversal to add
/// @param [in] cells Cells of the traversal
void TraversalManager::addTraversal(CellTraversal traversal, ExtArray<uint>* cells) {

  auto it=traversals.find(traversal);

  if (it!=traversals.end())
    return;

  Array<uint>* array = new Array<uint>(cells->size(), 0);
  for (uint i=0; i<cells->size(); ++i) (*array)[i] = (*cells)[i];

  traversals[traversal] = array;
  destroyMap[traversal] = true;
  
}


/// @brief Add a traversal and his cells (stored in a vector) to the map
/// @param [in] traversal Traversal to add
/// @param [in] cells Cells of the traversal
void TraversalManager::addTraversal(CellTraversal traversal, std::vector<uint>& cells) {

  auto it=traversals.find(traversal);

  if (it!=traversals.end())
    return;

  std::sort(cells.begin(), cells.end());
  auto jt = std::unique(cells.begin(), cells.end());
  cells.resize(std::distance(cells.begin(), jt));

  Array<uint>* array = new Array<uint>(cells.size(), 0);
  for (uint i=0; i<cells.size(); ++i) (*array)[i] = cells[i];

  traversals[traversal] = array;
  destroyMap[traversal] = true;
  
}


/// @brief Copy an existing traversal in a new traversal
/// @param [in] t1 New traversal
/// @param [in] t2 Copied traversal
void TraversalManager::addTraversal(CellTraversal t1, CellTraversal t2) {

  auto it=traversals.find(t1);

  if (it!=traversals.end()) 
    return;

  traversals[t1] = traversals[t2];
  destroyMap[t1] = false;
}


/// @brief Initialize all the traversals from the grid info (case AnyGrid)
/// @param [in] info Grid info
/// @param [in] coords Coordinates of the cells in the grid
/// @warning Defined for ghost thickness of 1 only
void TraversalManager::makeTraversals(AnyGridInfo* info, const std::vector< vec3<int> >& coords) {

  const int ghostThickness            = info->getGhostThickness();
  const bool test                     = (ghostThickness == 1);

  if (ghostThickness>1) exit(0);

  std::vector<uint> indexes(0);

  addTraversal(CellTraversal::EMPTY, indexes);
  indexes.clear();

  indexes.reserve(coords.size());
  for (uint i=0; i<coords.size(); ++i) 
    indexes.push_back(i);
  addTraversal(CellTraversal::ALL,   indexes);
  indexes.clear();

  std::vector<uint> rl(0);
  std::vector<uint> gt(0);
  for (uint i=0; i<coords.size(); ++i) {
    if (info->isReal(i))
      rl.push_back(i);
    else
      gt.push_back(i);
  }
  addTraversal(CellTraversal::REAL,  rl);
  addTraversal(CellTraversal::GHOST, gt);
  // rl.clear(); // do not clear, needed further
  gt.clear();

  if (test) {
    addTraversal(CellTraversal::ALL_BUT_ONE, CellTraversal::REAL);
  }
  else {
    // addTraversal(CellTraversal::ALL_BUT_ONE, indexes);
  }

  std::vector<uint> isd(0);
  std::vector<uint> edg(0);

  for (uint i=0; i<rl.size(); ++i) {

    bool test2 = true;
    for (uint8_t nbr=0; nbr<Neighbor::num_neighbors; ++nbr) {
      if ( !info->isReal( info->index(coords[rl[i]]+Neighbor::getCoords(nbr) ) ) /* if neighbor is not real */) {
	test2=false;
	break;
      }
    }

    if (test2) 
      isd.push_back(rl[i]);
    else
      edg.push_back(rl[i]);

  }
  addTraversal(CellTraversal::INSIDE_ONE, isd);
  addTraversal(CellTraversal::EDGE_ONE,   edg);
  isd.clear();
  edg.clear();

  if (test) {
    addTraversal(CellTraversal::INSIDE_GTHICK, CellTraversal::INSIDE_ONE);
  }
  else {
    // addTraversal(CellTraversal::INSIDE_GTHICK, indexes);
  }

  if (test) {
    addTraversal(CellTraversal::EDGE_GTHICK, CellTraversal::EDGE_ONE);
  }
  else {
    // addTraversal(CellTraversal::EDGE_GTHICK, indexes);
  }

  if (test) {
    addTraversal(CellTraversal::EDGE_MAKE_NBR, CellTraversal::EDGE_ONE);
  }
  else {
    // addTraversal(CellTraversal::EDGE_MAKE_NBR, indexes); 
  }

  if (test) {
    addTraversal(CellTraversal::EDGE_GTHICK_PLUS_ONE, CellTraversal::EMPTY);
  }
  else {
    // addTraversal(CellTraversal::EDGE_GTHICK_PLUS_ONE, indexes);
  }

}


/// @brief Initialize all the traversals from the grid info (case RectilinearGrid)
/// @param [in] info Grid info
/// @param [in] coords Coordinates of the cells in the grid
void TraversalManager::makeTraversals(RectilinearGridInfo* info, const std::vector< vec3<int> >& coords) {

  const vec3<int> origin              = info->getOrigin();
  const vec3<int> numberOfCellsPerDim = info->getNumberOfCellsPerDim();
  const int ghostThickness            = info->getGhostThickness();
  const bool test                     = (ghostThickness == 1);

  ExtArray<uint> indexes(0);

  auto createTraversal4 = [&] (const vec3<int>& a, const vec3<int>& b, const vec3<int>& c, const vec3<int>& d)
    -> void {

    indexes.clear();

    Zone zone(origin+a,origin+b,origin+c,origin+d);
    indexes.reserve(zone.estimatedSize());
    
    for (uint i=0; i<coords.size(); ++i) {
      if (zone.isIn(coords[i])) 
	indexes.push_back(i);
    }
    
  };
  
  auto createTraversal2 = [&] (const vec3<int>& a, const vec3<int>& b)
    -> void {

    indexes.clear();

    ZoneTwo zone(origin+a,origin+b);
    indexes.reserve(zone.estimatedSize());

    for (uint i=0; i<coords.size(); ++i) {
      if (zone.isIn(coords[i])) 
	indexes.push_back(i);
    }
    
  };

  addTraversal(CellTraversal::EMPTY, &indexes);

  createTraversal2(-ghostThickness, numberOfCellsPerDim+ghostThickness);
  addTraversal(CellTraversal::ALL, &indexes);

  createTraversal2(zeros(), numberOfCellsPerDim);
  addTraversal(CellTraversal::REAL, &indexes);

  createTraversal4(-ghostThickness*ones(), zeros(), numberOfCellsPerDim, numberOfCellsPerDim+ghostThickness);
  addTraversal(CellTraversal::GHOST, &indexes);

  if (test) {
    addTraversal(CellTraversal::ALL_BUT_ONE, CellTraversal::REAL);
  }
  else {
    createTraversal2(ones()-ghostThickness, numberOfCellsPerDim+ghostThickness-1);
    addTraversal(CellTraversal::ALL_BUT_ONE, &indexes);
  }

  createTraversal2(ones(), numberOfCellsPerDim-1);
  addTraversal(CellTraversal::INSIDE_ONE, &indexes);

  createTraversal4(zeros(), ones(), numberOfCellsPerDim-1, numberOfCellsPerDim);
  addTraversal(CellTraversal::EDGE_ONE, &indexes);

  if (test) {
    addTraversal(CellTraversal::INSIDE_GTHICK, CellTraversal::INSIDE_ONE);
  }
  else {
    createTraversal2(ghostThickness*ones(), numberOfCellsPerDim-ghostThickness);
    addTraversal(CellTraversal::INSIDE_GTHICK, &indexes);
  }

  if (test) {
    addTraversal(CellTraversal::EDGE_GTHICK, CellTraversal::EDGE_ONE);
  }
  else {
    createTraversal4(zeros(), ghostThickness*ones(), numberOfCellsPerDim-ghostThickness, numberOfCellsPerDim);
    addTraversal(CellTraversal::EDGE_GTHICK, &indexes);
  }

  if (test) {
    addTraversal(CellTraversal::EDGE_MAKE_NBR, CellTraversal::EDGE_ONE);
  }
  else {
    createTraversal4(ones()-ghostThickness, ghostThickness*ones(), numberOfCellsPerDim-ghostThickness, numberOfCellsPerDim+ghostThickness-1);
    addTraversal(CellTraversal::EDGE_MAKE_NBR, &indexes); 
  }

  if (test) {
    addTraversal(CellTraversal::EDGE_GTHICK_PLUS_ONE, CellTraversal::EMPTY);
  }
  else {
    createTraversal4(zeros()-1, ghostThickness*ones(), numberOfCellsPerDim-ghostThickness, numberOfCellsPerDim+1);
    addTraversal(CellTraversal::EDGE_GTHICK_PLUS_ONE, &indexes);
  }

}
