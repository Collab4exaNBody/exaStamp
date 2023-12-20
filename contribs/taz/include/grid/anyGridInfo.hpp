/// @file 
/// @brief Definition of info class for any grid

#ifndef __ANY_GRID_INFO_HPP_INCLUDED
#define __ANY_GRID_INFO_HPP_INCLUDED


#include <ostream>
#include <vector>

#include "globals.hpp"

#include "parallel/types/MPI_particle.hpp"

#include "utils/curve.hpp"
#include "utils/neighbor.hpp"
#include "utils/array/array.hpp"

#include "domain/domainInfo.hpp"


class AnyDecomposition;
template <class T> class Moduler;


/// @brief Class gathering all the data to build and use any grid
class AnyGridInfo {

public:

  AnyGridInfo(uint index);

  /// @brief Destructor
	///
	///
  virtual ~AnyGridInfo() {
    m_decomposition=nullptr;
  }

  /// @brief Accessor to m_numberOfCells
  int getNumberOfCells() const { return m_numberOfCells; }
  /// @brief Get number of cells (real + ghost)
  /// @return Total number of cells
  uint getTotalNumberOfCells() const { return m_numberOfCells + m_numberOfGhostCells; }

  /// @brief Accessor to m_ghostThickness
  int getGhostThickness() const { return m_ghostThickness; }

  /// @brief Accessor to m_cellInf
  const vec3<int>& getCellInf() { return m_cellInf; }
  /// @brief Accessor to m_cellSup
  const vec3<int>& getCellSup() { return m_cellSup; }

  /// @brief Accessor to m_limInf
  const vec3<double>& getLimitInf() { return m_limInf; }
  /// @brief Accessor to m_limSup
  const vec3<double>& getLimitSup() { return m_limSup; }

  int index(const vec3<int>& coord);

  int getNeighborDomain(const vec3<int>& mov, const vec3<int>& cell);

  /// @brief Accessor to m_decomposition
  AnyDecomposition* getDecomposition() { return m_decomposition; }

  bool inGrid(const vec3<int>& r);
  int inGridAndGetIndex(const vec3<int>& r);
  bool inGrid(const vec3<double>& r);
  void inGrid(const Particle* particles, bool* result, const uint size);

  template <class Q, typename F>
  void findCellIndex(const Q* q, uint* indexes, const uint size, F f);

  virtual void print(std::ostream& flux);

  vec3<bool> getGlobalEdges();

  /// @brief Accessor to m_neighbors
  const Array<int>& getNeighborArray() { return m_neighbors; }

  vec3<int> coords(uint index);
  

	/// @brief Get the index of an octree
	/// @param [in] coord Coordinates of the octree
	int indexForAMR(const vec3<int>& coord) 
	{ 
		 bool tmp;
		uint idx;
		find(coord-1, tmp, idx); ///> The cartesian positions of octree are shifted of 1 

		return (tmp && isReal(idx)? idx : -1 );
	}

	/// @brief Get the index of an octree
	/// @param [in] coord Coordinates of the octree
	int indexForAMRWithGhost(const vec3<int>& coord) { 

		bool tmp;
		uint idx;
		find(coord-1, tmp, idx);///> The cartesian positions of octree are shifted of 1 

		return (tmp ? idx : -1 );
	}


  bool isReal(const uint index);

  void setNeighbors(const Array<uint>& cells, const std::vector< Array<int> >& neighbors, const Moduler<int>& modI);

  /// @brief Get neighbor domains for specified cell (not used)
  /// @param [in] i Cell index
  const std::vector<int>& getNeighbors(const uint i) { return m_myNeighborDomains[i]; }

  /// @brief Get the global index of a cell from its coordinates
  /// @param [in] coords Cell coordinates
  uint getGlobalId(const vec3<int>& coords) { return m_globalCurveNoGhost.convert(coords); }
  /// @brief Get the coordinates of a cell from its global index (not used)
  /// @param [in] globId Cell index
  vec3<int> getGlobalCoords(const vec3<int>& globId) { return m_globalCurveNoGhost.convert(globId); }

  void find(const uint globalIndex, bool& found, uint& idx);

public:

  void find(const vec3<int>& coord, bool& found, uint& idx);

  uint m_domainIndex; ///< Domain index

  LCurve<int> m_globalCurve; ///< Tool to map the whole system + a ghost layer with a 1D index
  LCurve<int> m_globalCurveNoGhost; ///< Tool to map the whole system with a 1D index

  int m_ghostThickness; ///< Ghost thickness

  vec3<int> m_cellInf; ///< Lower limit of the grid (in term of cells)
  vec3<int> m_cellSup; ///< Upper limit of the grid (in term of cells)

  vec3<double> m_limInf; ///< Lower limit of the grid (in term of distances)
  vec3<double> m_limSup; ///< Upper limit of the grid (in term of distances)

  uint m_numberOfCells; ///< Number of real cells
  uint m_numberOfGhostCells; ///< Number of ghost cells

  AnyDecomposition* m_decomposition; ///< Link to the system decomposition

  Array<int> m_neighbors; ///< Grid neighbors

  std::vector<uint> m_myGlobalIndexes; ///< Global indexes of the cells of the grid
  std::vector< std::vector<int> > m_myNeighborDomains; ///< Neighbors of the cells of the grid

};


/// @brief Get local index of a cell
/// @param [in] coord Cell coordiantes
/// @return Index
inline int AnyGridInfo::index(const vec3<int>& coord) { 
  bool tmp;
  uint idx;
  find(coord, tmp, idx);

  return ( tmp ? idx : -1 );
}


/// @brief Get the neighbor domain of a cell in a direction
/// @param [in] mov Direction of the neighbor
/// @param [in] cell Cell coordinates
/// @return Neighbor index or null if no neighbor
inline int AnyGridInfo::getNeighborDomain(const vec3<int>& mov, const vec3<int>& cell) {
  
  uint idx = index(cell);

  if (m_myNeighborDomains[idx].size()>0) 
    return m_myNeighborDomains[idx][Neighbor::getIndex(mov)];
  else
    return Neighbor::null;

  return Neighbor::null;

}

/// @brief Check if a cell is in the grid
/// @param [in] r Cell coordinates (in term of cells)
/// @return True if cell is in the gri
inline int AnyGridInfo::inGridAndGetIndex(const vec3<int>& r) { 
  bool found;
  uint idx;
  find(r, found, idx);
  if(found && isReal(idx)) return idx;
  else return -1;
}


/// @brief Check if a cell is in the grid
/// @param [in] r Cell coordinates (in term of cells)
/// @return True if cell is in the gri
inline bool AnyGridInfo::inGrid(const vec3<int>& r) { 
  bool found;
  uint idx;
  find(r, found, idx);
  return found && isReal(idx);
}


/// @brief Check if a particle is in the grid
/// @param [in] r Particle coordinates (in term of distance)
/// @return True if particle is in the grid
inline bool AnyGridInfo::inGrid(const vec3<double>& r) {
  return inGrid( auxFloor<int>((r-Global::domainInfo.getMinBounds()) / Global::domainInfo.getCellLength()) );
}


/// @brief Check if some particles are in the grid (not used)
/// @param [in] particles Particles
/// @param [out] result Tab to store the result (true for each particle in the grid, false for the others)
/// @param [in] size Number of partiles to process
inline void AnyGridInfo::inGrid(const MPI__Particle* particles, bool* result, const uint size) {
  for (uint i=0; i<size; ++i) { 
    auto& p = particles[i];
    result[i] = inGrid( std::move( vec3<double>(p.r.x, p.r.y, p.r.z) ) );
  }
}


/// @brief Get the index of some particles cells
/// @tparam Q Particles type
/// @tparam F Function type
/// @param [in] q Particles
/// @param [out] indexes Returned indexes
/// @param [in] size Number of particles
/// @param [in] f Function to identify the cell
template <class Q, typename F>
inline void AnyGridInfo::findCellIndex(const Q* q, uint* indexes, const uint size, F f) { 
  for (uint i=0; i<size; ++i)
    indexes[i] = (uint) index( f(q[i]) );
}

/// @brief Get coordinates of a cell from its local index
/// @param [in] index Cell index
/// @return Cell coordinates
inline vec3<int> AnyGridInfo::coords(uint index  /* 0 to m_myGlobalIndexes.size() */) {
  return m_globalCurve.convert(m_myGlobalIndexes[index]) - m_ghostThickness;
}


/// @brief Check if cell is real and in the domain
/// @param [in] index Cell index
/// @return True if real
inline bool AnyGridInfo::isReal(const uint index /* 0 to m_myGlobalIndexes.size() */) {

	// Check if the cell is in the grid
  if (index>=m_myGlobalIndexes.size()) 
    return false;
  
  // Coordinates of the cell in the whole system with ghost
  auto globalCoords0 = m_globalCurve.convert(m_myGlobalIndexes[index]); 
  
  // Coordinates of the cell in the whole system without ghost
  auto globalCoords1 = globalCoords0 - m_ghostThickness;
  
  // Check if the cell is global ghost
  if (countNumberOf(true, globalCoords1<0 || globalCoords1>=Global::domainInfo.getNumberOfCellsPerDim()) > 0)
    return false;
  
  // Check in decomposition if cell is real
  return m_decomposition->isMine(m_globalCurveNoGhost.convert(globalCoords1));
  
}

#endif // __ANY_GRID_INFO_HPP_INCLUDED
