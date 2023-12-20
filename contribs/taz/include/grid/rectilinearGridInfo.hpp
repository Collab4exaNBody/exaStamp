/// @file 
/// @brief Definition of info class for rectilinear grid

#ifndef __RECTILINEAR_GRID_INFO_HPP_INCLUDED
#define __RECTILINEAR_GRID_INFO_HPP_INCLUDED


#include <ostream>

#include "globals.hpp"

#include "parallel/types/MPI_particle.hpp"

#include "utils/curve.hpp"
#include "utils/array/array.hpp"
	

#include<cassert>

class RectilinearDecomposition;



/// @brief Class gathering all the data to build and use a rectilinear grid
class RectilinearGridInfo {

public:

  RectilinearGridInfo(uint index);

  /// @brief Destructor (nothing to do)
  virtual ~RectilinearGridInfo() {}

  const vec3<int>& getOrigin() const;
  const vec3<int>& getNumberOfCellsPerDim() const;

  int getNumberOfCells() const;
  int getTotalNumberOfCells() const;

  int getGhostThickness() const;

  const vec3<double>& getMinBounds() const;
  const vec3<double>& getMaxBounds() const;
  const vec3<double>& getExtension() const;
  
  vec3<double> getLimitInf();
  vec3<double> getLimitSup();

  vec3<int> getCellInf();
  vec3<int> getCellSup();

	/// @brief Get the index of a cell
	/// @param [in] coord Coordinates of the cell
	int index(const vec3<int>& coord) 
	{ 
		assert(convert(m_shift+coord) >=0 );
		return convert(m_shift+coord); 
	}
  
	/// @brief Get the index of an octree
	/// @param [in] coord Coordinates of the octree
	int indexForAMR(const vec3<int>& coord) 
	{ 
		assert(index(coord-1) >= 0);
		assert(index(coord-1) < (int) ( 	
			(numberOfCellsPerDim.x+2)*
			(numberOfCellsPerDim.y+2)*
			(numberOfCellsPerDim.z+2))
		);
		return index(coord-1);  ///> The cartesian positions of octree are shifted of 1 
	}

	/// @brief Get the index of an octree 
	/// @param [in] coord Coordinates of the octree
	int indexForAMRWithGhost(const vec3<int>& coord) 
	{ 
		auto inGridsup = coord -1 + m_shift < numberOfCellsPerDim+2;
		auto inGridinf = coord -1 + m_shift  >= 0;
		if(inGridsup.x && inGridsup.y && inGridsup.z && inGridinf.x && inGridinf.y && inGridinf.z)
			return index(coord-1); ///> The cartesian positions of octree are shifted of 1 
		else return -1;
	}

  inline vec3<int> coords(const uint index) {return convert(index)- m_shift;}

  int getNeighborDomain(const vec3<int>& mov, const vec3<int>& cell);

  RectilinearDecomposition* getDecomposition();

  bool inGrid(const vec3<int>& r);
  int inGridAndGetIndex(const vec3<int>& r);
  bool inGrid(const vec3<double>& r);
  void inGrid(const Particle* particles, bool* result, const uint size);

  template <class Q, typename F>
  void findCellIndex(const Q* q, uint* indexes, const uint size, F f);

  virtual void print(std::ostream& flux);

  vec3<bool> getGlobalEdges();

  const Array<int>& getNeighborArray();

  // ONLY for rectilinear, above is everywhere !

  vec3<int> convert(const int& index);
  int convert(const vec3<int>& coord);

protected:

  vec3<int> origin; ///< Origin of the grid (in term of cells)
  vec3<int> numberOfCellsPerDim;	///< Number of cells per dim

  LCurve<int> curve; ///< Tool to map the grid with a 1D index

  int ghostThickness; ///< Ghost thickness

  vec3<double> minBounds; ///< Lower boundary of the grid in each dimension (in term of distance)
  vec3<double> maxBounds; ///< Upper boundary of the grid in each dimension (in term of distance)
  vec3<double> extension; ///< Extension of the grid in each dimension (in term of distance)

  vec3<double> invCellLength; ///< Inverse of the cell length in each dimension
  vec3<double> globalMin; ///< Lower boundary of the global system in each dimension

  uint numberOfCells; ///< Number of cells in the grid

  RectilinearDecomposition* decomposition; ///< Link to the system decomposition

  uint numberOfNeighbors; ///< Number of neighbor domains
  Array<int> neighbors; ///< Neighbor domains

  vec3<int> m_shift; ///< Shift to apply to cells local coordinates to obtain cells global coordinates

  vec3<int> freeOrNot;


};



/// @brief Accessor to origin
inline const vec3<int>& RectilinearGridInfo::getOrigin() const {
  return origin;
}


/// @brief Accessor to numberOfCellsPerDim
inline const vec3<int>& RectilinearGridInfo::getNumberOfCellsPerDim() const { 
  return numberOfCellsPerDim;
}


/// @brief Accessor to numberOfCells
inline int RectilinearGridInfo::getNumberOfCells() const { 
  return numberOfCells;
}


/// @brief Get number of cells (real + ghost)
/// @return Total number of cells
inline int RectilinearGridInfo::getTotalNumberOfCells() const { 
  return product(numberOfCellsPerDim + 2*ghostThickness);
}


/// @brief Accessor to ghostThickness
inline int RectilinearGridInfo::getGhostThickness() const {
  return ghostThickness;
}


/// @brief Accessor to minBounds
inline const vec3<double>& RectilinearGridInfo::getMinBounds() const {
  return minBounds;
}


/// @brief Accessor to maxBounds
inline const vec3<double>& RectilinearGridInfo::getMaxBounds() const {
 return maxBounds;
}


/// @brief Accessor to extension (not used)
inline const vec3<double>& RectilinearGridInfo::getExtension() const {
 return extension;
}


/// @brief Get the lower limits of the grid (in term of distance)
/// @return Lower limit in each dimension
inline vec3<double> RectilinearGridInfo::getLimitInf() {
  return getMinBounds();
}


/// @brief Get the upper limits of the grid (in term of distance)
/// @return Upper limit in each dimension
inline vec3<double> RectilinearGridInfo::getLimitSup() {
  return getMaxBounds();
}


/// @brief Get the lower limits of the grid (in term of cells)
/// @return Lower limit in each dimension
inline vec3<int> RectilinearGridInfo::getCellInf() {
  return origin;
}


/// @brief Get the upper limits of the grid (in term of cells)
/// @return Upper limit in each dimension
inline vec3<int> RectilinearGridInfo::getCellSup() {
  return origin+numberOfCellsPerDim;
}


/// @brief Accessor to decomposition
inline RectilinearDecomposition* RectilinearGridInfo::getDecomposition() {
  return decomposition;
}


/// @brief Convert index to coordinates in the grid
/// @param [in] index 1D index
/// @return Coordinates
inline vec3<int> RectilinearGridInfo::convert(const int& index) {
  return curve.convert(index);
}


/// @brief Convert coordinates in the grid to index
/// @param [in] coord Coordinates
/// @return 1D index
inline int RectilinearGridInfo::convert(const vec3<int>& coord) {
  return curve.convert(coord);
}


/// @brief Check if a cell is in the grid
/// @param [in] r Cell coordinates (in term of cells)
/// @return True if cell is in the grid
inline int RectilinearGridInfo::inGridAndGetIndex(const vec3<int>& r) {

  vec3<bool> min = r >= origin;
  vec3<bool> max = r <  origin + numberOfCellsPerDim;

  bool boolInGrid = min.x && min.y && min.z && max.x && max.y && max.z;
  
  if(boolInGrid) return convert(r+m_shift);
  else return -1;
}

/// @brief Check if a cell is in the grid
/// @param [in] r Cell coordinates (in term of cells)
/// @return True if cell is in the grid
inline bool RectilinearGridInfo::inGrid(const vec3<int>& r) {

  vec3<bool> min = r >= origin;
  vec3<bool> max = r <  origin + numberOfCellsPerDim;

  return min.x && min.y && min.z && max.x && max.y && max.z;
}

/// @brief Check if a particle is in the grid
/// @param [in] r Particle coordinates (in term of distance)
/// @return True if particle is in the grid
inline bool RectilinearGridInfo::inGrid(const vec3<double>& r) {

	return inGrid( auxFloor<int>( invCellLength*(r-globalMin) ) );
}


/// @brief Check if some particles are in the grid (not used)
/// @param [in] particles Particles
/// @param [out] result Tab to store the result (true for each particle in the grid, false for the others)
/// @param [in] size Number of partiles to process
inline void RectilinearGridInfo::inGrid(const MPI__Particle* particles, bool* result, const uint size) {

  for (uint i=0; i<size; ++i) {

    vec3<double> r(particles[i].r.x, particles[i].r.y, particles[i].r.z);

     result[i] = inGrid( auxFloor<int>( invCellLength*(r-globalMin) ) );

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
inline void RectilinearGridInfo::findCellIndex(const Q* q, uint* indexes, const uint size, F f) {

  for (uint i=0; i<size; ++i) {
    indexes[i] = (uint) index( f(q[i]) );
  }

}


/// @brief Accessor to neighbors
inline const Array<int>& RectilinearGridInfo::getNeighborArray() {
  return neighbors;
}

#endif // __RECTILINEAR_GRID_INFO_HPP_INCLUDED
