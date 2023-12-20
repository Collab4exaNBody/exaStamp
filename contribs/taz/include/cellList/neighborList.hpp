/// @file
/// @brief Definition of the class NeighborList

#ifndef __NEIGHBOR_LIST_HPP_INCLUDED
#define __NEIGHBOR_LIST_HPP_INCLUDED


#include <algorithm>
#include <ostream>
#include <vector>
#include <tuple>

#include "utils/auxMath.hpp"


#include "utils/array/extArray.hpp"
#include "simd/verlet.hpp"


/// @brief Base for the class to handle Neighbors
class NeighborList_base {

protected:

  static uint8_t numberOfTypes; ///< Number of particle types

  /// @brief Default constructor
  NeighborList_base() {}

public:

  /// @brief Shortcut for a type to store a neighbor identification (cell, index and type)
  typedef std::tuple<uint, uint16_t, uint8_t> nbr_id;

  static constexpr uint8_t CELL = 0; ///< Index of the cell in the storage of a neighbor
  static constexpr uint8_t INDX = 1; ///< Index of the index in the storage of a neighbor
  static constexpr uint8_t TYPE = 2; ///< Index of the type in the storage of a neighbor

  /// @brief Setter for the number of types
  static inline void setNumberOfTypes(const uint8_t numTypes) {
    numberOfTypes = numTypes;
  }

  /// @brief Destructor (nothing to do)
  virtual ~NeighborList_base() {}

};





/// @brief Class to handle Neighbors
class NeighborList : public NeighborList_base {

public:

  /// @brief Default constructor
  ///
  ///
  NeighborList() 
    : NeighborList_base(), maxNbr(0), sizePerType(), indexes() {}

  /// @brief Destructor (nothing to do)
  virtual ~NeighborList() {}

  uint maxNumberOfNeighbors() const;

  uint numberOfNeighbors(const uint16_t index) const;
  uint numberOfNeighbors(const uint16_t index, const uint8_t type) const;

  void getNeighbors(const uint16_t index, const uint8_t type, nbr_id*& start, uint& size);
  void getNeighbors(const uint16_t index, nbr_id*& start, uint& size);

  void addNeighbor(const uint16_t index, const nbr_id& nbr);
  void delNeighbor(const uint16_t index, const uint i);

  virtual void clear();

  void check(const uint numberOfParticles);

  template <typename F> void sort(F f);

  uint capacity(const uint16_t index) const;

  void __debug_print(std::ostream& flux);

protected:

  uint maxNbr; ///< Max number of neighbors observed in the lists

  std::vector<std::vector<uint>    > sizePerType;  ///< Number of neighbors of each type per particle
  std::vector<std::vector<nbr_id> > indexes;      ///< List of neighbors per particle

};



/// @brief Accessor to the maximum number of neighbors
inline uint NeighborList::maxNumberOfNeighbors() const {
  return maxNbr;
}


  
/// @brief Accessor to the number of neighbors of a specified particle
/// @param [in] index Particle index
/// @return Number of neighbors
inline uint NeighborList::numberOfNeighbors(const uint16_t index) const {
  return indexes[index].size();
}


 
/// @brief Accessor to the number of neighbors of a specified type from a specified particle
/// @param [in] index Particle index
/// @param [in] typeIndex Neighbors type
/// @return Number of neighbors
inline uint NeighborList::numberOfNeighbors(const uint16_t index, const uint8_t typeIndex) const {
  return sizePerType[index][typeIndex];
}



/// @brief Get the neighbors of a specified type for a particle
/// @param [in] index Particle index
/// @param [in] type Neighbors type
/// @param [out] start First element of the array containing the neighbors
/// @param [out] size Number of neighbors
inline void NeighborList::getNeighbors(const uint16_t index, const uint8_t type, nbr_id*& start, uint& size) {

  uint startIndex = 0;

  auto& spt = sizePerType[index];

  for (uint8_t i=0; i<type; ++i) startIndex += spt[i];

  start = &indexes[index][startIndex];
  size  = spt[type];

}


/// @brief Get the neighbors of all types for a particle
/// @param [in] index Particle index
/// @param [out] start First element of the array containing the neighbors
/// @param [out] size Number of neighbors
inline void NeighborList::getNeighbors(const uint16_t index, nbr_id*& start, uint& size) {

  start = indexes[index].data();
  size  = indexes[index].size();

}


/// @brief Add a neighbor for the specified particule
/// @param [in] index Particle index
/// @param [in] nbr Neighbor data
inline void NeighborList::addNeighbor(const uint16_t index, const nbr_id& nbr) {

  ++sizePerType[index][std::get<TYPE>(nbr)];
  indexes[index].push_back(std::move(nbr));
  
}


/// @brief Delete a neighbor identified
/// @param [in] index Particle index
/// @param [in] i Neighbor position
inline void NeighborList::delNeighbor(const uint16_t index, const uint i) {
	--sizePerType[index][std::get<TYPE>(indexes[index][i])];
	indexes[index][i]=indexes[index][indexes[index].size()-1];
	indexes[index].pop_back();
}


/// @brief Clear all neighbor lists
///
///
inline void NeighborList::clear() {
  
  for (uint index=0, size=indexes.size(); index<size; ++index) {
    sizePerType[index].assign(this->numberOfTypes, 0);
    indexes    [index].clear();
  }

}



/// @brief Reset the containers with a new size
/// @param [in] numberOfParticles New size
inline void NeighborList::check(const uint numberOfParticles) {

  // Store old size
  const uint oldSize = indexes.size();

  // Resize
  sizePerType.resize(numberOfParticles);
  indexes.resize(numberOfParticles);

  // Initialize the new elements of sizePerType
  for (uint i=oldSize; i<numberOfParticles; ++i) {
    sizePerType[i].assign(this->numberOfTypes, 0);
  }

}

/// @brief Sort the neighbors of each particle according to a specified function
/// @tparam F Type of the function
/// @param [in] f Sorting function
template <typename F> void NeighborList::sort(F f) {

  // Reset maxNbr
  maxNbr = 0;

  // For each particle
  for(uint i=0, size=indexes.size(); i<size; ++i) {

    // Update maxNbr
    maxNbr = auxMax<uint>(maxNbr, indexes[i].capacity());

    // Count the number of types
    uint8_t cnt = 0;
    for (uint t=0; t<this->numberOfTypes; ++t)
      if (numberOfNeighbors(i, t) > 0) ++cnt;

    // If there is more than one type, sort indexes according to f
    if (cnt>1) {
      auto& tmp = indexes[i];
      std::sort(tmp.begin(), tmp.end(), f);
    }

  }

}

/// @brief Get the neighbor capacity for a particle
/// @param [in] index Particle index
/// @return Neighbor capacity
inline uint NeighborList::capacity(const uint16_t index) const {
  return indexes[index].capacity();
}



/// @brief Debug print to given flux (not used)
/// @param [in] flux Print flux
inline void NeighborList::__debug_print(std::ostream& flux) {

}


/// @brief Class to handle Neighbors
class NeighborList_Verlet : public NeighborList {

	public:

	/// @brief Default constructor
	///
	///
	NeighborList_Verlet() 
	: store_x(), store_y(), store_z(), store_check() {}

	/// @brief Destructor (nothing to do)
	virtual ~NeighborList_Verlet() {}


	void addCoord(const uint16_t index, double x, double y, double z);
	void addCoord(size_t n, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);
	void getShiftedPointer(size_t shift, double* rx, double* ry, double* rz);
	void clear();
	size_t nbElements() const;

	bool checkVerlet( const double *x_, const double *y_, const double *z_, const double check, uint size);

	void resize_store(int n);

	protected:

	std::vector< double > store_x    ;      ///< store atom position x
	std::vector< double > store_y    ;      ///< store atom position y
	std::vector< double > store_z    ;      ///< store atom position z
	std::vector< double > store_check;      ///< temporary storage
};


/// @brief Clear all neighbor lists
///
///
inline void NeighborList_Verlet::clear() {
  
	for (size_t index=0, size=indexes.size(); index<size; ++index) {
		sizePerType[index].assign(this->numberOfTypes, 0);
		indexes    [index].clear();
	}

	store_x.clear();
	store_y.clear();
	store_z.clear();
	store_check.clear();
}

inline size_t NeighborList_Verlet::nbElements() const {
  return indexes.size();
}

/// @brief Add a neighbor for the specified particule
/// @param [in] index Particle index
/// @param [in] x X component of position of the atom
/// @param [in] y Y component of position of the atom
/// @param [in] z Y component of position of the atom
inline void NeighborList_Verlet::addCoord(const uint16_t index, double x, double y, double z) {
	store_x[index]=x;
	store_y[index]=y;
	store_z[index]=z;
	store_check[index]=0;
}

/// @brief Add a neighbor for the specified particule
/// @param [in] index Particle index
/// @param [in] x X component of position of the atom
/// @param [in] y Y component of position of the atom
/// @param [in] z Y component of position of the atom
inline void NeighborList_Verlet::addCoord(size_t n, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)  {
	std::copy(x.begin(), x.end(), store_x.data() );
	std::copy(y.begin(), y.end(), store_y.data() );
	std::copy(z.begin(), z.end(), store_z.data() );

	for(int i = 0; i < n ; i++)
		store_check[i]=0;
}

inline void NeighborList_Verlet::getShiftedPointer(size_t shift, double* rx, double* ry, double* rz)
{
	rx = store_x.data()+shift;
	ry = store_y.data()+shift;
	rz = store_z.data()+shift;
}

/// @brief Determines whether an atom has moved more than 1/2 of the verlet radius.
/// @param [in] x array of new position x  
/// @param [in] y array of new position y  
/// @param [in] z array of new position z  
/// @param [in] check value of 1/2 of the verlet radius
/// @param [in] size number of atoms in the cell
/// @return Indicates if at least one atom has moved more than 1/2 of the verlet radius.
inline bool NeighborList_Verlet::checkVerlet( const double *x, const double *y, const double *z, const double check, uint size)
{
	return simd::kernels::Verlet<double>(              
		x,              
		y,              
		z,
		store_x.data(), 
		store_y.data(), 
		store_z.data(),
		store_check.data(),           
		size,         
		check
	);
}

/// @brief reallocation of the storage if it is needed
///
/// @param [in] n new number of atom in the cell
inline void NeighborList_Verlet::resize_store(int n){
	store_x.resize(n);
	store_y.resize(n);
	store_z.resize(n);
	store_check.resize(n);
}


#endif // __NEIGHBORS_LIST_HPP_INCLUDED
