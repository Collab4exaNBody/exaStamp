/// @file 
/// @brief Tools to handle domain decomposition (class Decomposition)

#ifndef __DECOMPOSITION_HPP_INCLUDED
#define __DECOMPOSITION_HPP_INCLUDED


#include <vector>
#include <ostream>
//
#include "utils/array/array.hpp"
#include "utils/vec3/vec3.hpp"


template <class T> class Configuration;
class DomainInterface;


/// @brief Base class to handle domain decomposition on domains with a structure in cells
class Decomposition {

public:

  /// @brief Constructor
  /// @param [in] n Number of domains, default=0
  Decomposition(int n=0)
    : numberOfDomains(n) {}

  /// @brief Destructor (nothing to do)
  virtual ~Decomposition() {}

  /// @brief Accessor to number of domains
  uint getNumberOfDomains() { return numberOfDomains; }

  /// @brief Print the type of decomposition and number of domains in specified flux
  /// @param [in,out] flux Print flux
  virtual void print(std::ostream& flux) = 0;

protected:

  uint numberOfDomains; ///< Number of domains in the decomposition

};


/// @brief Decomposition into rectangular domains
class RectilinearDecomposition : public Decomposition {

public:

  RectilinearDecomposition(Configuration<DomainInterface>& configuration);

  virtual ~RectilinearDecomposition();

  /// @brief Get origin of specified domain
  /// @param [in] index Index of the domain
  /// @return Coordinates of domain origin
  const vec3<int>& origins(const int index) {
    return m_origins[index];
  }

  /// @brief Get size of specified domain
  /// @param [in] index Index of the domain
  /// @return Size of the domain in each direction
  const vec3<int>& sizes(const int index) {
    return m_sizes[index];
  }

  /// @brief Get one neighbor of specified domain
  /// @param [in] domainIndex Index of the domain
  /// @param [in] nbrIndex Number of the neighbor to get
  /// @return Neighbor index
  const int& neighbors(const int domainIndex, const int nbrIndex) {
    return m_neighbors[domainIndex][nbrIndex];
  }

  virtual void print(std::ostream& flux);
  /// @brief Print debug in specified flux (empty and not used)
  /// @param [in,out] flux Print flux
  void debug_print(std::ostream& flux){}

private:

  vec3<int> m_numberOfDomainsPerDim; ///<  Number of domain in each dimension

  Array< vec3<int>  > m_origins; ///< Origin of each domain
  Array< vec3<int>  > m_sizes; ///< Sizes of each domain
  Array< Array<int> > m_neighbors; ///< Neighbors list for each domain

};


/// @brief Decompostion into arbitrary domains (defined by user or load balancer)
class AnyDecomposition : public Decomposition {

public:

  AnyDecomposition(Configuration<DomainInterface>& configuration, int rank);

  /// @brief Destructor (nothing to do)
  virtual ~AnyDecomposition() {}

  /// @brief Get cells of the domain
  /// @param [out] ptr Pointer to an array containing the cells
  /// @param [out] size Number of cells
  void myGlobalIndexes(uint*& ptr, uint& size) { ptr=m_myGlobalIndexes.data(); size=m_myGlobalIndexes.size(); }

  bool isMine(const uint index);
  int neighborIndex(const uint index, const uint nbrDir);
  /// @brief Reset neighbors
  ///
  ///
  void releaseEnv() { m_neighborOwners.clear(); }

  void update_part_1(const std::vector<uint>& importIndexes, const std::vector<uint>& exportIndexes, const std::vector<uint32_t>& importEnv);
  void update_part_2(uint index, const Array<int>& neighbors, const Array<int>& allCellOwners);

  virtual void print(std::ostream& flux);

//private:
public:
  //
  int m_rank; ///< Domain rank
  std::vector<uint>                m_myGlobalIndexes; ///< Indexes of the cells in the domain, ghost not included
  // 1D global indexes -- ghost NOT included -- map on curve WITHOUT ghosts -- MUST BE SORTED !
  std::vector< std::vector<uint> > m_neighborOwners; ///< Ranks of domains where neighbors can be found for each cell
  // same size as above, if 2nd vector is empty, this is an inside cell. Else, array of size Neighbor::num_neighbors with owner rank of neighboring cells

};


Decomposition* Decompose(Configuration<DomainInterface>& configuration, int numberOfNodes, int rank);

#endif // __DECOMPOSITION_HPP_INCLUDED
