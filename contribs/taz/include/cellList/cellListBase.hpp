/// @file
/// @brief Definition of the class cellListBase

#ifndef __CELL_LIST_BASE_HPP_INCLUDED
#define __CELL_LIST_BASE_HPP_INCLUDED


#include "cellList/EAMStorage.hpp"
#include "cellList/neighborList.hpp"
#include "cellList/vectorizationBuffer.hpp"

#include "parallel/thread/thread.hpp"


/// @brief Fix a chunk for the ExtArray of cellList
#define CELL_LIST_CHUNK 32
#define align simd::vector_t<double>::align

/// @brief Base class to handle a cell : everything except the particles
/// @tparam align Alignment
class CellListBase {

  /// @brief Shortcut for the type of the cellListBase
  typedef CellListBase self_type;

public:

  static VectBuffer<align>& getVectBuffer();  
  static VectBuffer<align>& getVectBufferWithScreenterm(); // Meam 
  static VectBuffer<align>& getVectBuffer(const uint& n);   
  static VectBuffer<align>& getVectBufferWithScreenterm(const uint& n);  // Meam

  /// @brief Default constructor
  ///
  ///
  CellListBase()
    : mtx(), 
      marks(0),      
      neighborList(), 
      isEAM(Global::reference.isEAM()),
      m_eamStorage(),
      m_numPerType(Global::reference.getNumberOfTypes(), 0) {
  }

  /// @brief Destructor (nothing to do)
  virtual ~CellListBase() {}

  void lock  (uint index);
  void unlock(uint index);

  uint getMaxNumberOfNeighbors();
  void clearNeighborList();

  void mark(const uint& index);
  void clear_marked();

  void resetEAMData();

  void embCopy(self_type& to) const;

  double& rho(const uint i);
  double& emb(const uint i);

  const double& rho(const uint i) const;
  const double& emb(const uint i) const;


  /// @brief Accessor to the number of particles of a type
  /// @param [in] type Type
  /// @return Number of particules
  uint getNumberOfParticlesPerType(const uint8_t type) { return m_numPerType[type]; }

//protected:
  public:

  static Th_< VectBuffer<align> > m_shared_vectBuffer; ///< Local memory-aligned arrays to perform vectorized operations

  void checkLocals(const uint& numberOfParticles);

  MMutex mtx; ///< Array of mutexes

  ExtArray<uint, CELL_LIST_CHUNK, 0> marks;; ///< List of marked particles

  NeighborList_Verlet neighborList; ///< Store neighbouring particle with raduis cut-off = rcut potential + rverlet
  //NeighborList neighborList; ///< List of neighbors for each particle


  bool isEAM;              ///< Indicates if EAM potential is used and therefore if EAMStorage must be handled
  EAMStorage m_eamStorage; ///< Store the data specific to EAM & MEAM potential computation

  Array<uint> m_numPerType; ///< Number of particle per type

};




/// @brief Initialization of m_shared_vectBuffer
 
Th_< VectBuffer<align> > CellListBase::m_shared_vectBuffer = Th_< VectBuffer<align> >();


/// @brief Get the local version of the vectorization buffer
/// @return Local vectorization buffer
  
inline VectBuffer<align>& CellListBase::getVectBuffer() {
  return m_shared_vectBuffer.local();
}


/// @brief Get the local version of the vectorization buffer
/// @return Local vectorization buffer
  
inline VectBuffer<align>& CellListBase::getVectBufferWithScreenterm() {
  return m_shared_vectBuffer.local();
}


/// @brief Get a local version of the vectorization buffer of a specified size
/// @param [in] n Size
/// @return Local vectorization buffer
  
inline VectBuffer<align>& CellListBase::getVectBuffer(const uint& n) {
  // Get the local version
  auto& vb = m_shared_vectBuffer.local();
  // Resize to n
  vb.resize(n);
  return vb;
}


/// @brief Get a local version of the vectorization buffer of a specified size
/// @param [in] n Size
/// @return Local vectorization buffer

inline VectBuffer<align>& CellListBase::getVectBufferWithScreenterm(const uint& n) {
  // Get the local version
  auto& vb = m_shared_vectBuffer.local();
  // Resize to n
  vb.resizeMeam(n);
  return vb;
}


/// @brief Resize local storages
/// @param [in] numberOfParticles New size
  
inline void CellListBase::checkLocals(const uint& numberOfParticles) {

  this->mtx.check(numberOfParticles);

  if (isEAM) 
    m_eamStorage.check(numberOfParticles);
}


/// @brief Lock a mutex
/// @param [in] index Index of the mutex
  
inline void CellListBase::lock(uint index) {
  mtx.lock(index);
}


/// @brief Unlock a mutex
/// @param [in] index Index of the mutex
  
inline void CellListBase::unlock(uint index) {
  mtx.unlock(index);
}


/// @brief Get the maximum number of neighbor in the neighbor lists
/// @return Max number of neighbors
  
inline uint CellListBase::getMaxNumberOfNeighbors() {
  return neighborList.maxNumberOfNeighbors();
}


/// @brief Clear the neighbor lists
///
///
  
inline void CellListBase::clearNeighborList() {
  neighborList.clear();
}


/// @brief Mark a particle
/// @param [in] index Index of the particle to mark
  
inline void CellListBase::mark(const uint& index) {
  marks.push_back(index);
}


/// @brief Unmark all particles
///
///
  
inline void CellListBase::clear_marked() {
  marks.clear();
}


/// @brief Reset the EAM storage
///
///
  
inline void CellListBase::resetEAMData() {
  m_eamStorage.reset();
}


/// @brief Copy the embedding terms to another cell
/// @param [in] to Recipient cell
  
inline void CellListBase::embCopy(self_type& to) const {
  m_eamStorage.embCopy(to.m_eamStorage);
}


/// @brief Get a reference to rho term of a specified particle
/// @param [in] i Particle index
/// @return Reference to rho term
  
inline double& CellListBase::rho(const uint i) {
  return m_eamStorage.rho(i);
}


/// @brief Get a reference to embedding term of a specified particle
/// @param [in] i Particle index
/// @return Reference to embedding term
  
inline double& CellListBase::emb(const uint i) {
  return m_eamStorage.emb(i);
}


/// @brief Get a constant reference to rho term of a specified particle
/// @param [in] i Particle index
/// @return Constant reference to rho term
  
inline const double& CellListBase::rho(const uint i) const {
  return m_eamStorage.rho(i);
}


/// @brief Get a constant reference to embedding term of a specified particle
/// @param [in] i Particle index
/// @return Constant reference to embedding term
  
inline const double& CellListBase::emb(const uint i) const {
  return m_eamStorage.emb(i);
}


#endif // __CELL_LIST_BASE_HPP_INCLUDED
