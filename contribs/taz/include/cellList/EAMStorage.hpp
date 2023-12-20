/// @file
/// @brief Definition of the class EAMStorage

#ifndef __EAM_STORAGE_INCLUDED
#define __EAM_STORAGE_INCLUDED


#include <vector>

  #include<cassert>
/// @brief Class to store data specific to EAM potential
class EAMStorage {

public:

  /// @brief Default constructor
  EAMStorage() {}

  /// @brief Destructor (nothing to do)
  ~EAMStorage() {}

  /// @brief Get the size of the EAM storage
  /// @return Size
  inline uint size() { return m_rho.size(); }

  /// @brief Get a pointer to rho terms storage
  /// @return Pointer to rho terms storage
  inline double* rho() { return m_rho.data(); }

  /// @brief Get a pointer to embedding terms storage
  /// @return Pointer to embedding terms storage
  inline double* emb() { return m_emb.data(); }

  /// @brief Get a reference to rho terms of a specified particle
  /// @param [in] i Particle index
  /// @return Reference to rho terms
  inline double& rho(const uint i) { return m_rho[i]; }

  /// @brief Get a reference to embedding terms of a specified particle
  /// @param [in] i Particle index
  /// @return Reference to embedding terms
  inline double& emb(const uint i) { return m_emb[i]; }

  /// @brief Accessor to rho terms of a specified particle
  /// @param [in] i Particle index
  /// @return Constant reference to rho terms
  inline const double& rho(const uint i) const { return m_rho[i]; }

  /// @brief Accessor to embedding terms of a specified particle
  /// @param [in] i Particle index
  /// @return Constant reference to embedding terms
  inline const double& emb(const uint i) const { return m_emb[i]; }

  /// @brief Copy the embedding terms storage into another EAM storage
  /// @param [in] to Storage for the copy
  inline void embCopy(EAMStorage& to) const {
    to.m_emb = m_emb;
  }
   

    /// @brief Copy the embedding terms storage into another EAM storage
  /// @param [in] to Storage for the copy
  inline void embCopy_additional_storage(EAMStorage& to, int shift, int numberOfElements) const {
  
    assert(m_emb.size() >= shift + numberOfElements);
  
    const int oldSize = to.size();
    to.check(numberOfElements + oldSize);
    std::copy ( m_emb.begin()+shift, //first elem
                m_emb.begin()+shift+numberOfElements, //last elem
                to.emb()+oldSize ); // output 
  }

  /// @brief Resize the EAM storage
  /// @param [in] n New size
  inline void check(uint n) {
    m_rho.resize(n);
    m_emb.resize(n);
  }

  /// @brief Reset the EAM storage
  inline void reset() { 
    m_rho.clear(); // I added clear function for AMR
    m_emb.clear();
    m_rho.assign(m_rho.size(), 0);
    m_emb.assign(m_emb.size(), 0);
  }

//private:
  public :
  std::vector<double> m_rho; ///< Rho term for each particle of the cell
  std::vector<double> m_emb; ///< Embedding term for each particle of the cell

};

#endif // __EAM_STORAGE_INCLUDED
