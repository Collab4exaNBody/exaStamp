/// @file 
/// @brief Definition of the tools to initialize the particles

#ifndef __LATTICE_HPP_INCLUDED
#define __LATTICE_HPP_INCLUDED


#include <vector>
#include <string>

#include "referenceMap.hpp"

#include "domain/domainInfo.hpp"

#include "io/outputManager.hpp"

#include "parallel/types/MPI_mesoparticle.hpp"
#include "parallel/types/MPI_smoothparticle.hpp"

#include "saru/saru.hpp"
#include "utils/stampUnits.hpp"


template <class T> class Configuration;
class Input;


/// @brief Tool to initialize a lattice of particles
struct Lattice {

  /// @brief Constructor
  /// @param [in] nParameters Number of space parameters used to build the lattice
  /// @param [in] nAtoms Number of atom per cell of the lattice
  Lattice(int nParameters, int nAtoms)
    : numberOfParameters(nParameters),
      numberOfAtoms(nAtoms), 
      parameters(numberOfParameters),
      types(numberOfAtoms),
      base(numberOfAtoms) {}

  /// @brief Destructor (nothing to do)
  ~Lattice() {}

  int numberOfParameters; ///< Number of space parameters used to build the lattice
  int numberOfAtoms; ///< Number of atom per cell of the lattice

  Array<double> parameters; ///< Space parameters used to build the lattice
  Array<int> types; ///< Types of the atoms in the lattice
  Array< vec3<double> > base; ///< Position of the atoms in a cell of the lattice

};

Lattice* initFCCLattice(double parameter, Array<int> types_);
Lattice* initBCCLattice(double parameter, Array<int> types_);
Lattice* initDIAM100Lattice(double parameter, Array<int> types_);
Lattice* initSCLattice(double parameter, Array<int> types_);


/// @brief Temporary structure gathering all data to initialize particles
template <> class Configuration<Particle> {

public:

  /// @brief Default constructor
  Configuration() {}

  /// @brief Destructor
  ///
  ///
  ~Configuration() {
    if (lattice!=nullptr) 
      delete lattice;
  }

  Configuration(const Input& input);

  /// @brief Accessor to the initialization type
  InitType getInitType() { return m_initType; }
  /// @brief  Accessor to the name of the file containing particles positions
  const std::string& getFilename() { return filename; }
  /// @brief Accessor to the initial step
  int getInitStep() { return m_initStep; }
  /// @brief Accessor to the id of the file containing particles positions
  InputOutputManager::FileId& getFileId() { return fileId; }
  /// @brief Accessor to the initial internal temperature
  double getInitTint() { return initialTint; }

  uint64_t totalNumberOfParticles();

  /// @brief Get the total number of cells in the lattice
  /// @return Number of cells or -1 if there is no lattice
  uint64_t totalNumberOfLatticeCells() {
    if (m_initType==INIT_DEFAULT)
      return (uint64_t) NCELLS.x * (uint64_t) NCELLS.y * (uint64_t) NCELLS.z; //product(NCELLS);
    return -1;
  }

  /// @brief Get the number of atom per cell
  /// @return Number of atom per cell
  uint numberOfAtomPerCell() {
    return lattice->numberOfAtoms;
  }

  template <class Exchange_P>
  void buildSubset(const uint64_t begin, const uint64_t end, std::vector<Exchange_P>& particles, const vec3<double>& minBounds, const vec3<double>& maxBounds) const;

public:

  InitType m_initType; ///< Initialization type
  int m_initStep; ///< Initial step

  // Used for hard init
  vec3<int> numberOfLatticeCells; ///< Number of cells of the lattice in each dimension asked in the input
  double initialTemperature; ///< Initial temperature
  double initialTint; ///< Initial internal temperature
  Lattice* lattice; ///< Lattice
  vec3<int> NMAX; ///< Maximal number of cells of the lattice in each dimension if the extension is given
  vec3<int> NCELLS; ///< Effective number of cells of the lattice in each dimension

  vec3<double> wallMinBounds; ///< Walls lower bound 
  vec3<double> wallMaxBounds; ///< Walls upper bound 
  vec3<std::string> wallLowerTypes; ///< Walls lower types 
  vec3<std::string> wallUpperTypes; ///< Walls upper types 
  
  // Used for legacy init
  std::string filename; ///< Name of the file containing particles positions
  InputOutputManager::FileId fileId; ///< Id of the file containing particles positions

  uint64_t numberOfParticles; ///< Number of particles

};


#include "particle/lattice.hxx"


#endif // __LATTICE_HPP_INCLUDED
