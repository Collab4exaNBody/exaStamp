/// @file 
/// @brief Interface for GlobalInfo class and Configuration specialized to DomainInterface

#ifndef __DOMAIN_INFO_HPP_INCLUDED
#define __DOMAIN_INFO_HPP_INCLUDED


#include "domain/decomposition.hpp"

#include "utils/constants.hpp"
#include "utils/vec3/vec3.hpp"


class Decomposition;
class DomainInterface;
class Input;



template <class T> class Configuration;


/// @brief Temporary structure gathering all data to initialize the domains
template <> class Configuration<DomainInterface> {

public:

  /// @brief Enumeration of the possible "modes"
  enum Mode    {
	  SINGLE, ///< Only one subclass of particles
	  BONDED ///< One class of point paticles bonded into molecules
  };
  /// @brief Enumeration of the possible "submodes"
  enum SubMode {
	  ATOM, ///< Particles are atoms
	  ATOM_CHARGED, ///< Particles are atoms and have charges
	  MESO, ///< Particles are mesoparticles
	  SMOOTH, ///< Particles are smooth particles	
	  STIFF_MOLECULE ///< Particles are stiff molecules
  };

  /// @brief Enumeration of the various boundary conditions
  enum BoundaryCondition {
	  FREE, ///< Free boundary condition
	  PERIODIC, ///< Periodic boundary condition
	  WALL, ///< Wall boundary condition
	  FREE_WALL, ///< Free on domain lower bound and wall on upper bound
	  WALL_FREE ///< Wall on domain lower bound and free on upper bound
	  };
  /// @brief Shortcut for a vector of enum type BondaryCondition
  typedef vec3<BoundaryCondition> BoundaryConditions;

  /// @brief Enumeration of the splitting possibilities
  enum Decoupage {
	  RECTILINEAR, ///< Rectilinear splitting
	  ANY ///< Non-rectilinear splitting
  };


  /// @brief Default constructor
  ///
  ///
  Configuration() : boundaryConditions(PERIODIC) {}
  /// @brief Destructor (nothing to do)
  ~Configuration() {}

  Configuration(const Input& input);

  Mode mode; ///< Domains "mode" (single subclass of particle or other)
  SubMode submode; ///< Domains "submode" (subclass of particle used for the mode "single")

  vec3<BoundaryCondition> boundaryConditions; ///< Boundary condition for each dimension

  vec3<double> origin; ///< System origin
  vec3<double> extension; ///< System extension

  vec3<uint> numberOfDomainsPerDim; ///< Number of domain for each dimension

  Decoupage decoupage; ///< Space splitting type (rectilinear or any)

  bool symmetrization; ///< Allow symmetrization of force calculation
  
  // AMR 
  int Dmax;


};

/// @brief Shortcut for enum type BoudaryCondition (from Configuration<DomainInterface>)
typedef Configuration<DomainInterface>::BoundaryCondition  BoundaryCondition;
/// @brief Shortcut for type BoudaryConditions (from Configuration<DomainInterface>)
typedef Configuration<DomainInterface>::BoundaryConditions BoundaryConditions;


/// @brief Class to store data common to all domains
class GlobalInfo {

protected:

  /// @brief Shortcut for enum type Mode (from Configuration<DomainInterface>)
  typedef Configuration<DomainInterface>::Mode Mode;
  /// @brief Shortcut for enum type SubMode (from Configuration<DomainInterface>)
  typedef Configuration<DomainInterface>::SubMode SubMode;

public:

  GlobalInfo();
  ~GlobalInfo();

  void configure(Configuration<DomainInterface>& configuration, int numberOfNodes, int rank);

  Mode    getMode();
  SubMode getSubMode();

  const vec3<double>& getMinBounds();
  const vec3<double>& getMaxBounds();
  const vec3<double>& getExtension();

  double getVolume();
  
  const vec3<double>& getCellLength();

  double getCellVolume();
  double getWorkloadPrefactor();
  
  const vec3<int>& getNumberOfCellsPerDim();
  const BoundaryConditions& getBoundaryConditions();

  Decomposition* getDecomposition();  

  uint getNumberOfDomains();

  bool noFreeBoundaryConditions();

  vec3<bool> isFreeInf();
  vec3<bool> isFreeSup();

  void print(std::ostream& flux);


private:


  Mode mode; ///< Domains "mode" (single subclass of particle or other)
  SubMode submode; ///< Domains "submode" (subclass of particle used for the mode "single")

  vec3<double> minBounds; ///< System lower bound in each dimension
  vec3<double> maxBounds; ///< System upper bound in each dimension
  vec3<double> extension; ///< System extension in each dimension

  BoundaryConditions boundaryConditions; ///< Boundary conditions in each dimension

  vec3<double> cellLength; ///< Cell dimensions

  uint ghostThickness; ///< Ghost thickness
  uint numberOfCells; ///< Total number of cells

  vec3<int> numberOfCellsPerDim; ///< Number of cells in each dimension

  Decomposition* decomposition; ///< Specify the decomposition of the system in domains

};


/// @brief Accessor to mode
inline Configuration<DomainInterface>::Mode GlobalInfo::getMode() {
  return mode;
}


/// @brief Accessor to submode
inline Configuration<DomainInterface>::SubMode GlobalInfo::getSubMode() {
  return submode;
}


/// @brief Accessor to minBounds
inline const vec3<double>& GlobalInfo::getMinBounds() {
  return minBounds;
}


/// @brief Accessor to extension
inline const vec3<double>& GlobalInfo::getExtension() {
 return extension;
}


/// @brief Compute the domain volume
/// @return Total volume
inline double GlobalInfo::getVolume() {
  return product(extension);
}


/// @brief Accessor to maxBounds
inline const vec3<double>& GlobalInfo::getMaxBounds() {
 return maxBounds;
}


/// @brief Accessor to cellLength
inline const vec3<double>& GlobalInfo::getCellLength() { 
  return cellLength;
}


/// @brief Compute the cell volume
/// @return Cell volume
inline double GlobalInfo::getCellVolume() {
  return product(cellLength);
}


/// @brief Compute the prefactor for the potential induced workload
inline double GlobalInfo::getWorkloadPrefactor(){
	return 4*ConstantMath::pi/(3*product(cellLength));
}


/// @brief Accessor to numberOfCellsPerDim
inline const vec3<int>& GlobalInfo::getNumberOfCellsPerDim() { 
  return numberOfCellsPerDim;
}


/// @brief Accessor to boundaryConditions
inline const BoundaryConditions& GlobalInfo::getBoundaryConditions() { 
  return boundaryConditions;
}


/// @brief Accessor to decomposition
inline Decomposition* GlobalInfo::getDecomposition() {
  return decomposition;
}


/// @brief Return the number of domains (from decomposition)
/// @return The number of domains
inline uint GlobalInfo::getNumberOfDomains() { 
  return decomposition->getNumberOfDomains();
}


/// @brief Check if there is free boundary conditions
/// @return False if there is free boundary conditions
inline  bool GlobalInfo::noFreeBoundaryConditions() { 
  return (boundaryConditions[0]==Configuration<DomainInterface>::PERIODIC || boundaryConditions[0]==Configuration<DomainInterface>::WALL) && 
         (boundaryConditions[1]==Configuration<DomainInterface>::PERIODIC || boundaryConditions[1]==Configuration<DomainInterface>::WALL) &&
         (boundaryConditions[2]==Configuration<DomainInterface>::PERIODIC || boundaryConditions[2]==Configuration<DomainInterface>::WALL);
}


/// @brief For each lower boundary, check free conditions
/// @return True for each free lower boundary (in a vector)
inline vec3<bool> GlobalInfo::isFreeInf() {
  return vec3<bool>(boundaryConditions[0]==Configuration<DomainInterface>::FREE || boundaryConditions[0]==Configuration<DomainInterface>::FREE_WALL, 
		    boundaryConditions[1]==Configuration<DomainInterface>::FREE || boundaryConditions[1]==Configuration<DomainInterface>::FREE_WALL, 
		    boundaryConditions[2]==Configuration<DomainInterface>::FREE || boundaryConditions[2]==Configuration<DomainInterface>::FREE_WALL);
}


/// @brief For each upper boundary, check free conditions
/// @return True for each free upper boundary (in a vector)
inline vec3<bool> GlobalInfo::isFreeSup() {
  return vec3<bool>(boundaryConditions[0]==Configuration<DomainInterface>::FREE || boundaryConditions[0]==Configuration<DomainInterface>::WALL_FREE, 
		    boundaryConditions[1]==Configuration<DomainInterface>::FREE || boundaryConditions[1]==Configuration<DomainInterface>::WALL_FREE, 
		    boundaryConditions[2]==Configuration<DomainInterface>::FREE || boundaryConditions[2]==Configuration<DomainInterface>::WALL_FREE);
  
}

#endif // __DOMAIN_INFO_HPP_INCLUDED
