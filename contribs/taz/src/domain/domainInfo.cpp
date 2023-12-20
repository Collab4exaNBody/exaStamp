/// @file 
/// @brief Implementation of GlobalInfo class and Configuration specialized to DomainInterface


#include <iomanip>
#include <iostream>

#include "globals.hpp"
#include "referenceMap.hpp"

#include "io/input.hpp"

#include "domain/domainInfo.hpp"


/// @brief A constructor from the program input data
/// @param [in] input Input to convert
Configuration<DomainInterface>::Configuration(const Input& input)
  : boundaryConditions(PERIODIC) {

  // Simulation mode and submode
  switch (input.mode) {

  case (Input::SINGLE_SPEC_ATOM) :
    mode = SINGLE;
    submode = ATOM;
    break;

  case (Input::SINGLE_SPEC_MESO) :
    mode = SINGLE;
    submode = MESO;
    break;

  case (Input::SINGLE_SPEC_SMOOTH) :
    mode = SINGLE;
    submode = SMOOTH;
    break;

  case (Input::SINGLE_SPEC_STIFF_MOLEC) :
    mode = SINGLE;
    submode = STIFF_MOLECULE;
    break;

  case (Input::ATOM_IN_MOL) :
    mode = BONDED;
    submode = ATOM_CHARGED;
    break;

  }

  // Boundary conditions
  for (int dim=0; dim<VEC3_NDIMS; ++dim) {
    switch (input.boundaryConditions[dim]) {

    case (Input::FREE) :
      boundaryConditions[dim] = FREE;
      break;

    case (Input::PERIODIC) :
      boundaryConditions[dim] = PERIODIC;
      break;

    case (Input::WALL) :
      boundaryConditions[dim] = WALL;
      break;

    case (Input::FREE_WALL) :
      boundaryConditions[dim] = FREE_WALL;
      break;

    case (Input::WALL_FREE) :
      boundaryConditions[dim] = WALL_FREE;
      break;

    }
  }

  // Domain limits
  origin    = input.origin;
  extension = input.extension;

  // Decomposition type
  switch (input.decomposition) {

  case (Input::RECTILINEAR) :
    decoupage = RECTILINEAR;
    numberOfDomainsPerDim.x = input.decoupage.x;
    numberOfDomainsPerDim.y = input.decoupage.y;
    numberOfDomainsPerDim.z = input.decoupage.z;
    break;

  case (Input::ANY) :
    numberOfDomainsPerDim.x = input.decoupage.x;
    numberOfDomainsPerDim.y = input.decoupage.y;
    numberOfDomainsPerDim.z = input.decoupage.z;
    decoupage = ANY;
    break;

  }

  // 
  symmetrization = input.symmetrization;
  
  // AMR
  Dmax=input.Dmax;
}


/// @brief Default constructor
///
/// Set all boundary conditions to periodic
GlobalInfo::GlobalInfo() : boundaryConditions(Configuration<DomainInterface>::PERIODIC) {
}


/// @brief Destructor
///
///
GlobalInfo::~GlobalInfo() {
  delete decomposition;
}


/// @brief Configure the domains common data
/// @param [in] configuration Configuration of the domains
/// @param [in] numberOfNodes Number of nodes
/// @param [in] rank Node rank
void GlobalInfo::configure(Configuration<DomainInterface>& configuration, int numberOfNodes, int rank) {

  mode    = configuration.mode;
  submode = configuration.submode;

  minBounds = configuration.origin;
  extension = configuration.extension;
  maxBounds = minBounds + extension;

  boundaryConditions = configuration.boundaryConditions;

  ghostThickness = Global::reference.getGhostThickness();

  double maxRcut;
  if(Global::reference.isMEAM())  maxRcut = 1.*Global::reference.getMaxRcut();
  else  maxRcut = Global::reference.getMaxRcut();
  if (maxRcut<0) maxRcut=1; // This is arbitrary, but if there is no short potential, cell size should not be negative


// là je fait un truc moche en attendant pour l'AMR avec 2 niveaux de refinement soit possible;
  uint levelOfRefinementMax=configuration.Dmax;
  
  vec3<int> tmp = auxFloor<int>(vec3<double>(extension/( (1<<levelOfRefinementMax) *maxRcut)));

  numberOfCellsPerDim.x = tmp.x;
  numberOfCellsPerDim.y = tmp.y;
  numberOfCellsPerDim.z = tmp.z;

  numberOfCells = product(numberOfCellsPerDim);

  cellLength = extension/tmp;

  decomposition = Decompose(configuration, numberOfNodes, rank);

}


/// @brief Print global domains data in specified flux as a check
/// @param [in,out] flux Print flux
void GlobalInfo::print(std::ostream& flux) {
  
  using namespace std;

  const vec3<double>& min = getMinBounds();
  const vec3<double>& max = getMaxBounds();
  
  // Utility function to print boundary condition enum
  auto lambda = [&](BoundaryCondition bc) -> void { 
    switch (bc) {
    case Configuration<DomainInterface>::WALL     : flux<< setw(12) << right << "wall";      break;
    case Configuration<DomainInterface>::FREE     : flux<< setw(12) << right << "free";      break;
    case Configuration<DomainInterface>::PERIODIC : flux<< setw(12) << right << "periodic";  break;
    case Configuration<DomainInterface>::FREE_WALL: flux<< setw(12) << right << "free_wall"; break;
    case Configuration<DomainInterface>::WALL_FREE: flux<< setw(12) << right << "wall_free"; break;
    }
  };

  //
  flux<< "DOMAIN " << endl;

  flux<< "  " << setw(22) << left << "Min global bounds"   << " : " << "[" << scientific << setw(12) << min << "]"  << endl
      << "  " << setw(22) << left << "Max global bounds"   << " : " << "[" << scientific << setw(12) << max << "]"  << endl
      << "  " << setw(22) << left << "Boundary conditions" << " : " << "[";
 
  lambda(boundaryConditions.x); flux<< " "; 
  lambda(boundaryConditions.y); flux<< " "; 
  lambda(boundaryConditions.z); 
  flux<< "]" << endl << endl;

  //
  flux<< "DECOMPOSITION " << endl;

  decomposition->print(flux);

  flux<< "  " << setw(22) << left << "Cells"           << " : " << "[" << numberOfCellsPerDim << "] " << numberOfCells  << endl
      << "  " << setw(22) << left << "Cell size"       << " : " << "[" << scientific << cellLength << "]"  << endl
      << "  " << setw(22) << left << "Ghost thickness" << " : " << ghostThickness << endl
      << endl;
}
