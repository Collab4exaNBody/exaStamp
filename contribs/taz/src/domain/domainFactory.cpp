/// @file
/// @brief Implementation of domain building functions


#include "globals.hpp"

#include "parallel/types/MPI_inMol.hpp"
#include "parallel/types/MPI_mesoparticle.hpp"
#include "parallel/types/MPI_smoothparticle.hpp"
#include "parallel/types/MPI_atom.hpp"

#include "domain/domain.hpp"

#include "grid/AMRGrid/AMRGrid.hpp"


/// @brief Shortcut for class Configuration<DomainInterface>
typedef Configuration<DomainInterface> Config;


class Node;

/// @brief Build a domain from domains common info and domains configuration
/// @param [in,out] node The node (pointer)
/// @param [in] index Domain index
/// @param [in] config Domains configuration
/// @param [in] numberOfParticles Number of particles (not used)
/// @return Build domain
DomainInterface* buildDomain(Node* node, uint index, Config& config, uint64_t numberOfParticles) {

  DomainInterface* domain = nullptr;

  auto mode    = Global::domainInfo.getMode();
  auto subMode = Global::domainInfo.getSubMode();

  if(index == 0) std::cout << "                      with an Adaptive Mesh Refinement " ;

        switch (mode) {

  case Config::SINGLE :
      
    switch (subMode) {
	
    case Config::ATOM :
      
      switch (config.decoupage) {
	
      case (Config::RECTILINEAR) :

        domain = new Domain
	  < AMRGrid
	    < RectilinearGridInfo,
	      Atom >> (node, index, config, numberOfParticles);
	//
      	break; // end case rectilinear

      case (Config::ANY) :
	
	domain = new Domain
	  < AMRGrid
	    <  AnyGridInfo, 
	      Atom
	    >> (node, index, config, numberOfParticles);
	
      	break; // end case any
            
      } // END SWITCH DECOUPAGE
      break; // end case subMode = atom

      
    case Config::ATOM_CHARGED :
      std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": "
	       << "in function 'DomainInterface* buildDomain(Node*, uint, Configuration<DomainInterface>&, uint64_t)' : charged atoms are not handled yet. STOP."
	       << std::endl;
      exit(-1);
      break; // end case subMode = atom_charged

      
    case Config::MESO :
      std::cerr << "ERROR: You must compile ExaStamp with mesoparticles enabled!" << std::endl;
      exit(1);
      break; // end case subMode = meso


    case Config::SMOOTH :
      std::cerr << "ERROR: You must compile ExaStamp with smoothparticles enabled!" << std::endl;
      exit(1);
      break; // end case subMode = smooth
      

    case Config::STIFF_MOLECULE :
      std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": "
	       << "in function 'DomainInterface* buildDomain(Node*, uint, Configuration<DomainInterface>&, uint64_t)' : stiff molecules are not handled yet. STOP."
	       << std::endl;
      exit(-1);
      break; // end case subMode = stiff_molecule

    } // END SWITCH SUBMODE
    break; 



  case Config::BONDED :
    break; // end case mode = bonded
	
  } // END SWITCH MODE
  
  return domain;
      
}

