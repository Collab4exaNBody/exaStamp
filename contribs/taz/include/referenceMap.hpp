/// @file 
/// @brief Tools to manage the potentials and particle types

#ifndef __REFERENCE_MAP_HPP_INCLUDED
#define __REFERENCE_MAP_HPP_INCLUDED


#include <map>
#include <vector>

#include "particle/types/typeAtom.hpp"
#include "particle/types/typeMesoparticle.hpp"
#include "particle/types/typeSmoothParticle.hpp"
#include "particle/types/typeSmoothWallParticle.hpp"

#include "potential/potential.hpp"
#include "eos/eos.hpp"
#include "utils/kernelFunction/kernelFunction.hpp"
#include "chemistry/allKinetics.hpp"

class Input;
template <class T> class Configuration;
class MPI__InMol;




/// @brief Fix maximum ghost thickness
#define MAX_GHOST_THICKNESS 3


/// @brief Structure used to store a potential and two particle types
///
/// The force computation work with typed potentials
struct TypedPotential {

  /// @brief Constructor
	/// @param [in] typeIndexA First particle type
	/// @param [in] typeIndexB Second particle type
	/// @param [in] ipotential Potential
  TypedPotential(uint8_t typeIndexA, uint8_t typeIndexB, IPotential* ipotential)
    : first(typeIndexA), second(typeIndexB), 
      potential(ipotential),
      traversal(potential->getTraversal()), type(potential->getType()) {
  }

  /// @brief Destructor (nothing to do)
  ~TypedPotential() {}

  uint8_t first;         ///< Index of the first type
  uint8_t second;        ///< Index of the second type

  IPotential* potential; ///< Pointer to the potential

  IPotential::Traversal traversal; ///< Subclass of the potential
  IPotential::Type type; ///< Type of the potential

};



/// @brief Structure used to store an EOS and an associate type
/// 
struct TypedEOS {

  /// @brief Arguments constructor
  TypedEOS(uint8_t typeIndex, IEOS* ieos)
    : type(typeIndex),
      eos(ieos),
      typeEOS(eos->getType()) {
  }

  /// @brief Destructor
  ~TypedEOS() {}

  uint8_t type;         ///< index of the type

  IEOS* eos; ///< pointer to the eos

  IEOS::Type typeEOS; ///< EOS type

};



/// @brief Structure used to store a chemical reaction and an associate type
/// 
struct TypedReaction {

  /// @brief Arguments constructor
  TypedReaction(uint8_t typeIndex, IKinetics* ikinetics)
    : type(typeIndex),
      reaction(ikinetics),
      typeReaction(ikinetics->getType()) {
  }

  /// @brief Destructor
  ~TypedReaction() {}

  uint8_t type; ///< index of the type

  IKinetics* reaction; ///< pointer to the reaction kinetics

  IKinetics::Type typeReaction; ///< Kinetics type

};



/// @brief Class to manage potentials and types as global variables
class ReferenceMap {

public:

  ReferenceMap();
  ~ReferenceMap();

  void configure(const Configuration<TypeParticle>& configuration);
  void configure(const Configuration<IPotential>& configuration);
  void configure(Configuration<MPI__InMol>& configuration, const Configuration<IPotential>& potentials);
  void configure(const Configuration<IEOS>& configuration);
  void configure(const Configuration<IKinetics>& configuration);
  void configure(const Input& input);

  TypeParticle* find(uint8_t typeindex);
  TypeParticle* find(std::string name);

  uint8_t findTypeIndex(std::string name);
  std::string findSubtype(uint8_t id);

  uint    getGhostThickness();
  double  getMaxRcut();
  uint8_t getNumberOfTypes();

  bool    isEAM();
  bool    isMEAM();
  bool    isDPD();
  bool    isSDPD();
  bool    isReactive();
  bool    isTimers();
  bool    isVerlet();
  bool    isBlockVerlet();
  void setIsBlockVerlet(); 

  double getMass   (uint8_t typeindex);
  double getInvMass(uint8_t typeindex);

  IPotential* getPotential(uint8_t typeindex, uint8_t typejndex);

  double getRcut (uint8_t typeindex, uint8_t typejndex);
  double getCost (uint8_t typeindex, uint8_t typejndex);
  double getRcut2(uint8_t typeindex, uint8_t typejndex);
  void   getRcut (uint8_t typeindex, const uint8_t* typejndex, double* rcuts, uint n);
  void   getRcut2(uint8_t typeindex, const uint8_t* typejndex, double* rcuts, uint n);

  IEOS* getEOS(uint8_t typeindex);
  IKinetics* getReaction(uint8_t typeindex);
  double getExothermicity(uint8_t typeIndex);
  
  KernelFunction::Type getKernelType(uint8_t typeindex);
  KernelFunction* getKernel(uint8_t typeindex);

  double getGammaParallel(uint8_t typeindex);
  double getGammaOrthogonal(uint8_t typeindex);

  vec3<double> getVelocity(uint8_t typeindex);

  bool isWallType(uint8_t typeindex);
    
  /// @brief Loop through all potentials and run parameter function
  /// @tparam Func Function type
  /// @param [in] func Function to run
  template<typename Func> void for_all_potentials(Func func) {
  
    for (auto it=potentialReference.begin(); it!=potentialReference.end(); ++it) {
      func(TypedPotential(it->first.first, it->first.second, it->second));
    }
  }

  /// @brief Loop through all eos and run parameter function
  /// @tparam Func Function type
  /// @param [in] func Function to run
  template<typename Func> void for_all_eos(Func func) {
    for (auto it=eosReference.begin(); it!=eosReference.end(); ++it) {
      func(TypedEOS(it->first, it->second));
    }
  }

  /// @brief Loop through all dpd interactions and run parameter @c func
  /// @param func can be a function pointer, a functor or a lambda function
  template<typename Func> void for_all_dpd(Func func) {
    for (auto it=dpdReference.begin(); it!=dpdReference.end(); ++it) {
      func(std::pair<uint8_t,uint8_t>(it->first, it->second));
    }
  }

  /// @brief Loop through all dpde interactions and run parameter @c func
  /// @param func can be a function pointer, a functor or a lambda function
  template<typename Func> void for_all_dpde(Func func) {
    for (auto it=dpdeReference.begin(); it!=dpdeReference.end(); ++it) {
      func(std::pair<uint8_t,uint8_t>(it->first, it->second));
    }
  }

  /// @brief Loop through all sdpd interactions and run parameter @c func
  /// @param func can be a function pointer, a functor or a lambda function
  template<typename Func> void for_all_sdpd(Func func) {
    for (auto it=sdpdReference.begin(); it!=sdpdReference.end(); ++it) {
      func(std::pair<uint8_t,uint8_t>(it->first, it->second));
    }
  }

  /// @brief Loop through all sdpd wall interactions and run parameter @c func
  /// @param func can be a function pointer, a functor or a lambda function
  template<typename Func> void for_all_sdpd_wall(Func func) {
    for (auto it=sdpdWallReference.begin(); it!=sdpdWallReference.end(); ++it) {
      func(std::pair<uint8_t,uint8_t>(it->first, it->second));
    }
  }

  /// @brief Loop through all chemical reactions and run parameter function
  /// @tparam Func Function type
  /// @param [in] func Function to run
  template<typename Func> void for_all_reaction(Func func) {
    for (auto it=reactionReference.begin(); it!=reactionReference.end(); ++it) {
      func(TypedReaction(it->first, it->second));
    }
  }
  

  void print(std::ostream& flux);

  /// @brief Interaction types
  enum InteractionType {
    DPD, ///< DPD
    DPDE, ///< DPDE
    SDPD, ///< SDPD
    SDPD_WALL ///< SDPD with walls
  };

private:

  /// @brief Shortcut for a map that link an index to an atom type
  typedef std::map<uint8_t, TypeAtom> TypeAtomMap;

  /// @brief Shortcut for a map that link an index to a mesoparticle type
  typedef std::map<uint8_t, TypeMesoparticle> TypeMesoMap;

  /// @brief Shortcut for a map that link an index to a smooth mesoparticle type
  typedef std::map<uint8_t, TypeSmoothParticle> TypeSmoothMap;

  /// @brief Shortcut for a map that link an index to a virtual wall type
  typedef std::map<uint8_t, TypeSmoothWallParticle> TypeWallMap;

  /// @brief Shortcut for a map that link two atom types to a potential
  typedef std::map<std::pair<uint8_t, uint8_t>, IPotential* > PotentialMap; ///< Map linking pairs of atom types with a potential

  /// @brief Shortcut for a map that link an index to an equation of state
  typedef std::map<uint8_t, IEOS*> EOSMap;

  /// @brief Shortcut for a map that link an index to a kernel function
  typedef std::map<uint8_t, KernelFunction*> KernelMap;

  /// @brief Shortcut for a vector storing a pair of particle types 
  typedef std::vector< std::pair<uint8_t,uint8_t> > InteractionMap;

  /// @brief Shortcut for a map that link an index to a chemcial reaction
  typedef std::map<uint8_t, IKinetics*> ReactionMap;

  /// @brief Enumeration of the classes of particles
  enum ParticleType {
  	ATOM,		///< Atom
  	STIFF_MOLECULE,	///< Stiff molecule
  	SOFT_MOLECULE,	///< Soft molecule (not implemented yet)
  	MESOPARTICLE, ///< Mesoparticle
  	SMOOTH_PARTICLE, ///< Smooth particle
  	WALL_PARTICLE, ///< Smooth particle
  	NOT_FOUND=404	///< No class found
  };


  void setGhostThickness();

  ParticleType findType(uint8_t typeIndex);
  ParticleType findType(std::string name);

  void Register(TypeAtom& type);
  void RegisterNoWarn(TypeAtom& type);
  void UnRegister(TypeAtom& type);

  void Register  (TypeMesoparticle& type);
  void UnRegister(TypeMesoparticle& type);

  void Register  (TypeSmoothParticle& type);
  void UnRegister(TypeSmoothParticle& type);

  void Register  (TypeSmoothWallParticle& type);
  void UnRegister(TypeSmoothWallParticle& type);

  void Register(IPotential* potential, std::string nameA, std::string nameB);
  
  void Register(IEOS* eos, std::string name);
  
  void Register(KernelFunction* eos, std::string name);
  
  void Register(InteractionType itype, std::string nameA, std::string nameB);

  void Register(IKinetics* reaction, std::string name);

  
  bool checkReactiveEOS (Array< std::pair<IEOS*,IEOS*> >& store, IEOS* eos, const Array<std::string>& reactiveEOS0Name, const Array<std::string>& reactiveEOS1Name, const std::string& particleType, const Array<std::string>& reactiveParticleType);
  IEOS* createReactiveEOS(IEOS::Type EOS0, IEOS::Type EOS1, double size, double mass, IEOS* eos0, IEOS* eos1);
  template <class EOS0> IEOS* __createReactiveEOS(IEOS::Type EOS1, double size, double mass, EOS0* eos0, IEOS* eos1);
  
  uint8_t nextIndex; ///< Next index while registering the atom types then number of atom types
  std::map<std::string, uint8_t> indexReference; ///< Map that link the name of the atoms to the index of the atom types
  TypeAtomMap  atomReference; ///< Map that link an index to an atom type
  TypeMesoMap mesoReference; ///< Map that link an index to a mesoparticle type
  TypeSmoothMap smoothReference; ///< Map that link an index to a smooth particle type
  TypeWallMap wallReference; ///< Map that link an index to a virtual wall particle type
  PotentialMap potentialReference; ///<  Map that link two atom types to a potential
  std::map<uint8_t,std::string> subtypesReference; ///< Map that link a subtype and its index
  EOSMap eosReference; ///< Map that link an index to an equation of state
  KernelMap kernelReference; ///< Map that link an index to a kernelReference function
  InteractionMap dpdReference; ///< Vector storing the pair of particles interacting with the DPD FD forces
  InteractionMap dpdeReference; ///< Vector storing the pair of particles interacting with the DPDE FD forces
  InteractionMap sdpdReference; ///< Vector storing the pair of particles interacting with the SDPD forces (conservative + FD)
  InteractionMap sdpdWallReference; ///< Vector storing the pair of particles interacting with virtual SDPD wall particles
  ReactionMap reactionReference; ///< Map that link an index to a chemical reaction
  uint ghostThickness; ///< Ghost thickness for the system

  bool isPotentialEAM; ///< Indicates if some of the potentials are EAM
  bool isPotentialMEAM; ///< Indicates if some of the potentials are MEAM
  bool existsDPD; ///< Indicates if some interactions are of DPD or DPDE type
  bool existsSDPD; ///< Indicates if some interactions are of SDPD type

  bool isLocalTimers; ///< Indicate if the local timers are used

  bool existsReactive; ///< Indicates if there is a chemical reaction

  bool  isVERLET;  ///< Indicates if the verlet list method is used
  bool isBLOCKVERLET;
  double rVerlet; ///< raduis of Verlet
  double* rcutFastVerlet; ///< Cutoff radius for each pair of atom types
  double* rcutFast2Verlet; ///< Squared cutoff radius for each pair of atom types

  double  rcutMax; ///< Maximal cutoff radius
  double* rcutFast; ///< Cutoff radius for each pair of atom types
  double* rcutFast2; ///< Squared cutoff radius for each pair of atom types
  double* costFast; ///< Cost for each pair of atom types

  double* massFast; ///< Mass for each atom type
  double* massInvFast; ///< Inverse of the mass for each atom type

  double* gammaParaFast; ///< Friciton parameter in the parallel direction
  double* gammaOrthoFast; ///< Friction parameter in the orthogonal direction

  double* exothermicityFast; ///< Exothermicity for each mesoparticle type
  
  vec3<double>* velocityFast; ///< Velocity of the walls
  
  KernelFunction::Type* kernelType; ///< Kernel types

  public :

  void   neighbours_method(bool Verlet, double raduis);
  
  /// @brief return the valuis of raduis of Verlet
  double getrVerlet(){return rVerlet;};
  double getrVerlet (uint8_t typeindex, uint8_t typejndex);
  double getrVerlet2(uint8_t typeindex, uint8_t typejndex);
  void   getrVerlet (uint8_t typeindex, const uint8_t* typejndex, double* rcuts, uint n);
  void   getrVerlet2(uint8_t typeindex, const uint8_t* typejndex, double* rcuts, uint n);
  
#ifndef SingleGrid

  inline int getNumberOfPotential() {    
    return potentialReference.size(); 
 }
 
#endif

  // Adaptive Mesh Refinement
  
 public:
  double getDmax() {return Dmax;};
  int getAmrCriterion() {return r_amrCriterion;} 
  int Dmax;
  int r_amrCriterion;
};


// Some inlining to go faster ...


/// @brief Accessor to ghost thickness
inline uint ReferenceMap::getGhostThickness() {
  return ghostThickness;
}


/// @brief Accessor to maximum cutoff radius
inline double ReferenceMap::getMaxRcut() {
  return rcutMax;
}


/// @brief Get the cutoff radius between two types
/// @param [in] typeindex First type
/// @param [in] typejndex Second type
/// @return Cutoff radius
inline double ReferenceMap::getRcut(uint8_t typeindex, uint8_t typejndex) {

  // this is good programming but slow ...

  // IPotential* potential = this->getPotential(typeindex, typejndex);
  // ShortRangePotential* shortpot = dynamic_cast<ShortRangePotential*>(potential);
  // if (shortpot!=nullptr) return shortpot->getCutoffRadius();
  // else return 0;

  // this is bad but fast(er) : 
  return rcutFast[typeindex*nextIndex+typejndex];

}


/// @brief Get the cost between two types
/// @param [in] typeindex First type
/// @param [in] typejndex Second type
/// @return Cost
inline double ReferenceMap::getCost(uint8_t typeindex, uint8_t typejndex) {
  return costFast[typeindex*nextIndex+typejndex];
}


/// @brief Get the cutoff radii between a type and an array of types
/// @param [in] typeindex First type
/// @param [in] typejndex Array of second types
/// @param [out] rcuts Cutoff radii
/// @param [in] n Size of the array of types
inline void ReferenceMap::getRcut(uint8_t typeindex, const uint8_t* typejndex, double* rcuts, uint n) {
  double* tmp = &rcutFast[typeindex*nextIndex];
  for (uint j=0; j<n; ++j) 
    rcuts[j] = tmp[typejndex[j]];
}


/// @brief Get the square cutoff radius between two types
/// @param [in] typeindex First type
/// @param [in] typejndex Second type
/// @return Square cutoff radius
inline double ReferenceMap::getRcut2(uint8_t typeindex, uint8_t typejndex) {

  return rcutFast2[typeindex*nextIndex+typejndex];
}


/// @brief Get the square cutoff radii between a type and an array of types
/// @param [in] typeindex First type
/// @param [in] typejndex Array of second types
/// @param [out] rcuts Square cutoff raduis
/// @param [in] n Size of the array of types
inline void ReferenceMap::getRcut2(uint8_t typeindex, const uint8_t* typejndex, double* rcuts, uint n) {
  double* tmp = &rcutFast2[typeindex*nextIndex];

  for (uint j=0; j<n; ++j) 
    rcuts[j] = tmp[typejndex[j]];

}

/// @brief Get the cutoff radius between a type and an array of types with the Verlet lists method. If linked cells method is used, rcutFastVerlet = rcutFast
/// @param [in] typeindex First type
/// @param [in] typejndex Array of second types
/// @param [out] rcuts Cutoff raduis
/// @param [in] n Size of the array of types
inline void ReferenceMap::getrVerlet(uint8_t typeindex, const uint8_t* typejndex, double* rcuts, uint n) {
  double* tmp = &rcutFastVerlet[typeindex*nextIndex];
  for (uint j=0; j<n; ++j) 
    rcuts[j] = tmp[typejndex[j]];
}


/// @brief Get the square cutoff radius between two types with the Verlet lists method. If linked cells method is used, rcutFastVerlet2 = rcutFast2
/// @param [in] typeindex First type
/// @param [in] typejndex Second type
/// @return Square cutoff radius
inline double ReferenceMap::getrVerlet2(uint8_t typeindex, uint8_t typejndex) {

  return rcutFast2Verlet[typeindex*nextIndex+typejndex];
}

/// @brief Get the square cutoff radii between a type and an array of types with the Verlet lists method with the Verlet lists method. If linked cells method is used, rcutFastVerlet = rcutFast
/// @param [in] typeindex First type
/// @param [in] typejndex Array of second types
/// @param [out] rcuts Square cutoff radii
/// @param [in] n Size of the array of types
inline void ReferenceMap::getrVerlet2(uint8_t typeindex, const uint8_t* typejndex, double* rcuts, uint n) {
  double* tmp = &rcutFast2Verlet[typeindex*nextIndex];

  for (uint j=0; j<n; ++j) 
    rcuts[j] = tmp[typejndex[j]];

}



/// @brief Ge the mass for a type
/// @param [in] typeindex Index of the type
/// @return Mass
inline double ReferenceMap::getMass(uint8_t typeindex) {
  return massFast[typeindex];
}


/// @brief Get the mass inverse for a type
/// @param [in] typeindex Index of the type
/// @return Mass inverse
inline double ReferenceMap::getInvMass(uint8_t typeindex) {
  return massInvFast[typeindex];
}


/// @brief Get friction coefficient in the parallel direction \f$\gamma_{\parallel}\f$ for a type
/// @param [in] typeindex Index of the type
/// @return Friction coefficient in the parallel direction
inline double ReferenceMap::getGammaParallel(uint8_t typeindex) {
  return gammaParaFast[typeindex];
}


/// @brief Get friction coefficient in the orthogonal direction \f$\gamma_{\perp}\f$ for a type
/// @param [in] typeindex Index of the type
/// @return Friction coefficient in the orthogonal direction
inline double ReferenceMap::getGammaOrthogonal(uint8_t typeindex) {
  return gammaOrthoFast[typeindex];
}


/// @brief Get exothermicity for the particle
/// @param [in] typeindex Index of the type
/// @return Exothermicity of the associated chemical reaction
inline double ReferenceMap::getExothermicity(uint8_t typeindex) {
  return exothermicityFast[typeindex];
}


/// @brief Get velocity of the particle (for virtual particles only)
/// @param [in] typeindex Index of the type
/// @return Constant velocity of the particle
inline vec3<double> ReferenceMap::getVelocity(uint8_t typeindex) {
  return velocityFast[typeindex];
}


/// @brief Get the particle class from the name of it's type
/// @param [in] name Particle type name
/// @return Class of the particle (enum type)
inline ReferenceMap::ParticleType ReferenceMap::findType(std::string name) {
  auto it = indexReference.find(name);
  if (it!=indexReference.end()) return findType(it->second);
  else return NOT_FOUND;
}


/// @brief Get a subtype from it's index
/// @param [in] id Index of the subtype
/// @return Subtype
inline std::string ReferenceMap::findSubtype(uint8_t id) {
  return subtypesReference[id];
}


/// @brief Get a pointer to the particle type from it's name
/// @param [in] name Particle type name
/// @return Pointer to the type
inline TypeParticle* ReferenceMap::find(std::string name) {
  return find(indexReference.find(name)->second);
}


/// @brief Get a type index given a type name
/// @param [in] name Particle type name
/// @return Index of the type
inline uint8_t ReferenceMap::findTypeIndex(std::string name) {
  return indexReference.find(name)->second;
}


/// @brief Get address of a potential given two particle types
/// @param [in] typeindex First type
/// @param [in] typejndex Second type
/// @return Potential
inline IPotential* ReferenceMap::getPotential(uint8_t typeindex, uint8_t typejndex) {
  return potentialReference.find(std::make_pair(auxMin(typeindex, typejndex), auxMax(typeindex, typejndex)))->second;
}


/// @brief Get address of an eos given the type
inline IEOS* ReferenceMap::getEOS(uint8_t typeindex) {
  return eosReference.find(typeindex)->second;
}


/// @brief Get address of a kernel given the type
inline KernelFunction::Type ReferenceMap::getKernelType(uint8_t typeindex) {
  return kernelType[typeindex];
}


/// @brief Get address of a kernel given the type
inline KernelFunction* ReferenceMap::getKernel(uint8_t typeindex) {
  return kernelReference.find(typeindex)->second;
}



/// @brief Get address of a reaction given the type
inline IKinetics* ReferenceMap::getReaction(uint8_t typeindex) {
  return reactionReference.find(typeindex)->second;
}



/// @brief Get total number of types in the simulation
/// @return Number of types
inline uint8_t ReferenceMap::getNumberOfTypes() {
  return indexReference.size();
}


/// @brief Test if the simulation has (at least) one EAM potential
/// @return Boolean result
inline bool ReferenceMap::isEAM() {
  return isPotentialEAM;
}


/// @brief Test if the simulation has (at least) one MEAM potential
/// @return Boolean result
inline bool ReferenceMap::isMEAM() {
  return isPotentialMEAM;
}


/// @brief Test if the simulation has (at least) one mesoparticle
/// @return Boolean result
inline bool ReferenceMap::isDPD() {
  return existsDPD;
}


/// @brief Test if the simulation has (at least) one smoothparticle
/// @return Boolean result
inline bool ReferenceMap::isSDPD() {
  return existsSDPD;
}


/// @brief Test if the simulation has (at least) one reaction
/// @return Boolean result
inline bool ReferenceMap::isReactive() {
  return existsReactive;
}

/// @brief Test if the simulation included timers
/// @return Boolean result
inline bool ReferenceMap::isTimers() {
  return isLocalTimers;
}

/// @brief Test if the simulation included timers
/// @return Boolean result
inline bool ReferenceMap::isVerlet() {
  return isVERLET;
}


/// @brief Test if the simulation included timers
/// @return Boolean result
inline bool ReferenceMap::isBlockVerlet() {
  return isBLOCKVERLET && !isPotentialMEAM; // disable with a MEAM potential

}

/// @brief Test if the simulation included timers
/// @return Boolean result
inline void ReferenceMap::setIsBlockVerlet() {
  isBLOCKVERLET=true;
}

/// @brief Test if the type is part of a wall
inline bool ReferenceMap::isWallType(uint8_t typeindex) {
  return findType(typeindex)==ParticleType::WALL_PARTICLE;
}




#endif // __REFERENCE_MAP_HPP_INCLUDED

