/// @file 
/// @brief Implementation of class ReferenceMap


#include <iostream>

#include "globals.hpp"
#include "referenceMap.hpp"

#include "cellList/neighborList.hpp"

#include "eos/allEOS.hpp"

#include "potential/allPotentials.hpp"

#include "parallel/types/MPI_inMol.hpp"

#include "particle/types/allTypes.hpp"
#include "particle/molecule/configInMol.hpp"

#include "utils/kernelFunction/allKernelFunctions.hpp"


/// @brief Print some warnings if problem during type registration
///
///
void printRegisterTypeWarning() {
  std::cout<< "WARNING: in ReferenceMap::Register() : type to register is "
	   << "already here! " <<std::endl;
}


/// @brief Print some warnings if problem during type unregistration
///
///
void printUnregisterTypeWarning() {
  std::cout<< "WARNING: in ReferenceMap::Register() : type to unregister " 
	   << "not found! " <<std::endl;
}


/// @brief Default constructor
///
///
ReferenceMap::ReferenceMap()
  : nextIndex(0), indexReference(),
    atomReference(),mesoReference(),smoothReference(),
    potentialReference(), 
    eosReference(),kernelReference(),
    dpdReference(), dpdeReference(), sdpdReference(),
    reactionReference(),
    ghostThickness(), isPotentialEAM(false), isPotentialMEAM(false), existsDPD(false), existsSDPD(false), isLocalTimers(false), existsReactive(false),
    isVERLET(false), isBLOCKVERLET(false),
    rcutMax(-1.), rcutFast(nullptr), rcutFast2(nullptr), costFast(nullptr), massFast(nullptr), massInvFast(nullptr),
    gammaParaFast(nullptr), gammaOrthoFast(nullptr), exothermicityFast(nullptr), velocityFast(nullptr), kernelType(nullptr), Dmax(0), r_amrCriterion(0) {
}


/// @brief Destructor
///
///
ReferenceMap::~ReferenceMap() {

  // Delete all potentials
  for (auto it=potentialReference.begin(); it!=potentialReference.end(); ++it) {
    delete it->second;
  }

  // Delete all eos
  for (auto it=eosReference.begin(); it!=eosReference.end(); ++it) {
    delete it->second;
  }

  // Delete all kernels
  for (auto it=kernelReference.begin(); it!=kernelReference.end(); ++it) {
    delete it->second;
  }

  // Delete all reactions
  for (auto it=reactionReference.begin(); it!=reactionReference.end(); ++it) {
    delete it->second;
  }

  // Erase all maps
  indexReference.clear();
  atomReference.clear();
  mesoReference.clear();
  smoothReference.clear();
  wallReference.clear();
  
  dpdReference.clear();
  dpdeReference.clear();
  sdpdReference.clear();
  sdpdWallReference.clear();

  reactionReference.clear();
  
  // Reset nextIndex
  nextIndex = 0;

  // 
  if (rcutFast !=nullptr) delete [] rcutFast;
  if (rcutFast2!=nullptr) delete [] rcutFast2;
  if (costFast !=nullptr) delete [] costFast;
  if (massFast !=nullptr) delete [] massFast;
  if (massInvFast !=nullptr) delete [] massInvFast;
  if (gammaParaFast != nullptr) delete [] gammaParaFast;
  if (gammaOrthoFast != nullptr) delete [] gammaOrthoFast;
  if (exothermicityFast != nullptr) delete [] exothermicityFast;
  if (velocityFast != nullptr) delete [] velocityFast;
  if (kernelType != nullptr) delete [] kernelType;

}


/// @brief Create and reference all types
/// @param [in] configuration Input configuration for the particles
void ReferenceMap::configure(const Configuration<TypeParticle>& configuration) {

  // Create temporary TypeAtom and register them
  int numAtoms = configuration.atomNames.size();
  for (int i=0; i<numAtoms; ++i) {
    TypeAtom tmp( configuration.atomNames[i], 
		  configuration.atomMasses[i], 
		  configuration.atomAtomicNumbers[i] );
    Register(tmp);
  }

  // Create temporary TypeMesoparticle and register them
  int numMesos = configuration.mesoNames.size();
  existsDPD |= numMesos>0;
  for (int i=0; i<numMesos; ++i) {
    TypeMesoparticle tmp( configuration.mesoNames[i], 
			  configuration.mesoMasses[i], 
			  configuration.mesoUnitMasses[i],
			  configuration.mesoGammaPara[i],
			  configuration.mesoGammaOrtho[i]);
    Register(tmp);
  }

  // Create temporary TypeSmoothParticle and register them
  int numSmooths = configuration.smoothNames.size();
  existsDPD |= numSmooths>0;
  existsSDPD |= numSmooths>0;
  for (int i=0; i<numSmooths; ++i) {
    TypeSmoothParticle tmp( configuration.smoothNames[i], 
			    configuration.smoothMasses[i], 
			    configuration.smoothUnitMasses[i],
			    configuration.smoothBulkViscosity[i],
			    configuration.smoothShearViscosity[i],
			    configuration.smoothSmoothingLength[i]);
    Register(tmp);

    // Register corresponding kernel function
    if(configuration.smoothKernel[i] == "lucy") {
      KernelFunction* tmpk = new Lucy( tmp.getSmoothingLength() );
      Register(tmpk,configuration.smoothNames[i]);
    }
    else { 
      std::cout << "WARNING: Unknown kernel function " << configuration.smoothKernel[i] << std::endl;
      std::cout << "Lucy kernel chosen by default" << std::endl;
      KernelFunction* tmpk = new Lucy( tmp.getSmoothingLength() );
      Register(tmpk,configuration.smoothNames[i]);
    }
  }

  // For wall particles
  // Create temporary TypeSmoothWallParticle and register them
  int numWalls = configuration.wallNames.size();
  for (int i=0; i<numWalls; ++i) {
    TypeSmoothWallParticle tmp( configuration.wallNames[i], 
				configuration.wallMasses[i], 
				configuration.wallUnitMasses[i],
				configuration.wallSmoothingLength[i],
				configuration.wallVelocity[i]);
    Register(tmp);

    // Register corresponding kernel function
    if(configuration.wallKernel[i] == "lucy") {
      KernelFunction* tmpk = new Lucy( tmp.getSmoothingLength() );
      Register(tmpk,configuration.wallNames[i]);
    }
    else { 
      std::cout << "WARNING: Unknown kernel function " << configuration.wallKernel[i] << std::endl;
      std::cout << "Lucy kernel chosen by default" << std::endl;
      KernelFunction* tmpk = new Lucy( tmp.getSmoothingLength() );
      Register(tmpk,configuration.wallNames[i]);
    }
  }

  // Later : other types
  // ...

  // Initialize the variables that depend on the type
  uint numTypes = indexReference.size();
  NeighborList_base::setNumberOfTypes(numTypes);
  rcutFast          = new double[numTypes*numTypes];
  costFast          = new double[numTypes*numTypes];
  rcutFast2         = new double[numTypes*numTypes];
  
  // Used to store the max between rcut and rverlet.
  // If the linked cell method is used, the two following array are equal to rcut.
  rcutFastVerlet    = new double[numTypes*numTypes];
  rcutFast2Verlet   = new double[numTypes*numTypes];

  massFast          = new double[numTypes];
  massInvFast       = new double[numTypes];
  gammaParaFast     = new double[numTypes];
  gammaOrthoFast    = new double[numTypes];
  exothermicityFast = new double[numTypes];
  velocityFast      = new vec3<double>[numTypes];
  kernelType        = new KernelFunction::Type[numTypes];

  // Set the masses and mass inverses
  for (uint t=0; t<numTypes; ++t) {
    massFast[t] = find(t)->getMass();
    massInvFast[t] = find(t)->getInvMass();
    gammaParaFast[t] = find(t)->getGammaParallel();
    gammaOrthoFast[t] = find(t)->getGammaOrthogonal();
    exothermicityFast[t] = 0;
    velocityFast[t] = find(t)->getVelocity();
    kernelType[t] = KernelFunction::NONE;

    
    for (uint tt=0; tt<numTypes; ++tt) {
      costFast[t*numTypes+tt] = 0;
      rcutFast[t] = -1;
      rcutFast2[t] = -1;
      rcutFastVerlet[t] = -1;
      rcutFast2Verlet[t] = -1;
    }
  }

  for (auto it=kernelReference.begin(); it!=kernelReference.end(); ++it)
    kernelType[it->first] = it->second->getType();
}


/// @brief Create all potentials and link them to types
/// @param [in] configuration Input configuration for the potentials
void ReferenceMap::configure(const Configuration<IPotential>& configuration) {

  // Ideal Gases
  int numIG = configuration.IGtypeA.size();
  for (int i=0; i<numIG; ++i) {
    IPotential* tmpPot = new IdealGas();
    Register(tmpPot, configuration.IGtypeA[i], configuration.IGtypeB[i]);
  }

  // Lennard-Jones potentials
  int numLJ = configuration.LJtypeA.size();
  for (int i=0; i<numLJ; ++i) {
    IPotential* tmpPot = new LennardJonesPotential(configuration.LJrcut[i], configuration.LJparameters[i]);
    Register(tmpPot, configuration.LJtypeA[i], configuration.LJtypeB[i]);
  }

  // // Exp6 potentials
  // int numExp6 = configuration.Exp6typeA.size();
  // for (int i=0; i<numExp6; ++i) {
  //   IPotential* tmpPot = new Exp6Potential(configuration.Exp6rcut[i], configuration.Exp6parameters[i]);
  //   Register(tmpPot, configuration.Exp6typeA[i], configuration.Exp6typeB[i]);
  // }

  // Sutton-Chen potentials (with or without symmetrization)
  int numSC = configuration.SCtypeA.size();
  isPotentialEAM |= numSC>0;
  for (int i=0; i<numSC; ++i) {
    IPotential* tmpPot = new SuttonChenPotential(configuration.SCrcut[i], configuration.SCparameters[i]);
    Register(tmpPot, configuration.SCtypeA[i], configuration.SCtypeB[i]);
  }

  // Eam Vniitf potentials (with or without symmetrization)
  int numEamVniitf = configuration.EamVNIITFtypeA.size();
  isPotentialEAM |= numEamVniitf>0;
  for (int i=0; i<numEamVniitf; ++i) {
    IPotential* tmpPot = new EamVniitfPotential(configuration.EamVNIITFrcut[i], configuration.EamVNIITFparameters[i]);
    Register(tmpPot, configuration.EamVNIITFtypeA[i], configuration.EamVNIITFtypeB[i]);
  }

 // meam potentials 
  int numMeam = configuration.MeamTypeA.size();
  isPotentialMEAM |= numMeam>0;
  for (int i=0; i<numMeam; ++i) {
    if(!configuration.symmetrize){
   	 IPotential* tmpPot;
	 tmpPot= new MeamPotential(configuration.MeamRcut[i], configuration.MeamParam[i]);
         Register(tmpPot, configuration.MeamTypeA[i], configuration.MeamTypeB[i]);
	}
	else{
	std::cout<<std::endl<<"                ****************************************************"<<std::endl<<"                ** WARNING: choice parameter symetrisation = true **"<<std::endl<<"                ** It's impossible to use symetrisation with MEAM **"<<std::endl<<"                ** Try with symetrisation = false                 **"<<std::endl<<"                ****************************************************"<<std::endl<<std::endl ;
	exit(0);
	}

  }

  // Gaussian potentials
  int numGauss = configuration.GaussTypeA.size();
  for (int i=0; i<numGauss; ++i) {
    IPotential* tmpPot = new GaussianPotential(configuration.GaussRcut[i], configuration.GaussParameters[i]);
    Register(tmpPot, configuration.GaussTypeA[i], configuration.GaussTypeB[i]);
  }



  // For long range potentials set rcut to rcutMax
  for (unsigned int i=0; i<auxSq(indexReference.size()); ++i) {
    if (rcutFast[i] < 0.) {
      rcutFast [i] = rcutMax;
      rcutFast2[i] = auxSq(rcutMax);
    }
  }
  

  if(isVerlet()){
  
    rcutMax+=rVerlet;
    for (unsigned int i=0; i<auxSq(indexReference.size()); ++i) {
      if (rcutFastVerlet[i] < 0.) {
        rcutFastVerlet [i] = rcutMax;
        rcutFast2Verlet[i] = auxSq(rcutMax);
      }
    }
  }


  // Set ghost thickness
  setGhostThickness();

}


/// @brief Build a type for each atom and register one of each type
void ReferenceMap::configure(Configuration<MPI__InMol>& configuration, const Configuration<IPotential>& potentials) {
  //
  // Create and register atom subtypes and types, then build a temporary storage of the type for each subtype
  //
  std::vector<uint8_t> theTypes; // Temporary storage of the type index for each subtype
  // For each subtype from the configuration list
  for(uint iTy(0); iTy<configuration.m_listOfSubtypes.size(); ++iTy) {
    std::string subtype,data; // Data from the reader and real subtype
    int atomicNumber; // Atomic number
    double mass; // Mass
    // Get the data
    data=configuration.m_listOfSubtypes[iTy];
    // Extract the subtype
    subtype=data.substr(0,data.find('$'));
    data=data.substr(data.find('$')+1);
    // Extract the atomic number
    atomicNumber=std::stoi(data.substr(0,data.find('$')));
    data=data.substr(data.find('$')+1);
    // Extract the mass
    mass=std::stod(data);
    // Register the subtype
    subtypesReference[iTy]=subtype;
    // Get the type
    std::string type=Global::ffield->getAtomType(subtype);
    // Make and register an atom type from the data
    TypeAtom typeAtom(type, mass, atomicNumber);
    RegisterNoWarn(typeAtom);
    // Put the type index into the temporary type indexes storage
    theTypes.push_back(indexReference[type]);
  }
  configuration.m_listOfSubtypes.clear();
  //
  // Set the atom types
  //
  for(uint iAt(0); iAt<configuration.m_particles.size(); ++iAt) {
    // Get the index corresponding to the atom type
    MPI__InMol* ptrInMol=dynamic_cast<MPI__InMol*>(configuration.m_particles[iAt]);
    configuration.m_particles[iAt]->type()=theTypes[ptrInMol->subtype()];
  }
  //
  // Things to do once the number of types is known
  //
  uint numTypes = indexReference.size();
  NeighborList_base::setNumberOfTypes(numTypes); // Pass the number of type to the neighbor list
  // Initialize the arrays whose size depend on the number of types
  rcutFast  = new double[numTypes*numTypes];
  costFast  = new double[numTypes*numTypes];
  rcutFast2 = new double[numTypes*numTypes];
  massFast  = new double[numTypes];
  massInvFast  = new double[numTypes];
  gammaParaFast  = new double[numTypes];
  gammaOrthoFast = new double[numTypes];
  velocityFast   = new vec3<double>[numTypes];
  kernelType     = new KernelFunction::Type[numTypes];
  
  if(isVerlet()){
  
    std::cout << " Verlet Lists are not implemented for molecules " << std::endl;
    exit(0);
     
  }

  //
  // Set the masses and mass inverses
  //
  for (uint t=0; t<numTypes; ++t) {
    massFast[t] = find(t)->getMass();
    massInvFast[t] = find(t)->getInvMass();
    gammaParaFast[t] = find(t)->getGammaParallel();
    gammaOrthoFast[t] = find(t)->getGammaOrthogonal();
    velocityFast[t] = find(t)->getVelocity();
    kernelType[t] = KernelFunction::NONE;
    rcutFast[t] = -1;
  }
  //
  // Create and register the potentials
  //
  // Force field potentials
  for(auto it1=indexReference.begin(); it1!=indexReference.end(); ++it1) {
    for(auto it2=indexReference.begin(); it2!=indexReference.end(); ++it2) {
  		if(it1->first.compare(it2->first)>=0) {
    		IPotential* tmpPot=Global::ffield->makePotential(it1->first,it2->first);
    		Register(tmpPot, it1->first, it2->first);
  		}
  	}
  }
  // Overloaded potentials
  // Ideal Gases
  int numIG = potentials.IGtypeA.size();
  for (int i=0; i<numIG; ++i) {
    IPotential* tmpPot = new IdealGas();
    Register(tmpPot, potentials.IGtypeA[i], potentials.IGtypeB[i]);
  }
  // Lennard-Jones potentials
  int numLJ = potentials.LJtypeA.size();
  for (int i=0; i<numLJ; ++i) {
    IPotential* tmpPot = new LennardJonesPotential(potentials.LJrcut[i], potentials.LJparameters[i]);
    Register(tmpPot, potentials.LJtypeA[i], potentials.LJtypeB[i]);
  }
  // Sutton-Chen potentials (with or without symmetrization)
  int numSC = potentials.SCtypeA.size();
  isPotentialEAM |= numSC>0;
  for (int i=0; i<numSC; ++i) {
    IPotential* tmpPot = new SuttonChenPotential(potentials.SCrcut[i], potentials.SCparameters[i]);
    Register(tmpPot, potentials.SCtypeA[i], potentials.SCtypeB[i]);
  }
  // Eam Vniitf potentials (with or without symmetrization)
  int numEamVniitf = potentials.EamVNIITFtypeA.size();
  isPotentialEAM |= numEamVniitf>0;
  for (int i=0; i<numEamVniitf; ++i) {
    IPotential* tmpPot = new EamVniitfPotential(potentials.EamVNIITFrcut[i], potentials.EamVNIITFparameters[i]);
    Register(tmpPot, potentials.EamVNIITFtypeA[i], potentials.EamVNIITFtypeB[i]);
  }
  // Consider the molecules in rcutMax
  rcutMax=auxMax(rcutMax,ForceField::s_rcut);
  // For long range potentials set rcut to rcutMax
  for (unsigned int i=0; i<auxSq(indexReference.size()); ++i) {
    if (rcutFast[i] < 0.) {
      rcutFast [i] = rcutMax;
      rcutFast2[i] = auxSq(rcutMax);
    }
  }
  //
  // Set ghost thickness
  //
  setGhostThickness();

}


/// @brief Create all eos and link them to types
/// @param [in] configuration Configuration data for EOS
void ReferenceMap::configure(const Configuration<IEOS>& configuration) {

  // Check whether there is a reactive EOS and store the corresponding type
  int numReactive = configuration.Reactivetype.size();
  Array< std::pair<IEOS*,IEOS*> > storeReactiveEOS(numReactive,std::pair<IEOS*,IEOS*>(nullptr,nullptr));
  
  
  // DPDE
  int numDPDE = configuration.DPDEtype.size();
  for (int i=0; i<numDPDE; ++i) {
    IEOS* tmpEos = new DPDEEOS(configuration.DPDEcv[i]);
    if (!checkReactiveEOS(storeReactiveEOS,tmpEos,configuration.ReactiveEOS0,configuration.ReactiveEOS1,
			  configuration.DPDEtype[i],configuration.Reactivetype))
      Register(tmpEos, configuration.DPDEtype[i]);
  }

  // Ideal Gas
  int numIG = configuration.IGtype.size();
  for (int i=0; i<numIG; ++i) {
    TypeMesoparticle* ptype = static_cast<TypeMesoparticle*>(find(configuration.IGtype[i]));
    IEOS* tmpEos = new IdealGasEOS(ptype->getMass()/ptype->getUnitMass(),ptype->getMass());
    if (!checkReactiveEOS(storeReactiveEOS,tmpEos,configuration.ReactiveEOS0,configuration.ReactiveEOS1,
			  configuration.IGtype[i],configuration.Reactivetype))
      Register(tmpEos, configuration.IGtype[i]);
  }

  // Mie-Gruneisen
  int numMG = configuration.MGtype.size();
  for(int i=0; i<numMG; ++i) {
    TypeMesoparticle* ptype = static_cast<TypeMesoparticle*>(find(configuration.MGtype[i]));
    IEOS* tmpEos = new MieGruneisenEOS(ptype->getMass()/ptype->getUnitMass(),ptype->getMass(),
				       configuration.MGgamma0[i],configuration.MGgammaInf[i],configuration.MGtheta0[i],
				       configuration.MGq[i],configuration.MGrho0[i],configuration.MGks[i],configuration.MGns[i],
				       configuration.MGrhos[i],configuration.MGur[i],configuration.MGcvr[i]);
    if (!checkReactiveEOS(storeReactiveEOS,tmpEos,configuration.ReactiveEOS0,configuration.ReactiveEOS1,
			  configuration.MGtype[i],configuration.Reactivetype))
      Register(tmpEos, configuration.MGtype[i]);
  }

  // HZ
  int numHZ = configuration.HZtype.size();
  for(int i=0; i<numHZ; ++i) {
    TypeMesoparticle* ptype = static_cast<TypeMesoparticle*>(find(configuration.HZtype[i]));
    IEOS* tmpEos = new HZEOS(ptype->getMass()/ptype->getUnitMass(),ptype->getMass(),
			     configuration.HZgamma0[i],configuration.HZrho0[i],configuration.HZc0[i],
			     configuration.HZcv[i],configuration.HZs[i]);
    if (!checkReactiveEOS(storeReactiveEOS,tmpEos,configuration.ReactiveEOS0,configuration.ReactiveEOS1,
			  configuration.HZtype[i],configuration.Reactivetype))
      Register(tmpEos, configuration.HZtype[i]);
  }

  // JWL
  int numJWL = configuration.JWLtype.size();
  for(int i=0; i<numJWL; ++i) {
    TypeMesoparticle* ptype = static_cast<TypeMesoparticle*>(find(configuration.JWLtype[i]));
    IEOS* tmpEos = new JWLEOS(ptype->getMass()/ptype->getUnitMass(),ptype->getMass(),
			      configuration.JWLgamma0[i],configuration.JWLrho0[i],configuration.JWLe0[i],
			      configuration.JWLdcj[i],configuration.JWLpcj[i],configuration.JWLtcj[i],configuration.JWLcv[i],
			      configuration.JWLa[i],configuration.JWLb[i],configuration.JWLr1[i],configuration.JWLr2[i]);
    if (!checkReactiveEOS(storeReactiveEOS,tmpEos,configuration.ReactiveEOS0,configuration.ReactiveEOS1,
			  configuration.JWLtype[i],configuration.Reactivetype))
      Register(tmpEos, configuration.JWLtype[i]);
  }

  // Reactive
  for(int i=0; i<numReactive; ++i) {
    TypeMesoparticle* ptype = static_cast<TypeMesoparticle*>(find(configuration.Reactivetype[i]));

    if (storeReactiveEOS[i].first == nullptr || storeReactiveEOS[i].second == nullptr) {
      std::cerr << "Error: At least one equation of state is missing for reactive EOS" << std::endl;
      exit(0);
    }

    IEOS* tmpEos = createReactiveEOS(storeReactiveEOS[i].first->getType(),storeReactiveEOS[i].second->getType(),ptype->getSize(),ptype->getMass(),storeReactiveEOS[i].first,storeReactiveEOS[i].second);
    Register(tmpEos, configuration.Reactivetype[i]);
  }

  // Check if there is more than one EOS per type !
  // TODO : Check where to delete EOS pointers (are they even deleted now ?)
  
}




/// @brief Create all reaction kinetics and link them to types
/// @param [in] configuration Configuration data for chemical reactions
void ReferenceMap::configure(const Configuration<IKinetics>& configuration) {
  
  // Second order
  int numSO = configuration.SOtype.size();
  existsReactive |= numSO > 0;
  for (int i=0; i<numSO; ++i) {
    IKinetics* tmpReac = new SecondOrder(configuration.SOzab[i],configuration.SOzba[i],configuration.SOeab[i],configuration.SOeba[i]);
    Register(tmpReac, configuration.SOtype[i]);
  }

}




/// @brief Create all (s)dpd(e) interactions and link them to types
/// @param [in] input Input data
void ReferenceMap::configure(const Input& input) {

  // AMR
  isLocalTimers=input.isTimers; //Rappel trouver le bon endroit
  Dmax = input.Dmax;
  r_amrCriterion = input.amrCriterion;

  // DPD
  int numDPD = input.DPDtypeA.size();
  for (int i=0; i<numDPD; ++i) {
    Register(DPD, input.DPDtypeA[i], input.DPDtypeB[i]);
  }

  // DPD
  int numDPDE = input.DPDEtypeA.size();
  for (int i=0; i<numDPDE; ++i) {
    Register(DPDE, input.DPDEtypeA[i], input.DPDEtypeB[i]);
  }

  // SDPD
  int numSDPD = input.SDPDtypeA.size();
  for (int i=0; i<numSDPD; ++i) {
    Register(SDPD, input.SDPDtypeA[i], input.SDPDtypeB[i]);
  }

  // SDPD
  int numSDPD_WALL = input.SDPDWalltypeReal.size();
  for (int i=0; i<numSDPD_WALL; ++i) {
    Register(SDPD_WALL, input.SDPDWalltypeReal[i], input.SDPDWalltypeVirtual[i]);
  }

  for (auto it=kernelReference.begin(); it!=kernelReference.end(); ++it)
    kernelType[it->first] = it->second->getType();

}



/// @brief Set ghost thickness
///
///
void ReferenceMap::setGhostThickness() {

  ghostThickness = 0;

  // Loop on all potentials, get the max
  for (auto it=potentialReference.begin(); it!=potentialReference.end(); ++it) {

    ShortRangePotential* potential = dynamic_cast<ShortRangePotential*>(it->second);

    if (potential!=nullptr) 
      ghostThickness = auxMax(ghostThickness, potential->getGhostThickness());

  }

  // if there is SDPD, ensure that ghostThickness >= 1
  if(isSDPD())
    ghostThickness = auxMax(ghostThickness, uint(1));

}


/// @brief Get the particle class from the index of it's type
/// @param [in] typeIndex Particle type index
/// @return Class of the particle (enum type)
ReferenceMap::ParticleType ReferenceMap::findType(uint8_t typeIndex) {
  
  auto it1 = atomReference.find(typeIndex);
  if (it1!=atomReference.end()) return ATOM;

  auto it2 = mesoReference.find(typeIndex);
  if (it2!=mesoReference.end()) return MESOPARTICLE;

  auto it3 = smoothReference.find(typeIndex);
  if (it3!=smoothReference.end()) return SMOOTH_PARTICLE;

  auto it4 = wallReference.find(typeIndex);
  if (it4!=wallReference.end()) return WALL_PARTICLE;

  return NOT_FOUND;
}


/// @brief Get the adress of a type given its index
/// @param [in] typeIndex Particle type index
/// @return Pointer to the type
TypeParticle* ReferenceMap::find(uint8_t typeIndex)  {
  
  TypeAtomMap::iterator it1 = atomReference.find(typeIndex);
  if (it1!=atomReference.end()) return &(it1->second);

  TypeMesoMap::iterator it2 = mesoReference.find(typeIndex);
  if (it2!=mesoReference.end()) return &(it2->second);
  
  TypeSmoothMap::iterator it3 = smoothReference.find(typeIndex);
  if (it3!=smoothReference.end()) return &(it3->second);

  TypeWallMap::iterator it4 = wallReference.find(typeIndex);
  if (it4!=wallReference.end()) return &(it4->second);

  return nullptr;
}


/// @brief Add an atom type to the maps
/// @param [in] type Atom type
void ReferenceMap::Register(TypeAtom& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==NOT_FOUND) {
    indexReference.insert(std::make_pair(type.getName(), nextIndex++));
    atomReference.insert(std::make_pair(indexReference[type.getName()], type));
  }
  else printRegisterTypeWarning();
}


/// @brief Add an atom type to the maps, no warning if the type is already there
/// @param [in] type Atom type
void ReferenceMap::RegisterNoWarn(TypeAtom& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==NOT_FOUND) {
  	indexReference.insert(std::make_pair(type.getName(), nextIndex++));
  	atomReference.insert(std::make_pair(indexReference[type.getName()], type));
  }
}


/// @brief Remove an atom type to the maps (not used)
/// @param [in] type Atom type
void ReferenceMap::UnRegister(TypeAtom& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==ATOM) {
    atomReference.erase( indexReference[type.getName()] );
    indexReference.erase(type.getName());
  }
  else printUnregisterTypeWarning();
}


/// @brief  Register a type a typeMesoparticle
/// @param [in] type Mesoparticle type
void ReferenceMap::Register(TypeMesoparticle& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==NOT_FOUND) {
    indexReference.insert(std::make_pair(type.getName(), nextIndex++));
    mesoReference.insert(std::make_pair(indexReference[type.getName()], type));
  }
  else printRegisterTypeWarning();
}


/// @brief Unregister a typeMesoparticle (not used)
/// @param [in] type Mesoparticle type
void ReferenceMap::UnRegister(TypeMesoparticle& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==MESOPARTICLE) {
    mesoReference.erase( indexReference[type.getName()] );
    indexReference.erase(type.getName());
  }
  else printUnregisterTypeWarning();
}


/// @brief Register a type a typeSmoothParticle 
/// @param [in] type Smooth particle type
void ReferenceMap::Register(TypeSmoothParticle& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==NOT_FOUND) {
    indexReference.insert(std::make_pair(type.getName(), nextIndex++));
    smoothReference.insert(std::make_pair(indexReference[type.getName()], type));
  }
  else printRegisterTypeWarning();
}


/// @brief Unregister a typeSmoothParticle (not used)
/// @param [in] type Smooth particle type
void ReferenceMap::UnRegister(TypeSmoothParticle& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==SMOOTH_PARTICLE) {
    smoothReference.erase( indexReference[type.getName()] );
    indexReference.erase(type.getName());
  }
  else printUnregisterTypeWarning();
}


/// @brief Register a type a typeSmoothWallParticle
/// @param [in] type Smooth particle type
void ReferenceMap::Register(TypeSmoothWallParticle& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==NOT_FOUND) {
    indexReference.insert(std::make_pair(type.getName(), nextIndex++));
    wallReference.insert(std::make_pair(indexReference[type.getName()], type));
  }
  else printRegisterTypeWarning();
}


/// @brief Unregister a typeSmoothWallParticle (not used)
/// @param [in] type Smooth particle type
void ReferenceMap::UnRegister(TypeSmoothWallParticle& type) {
  ParticleType particleType = findType(type.getName());
  if (particleType==WALL_PARTICLE) {
    smoothReference.erase( indexReference[type.getName()] );
    indexReference.erase(type.getName());
  }
  else printUnregisterTypeWarning();
}


/// @brief Add a potential linked to two particle types
///
/// If there is already a potential for those types, this will replace it
/// @param [in] potential Potential to add
/// @param [in] nameA Name of the first particle type
/// @param [in] nameB Name of the second particle type
void ReferenceMap::Register(IPotential* potential, std::string nameA, 
			 std::string nameB) {

  // Find the particle types from the names
  ParticleType particleTypeA = findType(nameA);
  ParticleType particleTypeB = findType(nameB);

  // Check if the types match in the database
  if (particleTypeA==NOT_FOUND || particleTypeB==NOT_FOUND) {
    std::cout<< "WARNING: in ReferenceMap::Register() : one of the types ("
	     << nameA <<" and " << nameB <<") not found (at least)! " << std::endl;
    return;
  }
  // This seem to be a wrong assumption
  //	  else if (particleTypeA!=particleTypeB) {
  // 	   std::cout<< "WARNING: in ReferenceMap::Register() : a potential is an "
  //		     << "interaction between two particles of the same type! "
  //	  	   << std::endl;
  //   	 return;
  // 	 }

  // Get the index
  uint8_t A = indexReference[nameA], B = indexReference[nameB];

  // Put the potential in the map
  if (A<B)
    potentialReference[std::make_pair(A,B)]=potential;
  else
    potentialReference[std::make_pair(B,A)]=potential;

  // If it is a short-range potential, update the max cutoff radius
  ShortRangePotential* shortpot = dynamic_cast<ShortRangePotential*>(potential);
  double rcut;

  if (shortpot!=nullptr) {
    rcut = shortpot->getCutoffRadius();
    rcutMax = auxMax(rcutMax, rcut);
  }
  else {
    rcut = -1.;
  }

  // Fill potential dependent data
  rcutFast [A*nextIndex+B] = rcut;
  rcutFast [B*nextIndex+A] = rcut;
  rcutFast2[A*nextIndex+B] = auxSq(rcut);
  rcutFast2[B*nextIndex+A] = auxSq(rcut);



  costFast[A*nextIndex+B] = potential->cost()*auxPow(rcut,3); // Include dependence to rcut
  costFast[B*nextIndex+A] = potential->cost()*auxPow(rcut,3);


  // This part is added to take into account the list of verlet method.
  // To do so, we add a storage for the Verlet raduis.
  // If this method is not used, the verlet radius is equal to the cut-off radius.
  rVerlet = rVerlet < 0 ? 0 : rVerlet;
  
  rcut+=rVerlet;
  rcutFastVerlet [A*nextIndex+B] = rcut;
  rcutFastVerlet [B*nextIndex+A] = rcut;
  rcutFast2Verlet[A*nextIndex+B] = auxSq(rcut);
  rcutFast2Verlet[B*nextIndex+A] = auxSq(rcut);


}


/// @brief Add an equation of state linked to a particle type
/// @param [in] eos Equation of state to add
/// @param [in] name Related particle type
void ReferenceMap::Register(IEOS* eos, std::string name) {

  ParticleType particleType = findType(name);

  // Check if the types match in the database

  if (particleType==NOT_FOUND) {
    std::cout<< "WARNING: in ReferenceMap::Register() : type not "
	     << "found ! " << std::endl;
    return;
  }

  uint8_t A = indexReference[name];

  // Put the eos in the map

  eosReference.insert( std::make_pair(A, eos) );

}


/// @brief add a smoothing kernel linked to a particle type
/// @param [in] kernel Kernel function to add
/// @param[in] name Related particle type
void ReferenceMap::Register(KernelFunction* kernel, std::string name) {

  ParticleType particleType = findType(name);

  // Check if the types match in the database

  if (particleType==NOT_FOUND) {
    std::cout<< "WARNING: in ReferenceMap::Register() : type not "
	     << "found ! " << std::endl;
    return;
  }

  uint8_t A = indexReference[name];

  // Put the kernel in the map

  kernelReference.insert( std::make_pair(A, kernel) );

}


/// @brief  Register a FD interaction
void ReferenceMap::Register(InteractionType itype, std::string nameA, 
			    std::string nameB) {

  ParticleType particleTypeA = findType(nameA);
  ParticleType particleTypeB = findType(nameB);

  // Check if the types match in the database

  if (particleTypeA==NOT_FOUND || particleTypeB==NOT_FOUND) {
    std::cout<< "WARNING: in ReferenceMap::Register() : one of the types not "
	     << "found (at least)! " << std::endl;
    return;
  }

  uint8_t A = indexReference[nameA], B = indexReference[nameB];

  // Put the interaction in the correct map

  switch(itype) {
  case DPD:
    if (A<B)
      dpdReference.push_back( std::make_pair(A, B) );
    else	   
      dpdReference.push_back( std::make_pair(B, A) );

    {
      KernelFunction* tmpk = new DPDKernel( rcutFast[A*nextIndex+B] );
      Register(tmpk, nameA);
    }
    costFast[A*nextIndex+B] += 2;
    costFast[B*nextIndex+A] += 2;
    break;	   

  case DPDE:	   
    if (A<B)	   
      dpdeReference.push_back( std::make_pair(A, B) );
    else	   
      dpdeReference.push_back( std::make_pair(B, A) );
    {
      KernelFunction* tmpk = new DPDKernel( rcutFast[A*nextIndex+B] );
      Register(tmpk, nameA);
    }
    costFast[A*nextIndex+B] += 3;
    costFast[B*nextIndex+A] += 3;

    break;	   
  case SDPD:	   
    if (A<B)
      sdpdReference.push_back( std::make_pair(A, B) );
    else	   
      sdpdReference.push_back( std::make_pair(B, A) );

    rcutFast [A*nextIndex+B] = find(A)->getSmoothingLength();
    rcutFast [B*nextIndex+A] = find(B)->getSmoothingLength();
    rcutFast2[A*nextIndex+B] = rcutFast[A*nextIndex+B]*rcutFast[A*nextIndex+B];
    rcutFast2[B*nextIndex+A] = rcutFast[B*nextIndex+A]*rcutFast[B*nextIndex+A];
    rcutMax = auxMax(rcutMax, auxMax(rcutFast [A*nextIndex+B],rcutFast [B*nextIndex+A]));   

    costFast[A*nextIndex+B] += 6;
    costFast[B*nextIndex+A] += 6;

    break;
  case SDPD_WALL:	   
    sdpdWallReference.push_back( std::make_pair(A, B) );

    rcutFast [A*nextIndex+B] = find(A)->getSmoothingLength();
    rcutFast [B*nextIndex+A] = find(B)->getSmoothingLength();
    rcutFast2[A*nextIndex+B] = rcutFast[A*nextIndex+B]*rcutFast[A*nextIndex+B];
    rcutFast2[B*nextIndex+A] = rcutFast[B*nextIndex+A]*rcutFast[B*nextIndex+A];
    rcutMax = auxMax(rcutMax, auxMax(rcutFast [A*nextIndex+B],rcutFast [B*nextIndex+A]));

    costFast[A*nextIndex+B] += 3;
    costFast[B*nextIndex+A] += 3;

    break;
  default:
    break;
  }
}



/// @brief Add a reaction kinetics linked to a particle type
/// @param [in] reaction Kinetics to add
/// @param [in] name Related particle type
void ReferenceMap::Register(IKinetics* reaction, std::string name) {

  ParticleType particleType = findType(name);

  // Check if the types match in the database

  if (particleType==NOT_FOUND) {
    std::cout<< "WARNING: in ReferenceMap::Register() : type not "
	     << "found ! " << std::endl;
    return;
  }

  uint8_t A = indexReference[name];

  // Put the eos in the map

  reactionReference.insert( std::make_pair(A, reaction) );

  exothermicityFast[A] = reaction->getExothermicity()*find(name)->getSize();
  
}




/// @brief extract data from input file to know which method was chosen
/// @param [in] Verlet currently two possibilities are implemented, Verlet list (true) or linked cell method (false)
/// @param [in] raduis is the rVerlet
void ReferenceMap::neighbours_method(bool Verlet, double raduis) {

    rVerlet=raduis;
    isVERLET=Verlet;
 
  }




/// @brief Check and build a reactive equation of state
/// @param [in,out] store Array of reactive equations of state
/// @param [in] eos Equation of state to be added
/// @param [in] reactiveEOS0Name Array of names of the equation of state for the reactants
/// @param [in] reactiveEOS1Name Array of names of the equation of state for the products
/// @param [in] particleType Type of the particle
/// @param [in] reactiveParticleType Array of particle types involved in a chemical reaction
bool ReferenceMap::checkReactiveEOS(Array< std::pair<IEOS*,IEOS*> >& store, IEOS* eos, const Array<std::string>& reactiveEOS0Name, const Array<std::string>& reactiveEOS1Name, const std::string& particleType, const Array<std::string>& reactiveParticleType) {

  for (uint i = 0; i < reactiveParticleType.size(); i++) {
    if (particleType != reactiveParticleType[i])
      return false;
    if (eos->getName() == reactiveEOS0Name[i] && store[i].first == nullptr) {
      store[i].first = eos;
      return true;
    }
    if (eos->getName() == reactiveEOS1Name[i]) {
      store[i].second = eos;
      return true;
    }
    std::cerr << "Equation of state " << eos->getName() << " does not match the reactive EOS for particle " << particleType << std::endl;
    exit(1);
  }

  return false;
}


/// @brief Create reactive equation of states
/// @param [in] EOS0 Type of first EOS
/// @param [in] EOS1 Type of second EOS
/// @param [in] size Size of a mesoparticle
/// @param [in] mass Mass of a mesoparticle
/// @param [in] eos0 First equation of state
/// @param [in] eos1 Second equation of state
/// @return Pointer to new reactive EOS
IEOS* ReferenceMap::createReactiveEOS(IEOS::Type EOS0, IEOS::Type EOS1,double size, double mass, IEOS* eos0, IEOS* eos1) {

  switch(EOS0) {
  case IEOS::DPDE_ET:
  case IEOS::IDEAL_GAS:
  case IEOS::MIE_GRUNEISEN:
  case IEOS::REACTIVE:
    std::cerr << "Error: Equation of state " << eos0->getName() << " cannot be used in reactive EOS" <<std::endl;
    exit(0);
    break;
    
  case IEOS::HZ:
    return __createReactiveEOS<HZEOS>(EOS1,size,mass,static_cast<HZEOS*>(eos0),eos1);
    break;
    
  case IEOS::JWL:
    return __createReactiveEOS<JWLEOS>(EOS1,size,mass,static_cast<JWLEOS*>(eos0),eos1);
    break;
  }

  return nullptr;
  
}


/// @brief Create reactive equation of states
/// @tparam EOS0 Type of first EOS
/// @param [in] EOS1 Type of second EOS
/// @param [in] size Size of a mesoparticle
/// @param [in] mass Mass of a mesoparticle
/// @param [in] eos0 First equation of state
/// @param [in] eos1 Second equation of state
/// @return Pointer to new reactive EOS
template <class EOS0>
IEOS* ReferenceMap::__createReactiveEOS(IEOS::Type EOS1, double size, double mass, EOS0* eos0, IEOS* eos1) {

  IEOS* tmpEOS = nullptr;
  
  switch(EOS1) {
  case IEOS::DPDE_ET:
  case IEOS::IDEAL_GAS:
  case IEOS::MIE_GRUNEISEN:
  case IEOS::REACTIVE:
    std::cerr << "Error: Equation of state " << eos0->getName() << " cannot be used in reactive EOS" <<std::endl;
    exit(0);
    break;
    
  case IEOS::HZ:
    tmpEOS = new ReactiveEOS<EOS0,HZEOS>(size,mass,eos0,static_cast<HZEOS*>(eos1));
    return tmpEOS;
    break;
    
  case IEOS::JWL:
    tmpEOS = new ReactiveEOS<EOS0,JWLEOS>(size,mass,eos0,static_cast<JWLEOS*>(eos1));
    return tmpEOS;
    break;
  }

  return nullptr;
}




/// @brief  Print the particle types and potentials to check
/// @param [in,out] flux Print flux
void ReferenceMap::print(std::ostream& flux) {

  using namespace std;

  uint maxSize = 0;
  for (auto& it : atomReference) maxSize = auxMax(maxSize, (uint)it.second.getName().size());
  for (auto& it : mesoReference) maxSize = auxMax(maxSize, (uint)it.second.getName().size());

  if (atomReference.size()>0) {
    flux<< "  Atom type(s)           : ";
    for (auto& it : atomReference) flux<< it.second.getName() << " ";
    flux<< endl;
  }

  if (mesoReference.size()>0) {
    flux<< "  Mesoparticle type(s)   : ";
    for (auto& it : mesoReference) flux<< it.second.getName() << " ";
    flux<< endl;
  }

  if (smoothReference.size()>0) {
    flux<< "  Smoothparticle type(s) : ";
    for (auto& it : smoothReference) flux<< it.second.getName() << " ";
    flux<< endl;
  }

  if (wallReference.size()>0) {
    flux<< "  Wall particle type(s)  : ";
    for (auto& it : wallReference) flux<< it.second.getName() << " ";
    flux<< endl;
  }

  flux<< "  Interaction(s)         : " << endl << endl;
  flux<< "    * " << endl;
  for (auto& it : potentialReference) {
    flux << "    * " 
	 << setw(maxSize) << left << find(it.first.first)->getName() 
	 << " / " 
	 << setw(maxSize) << left << find(it.first.second)->getName() 
	 << " -> " 
	 << it.second->getName()
	 << ", with rcut : " << getRcut(it.first.first, it.first.second)
	 << endl;
  }
  for (auto& it : dpdReference) {
    flux << "    * "
	 << setw(maxSize) << left << find(it.first)->getName()
	 << " / "
	 << setw(maxSize) << left << find(it.second)->getName()
	 << " -> DPD "
	 << endl;
  }
  for (auto& it : dpdeReference) {
    flux << "    * "
	 << setw(maxSize) << left << find(it.first)->getName()
	 << " / "
	 << setw(maxSize) << left << find(it.second)->getName()
	 << " -> DPDE "
	 << endl;
  }
  for (auto& it : sdpdReference) {
    flux << "    * "
	 << setw(maxSize) << left << find(it.first)->getName()
	 << " / "
	 << setw(maxSize) << left << find(it.second)->getName()
	 << " -> SDPD "
	 << endl;
  }

  for (auto& it : sdpdWallReference) {
    flux << "    * "
	 << setw(maxSize) << left << find(it.first)->getName()
	 << " / "
	 << setw(maxSize) << left << find(it.second)->getName()
	 << " -> Wall (SDPD) "
	 << endl;
  }
  flux<< "    * " << endl;

  if (isDPD() || isSDPD()) {
    flux<< "  Equation(s) of state   : " << endl << endl;
    flux<< "    * " << endl;
    for (auto& it : eosReference) {
      flux << "    * "
	   << setw(maxSize) << left << find(it.first)->getName()
	   << " ->  " << setw(maxSize) << left << it.second->getName()
	   << endl;
    }
    flux<< "    * " << endl;
  }

  if (isReactive()) {
    flux<< "  Chemical reaction(s)   : " << endl << endl;
    flux<< "    * " << endl;
    for (auto& it : reactionReference) {
      flux << "    * "
	   << setw(maxSize) << left << find(it.first)->getName()
	   << " ->  " << setw(maxSize) << left << it.second->getName()
	   << endl;
    }
    flux<< "    * " << endl;
  }
  
}
