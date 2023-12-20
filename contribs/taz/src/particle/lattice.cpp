/// @file 
/// @brief Implementation of the tools to initialize the particles


#include "globals.hpp"

#include "particle/lattice.hpp"


/// @brief Build a face-centered cubic lattice
/// @param [in] parameter Space parameter of the lattice
/// @param [in] types_ Types of the atoms in the lattice
Lattice* initFCCLattice(double parameter, Array<int> types_) {

  Lattice* fcc = new Lattice(1, 4);

  fcc->parameters[0] = parameter;

  for (int i=0; i<fcc->numberOfAtoms; ++i) 
    fcc->types[i] = types_[i];

  fcc->base[0] = vec3<double>(.25, .25, .25);
  fcc->base[1] = vec3<double>(.25, .75, .75);
  fcc->base[2] = vec3<double>(.75, .25, .75);
  fcc->base[3] = vec3<double>(.75, .75, .25);

  return fcc;

}


/// @brief Build a body centered cubic lattice
/// @param [in] parameter Space parameter of the lattice
/// @param [in] types_ Types of the atoms in the lattice
Lattice* initBCCLattice(double parameter, Array<int> types_) {

  Lattice* bcc = new Lattice(1, 2);

  bcc->parameters[0] = parameter;

  for (int i=0; i<bcc->numberOfAtoms; ++i) 
    bcc->types[i] = types_[i];

  bcc->base[0] = vec3<double>(.25, .25, .25);
  bcc->base[1] = vec3<double>(.75, .75, .75);

  return bcc;

}


/// @brief Build a diam 100 lattice
/// @param [in] parameter Space parameter of the lattice
/// @param [in] types_ Types of the atoms in the lattice
Lattice* initDIAM100Lattice(double parameter, Array<int> types_) {
  Lattice* diam100 = new Lattice(1, 8);

  diam100->parameters[0] = parameter;

  for (int i=0; i<diam100->numberOfAtoms; ++i) 
    diam100->types[i] = types_[i];

  diam100->base[0] = vec3<double>(.20, .20, .20);
  diam100->base[1] = vec3<double>(.70, .70, .20);
  diam100->base[2] = vec3<double>(.70, .20, .70);
  diam100->base[3] = vec3<double>(.20, .70, .70);
  diam100->base[4] = vec3<double>(.45, .45, .45);
  diam100->base[5] = vec3<double>(.95, .95, .45);
  diam100->base[6] = vec3<double>(.95, .45, .95);
  diam100->base[7] = vec3<double>(.45, .95, .95);
  return diam100;

}

/// @brief Build a simple cubic lattice
/// @param [in] parameter Space parameter of the lattice
/// @param [in] types_ Types of the atoms in the lattice
Lattice* initSCLattice(double parameter, Array<int> types_) {

  Lattice* sc = new Lattice(1, 1);

  sc->parameters[0] = parameter;

  for (int i=0; i<sc->numberOfAtoms; ++i) 
    sc->types[i] = types_[i];

  sc->base[0] = vec3<double>(.5, .5, .5);

  return sc;

}

/// @brief Constructor from some input data
/// @param [in] input Input data
Configuration<Particle>::Configuration(const Input& input)
  : m_initType(),
    m_initStep(-1),
    lattice(nullptr),
    wallLowerTypes(""),
    wallUpperTypes("") {

  // Get the number of lattice cells and initial temperature from the input
  numberOfLatticeCells = input.numberOfCells;
  initialTemperature   = input.initialTemperature;
  initialTint          = input.initialTint;

  // Get the atoms type to put in the lattice and convert them to type indexes
  Array<int> latticeTypes(input.latticeAtoms.size());
  for (uint i=0; i<latticeTypes.size(); ++i) {
    latticeTypes[i] = Global::reference.findTypeIndex(input.latticeAtoms[i]);
  }

  // Switch on the initialization type
  switch (input.m_initType) {
    
    // Case of an initialization from a Stamp output file
    case InitType::INIT_STAMP_LEGACY_DUMP :
    {
      // std::cerr << "PreInit createDomains legacy ... " <<  std::endl;
      m_initType = INIT_STAMP_LEGACY_DUMP;
      // Get the file name and initial step
      filename = input.initFile;
      m_initStep = input.initStep;
      // Reconstruct the file name if necessary
      if (filename == ""){
      	std::string tmp = std::to_string(m_initStep);
      	while (tmp.size() < 9) tmp.insert(0,"0");
      	filename="StampV3prot_"+tmp+".MpiIO";
      }
      // Get the number of particles
      numberOfParticles = input.legacyHeader->particlesTotalNumber;
      break;
    }
    
    // Case of an initialization from a Hercule output file
    case InitType::INIT_HERCULE_DUMP:
    {
      // std::cerr << "PreInit createDomains hercule ... " <<  std::endl;
      m_initType = INIT_HERCULE_DUMP;
      // Case of an initialization from a Stamp output file
      filename = input.initFile;
      m_initStep = input.initStep;
      // Get the number of particles
      numberOfParticles = input.legacyHeader->particlesTotalNumber;
      break;
    }
    
    // Case of an initialization from a stampv4 file
    case InitType::INIT_STAMPV4 :
      m_initType = INIT_STAMPV4;
      break;
      
    // Default initialization from a lattice
    case InitType::INIT_DEFAULT :
    default :
    {
      m_initType = INIT_DEFAULT;
      switch(input.latticeType) {
      // Case of a face-centered cubic lattice
      case Input::FCC : 
      	lattice = initFCCLattice(input.latticeParameter, latticeTypes);
      	NMAX = auxFloor<int>(Global::domainInfo.getExtension()/lattice->parameters[0]);
      	NCELLS = auxMin(numberOfLatticeCells, NMAX);
      	break;
      // Case of a body centered cubic lattice
      case Input::BCC : 
      	lattice = initBCCLattice(input.latticeParameter, latticeTypes);
      	NMAX = auxFloor<int>(Global::domainInfo.getExtension()/lattice->parameters[0]);
      	NCELLS = auxMin(numberOfLatticeCells, NMAX);
      	break;
      // Case of a body centered cubic lattice
      case Input::DIAM100 : 
      	lattice = initDIAM100Lattice(input.latticeParameter, latticeTypes);
      	NMAX = auxFloor<int>(Global::domainInfo.getExtension()/lattice->parameters[0]);
      	NCELLS = auxMin(numberOfLatticeCells, NMAX);
      	break;
      // Case of a simple cubic lattice
      case Input::SC : 
      	lattice = initSCLattice(input.latticeParameter, latticeTypes);
      	NMAX = auxFloor<int>(Global::domainInfo.getExtension()/lattice->parameters[0]);
      	NCELLS = auxMin(numberOfLatticeCells, NMAX);
      	break;
      }

      // Manage walls
      // Left wall
      vec3<bool> isWall = input.boundaryConditions >= Input::WALL && input.boundaryConditions <= Input::WALL;
      isWall = isWall || (input.boundaryConditions >= Input::WALL_FREE && input.boundaryConditions <= Input::WALL_FREE);
      wallMinBounds = input.origin + input.wallWidth*isWall;

      // Right wall
      isWall = input.boundaryConditions >= Input::WALL && input.boundaryConditions <= Input::WALL;
      isWall = isWall || (input.boundaryConditions >= Input::FREE_WALL && input.boundaryConditions <= Input::FREE_WALL);
      wallMaxBounds = input.origin + input.extension - input.wallWidth*isWall;

      wallLowerTypes = input.wallLowerTypes;
      wallUpperTypes = input.wallUpperTypes;

      break;
    }
    
  }

}


/// @brief Get the total number of particles in the system
/// @return Number of particles
uint64_t Configuration<Particle>::totalNumberOfParticles() {

  switch (m_initType) {

  case INIT_DEFAULT :
    return static_cast<uint64_t>(lattice->numberOfAtoms * (uint64_t) NCELLS.x * (uint64_t) NCELLS.y * (uint64_t) NCELLS.z);
    break;

  case INIT_STAMP_LEGACY_DUMP :
    return numberOfParticles;
    break;

  case INIT_HERCULE_DUMP :
    return numberOfParticles;
    break;

  case INIT_STAMPV4 :
  	std::cout << "Error ! This configuration should not be used in a molecular initialization." << std::endl;
  	return -1;
  	break;

  }

  return -1;

}



// Specializations of buildSubset

/// @brief Initialize the mesoparticles in a subset of the lattice
/// @param [in] begin First cell of the subset
/// @param [in] end Last cell of the subset
/// @param [out] particles Initialized particles
/// @param [in] minBounds Lower bounds of the system
/// @param [in] maxBounds Upper bounds of the system
template <>
void Configuration<Particle>::buildSubset<MPI__Mesoparticle>(const uint64_t begin, const uint64_t end, std::vector<MPI__Mesoparticle>& particles, const vec3<double>& minBounds, const vec3<double>& maxBounds) const {

  // Get the center of the system
  vec3<double> center = Global::domainInfo.getMinBounds() + 0.5*Global::domainInfo.getExtension() - 0.5*NCELLS*lattice->parameters[0];

  // Function to check if a coordinate is inside the system
  auto inBounds = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= minBounds;
    vec3<bool> M = r <  maxBounds;
    return m.x && m.y && m.z && M.x && M.y && M.z;
  };

  // Get the boundary conditions
  auto BC = Global::domainInfo.getBoundaryConditions();

  const uint64_t nyz = (uint64_t) NCELLS.y * (uint64_t) NCELLS.z;
  const uint64_t ny  = (uint64_t) NCELLS.y;
  const uint64_t nz  = (uint64_t) NCELLS.z;

  // Function to convert the index of a cell into a coordinate
  auto index_convert = [&] (const uint64_t& c) -> vec3<double> {
    uint64_t x = c/nyz;
    uint64_t y = (c/nz) % ny;
    uint64_t z = c%nz;
    return vec3<double>( (double) x, (double) y, (double) z);
  };

  // Set the sigma parameter necessary to initialize the velocity
  Array<double> sigma(Global::reference.getNumberOfTypes());

  for (uint i=0; i<sigma.size(); ++i)
    sigma[i] = auxSqrt( Stamp_Constant::boltzmann * initialTemperature / Global::reference.find(i)->getMass());
  

  // Some special case where only a part of the system is initialized (70% in each dimension)
#ifdef INIT_CROP

  vec3<double> gMin = Global::domainInfo.getMinBounds() + Global::domainInfo.getExtension()*vec3<double>(0.10, 0.10, 0.10);
  vec3<double> gMax = gMin + Global::domainInfo.getExtension()*vec3<double>(0.8, 0.8, 0.8);

  auto inCube = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= gMin;
    vec3<bool> M = r <  gMax;
    return m.x && m.y && m.z && M.x && M.y && M.z;
  };

#endif
  

  // Some special case where only a part of the system is initialized (approximatively 5% of x and 97% of y and z)
#ifdef BIG_CRASH

  vec3<double> gMin1 = Global::domainInfo.getMinBounds()+Global::domainInfo.getExtension()*vec3<double>(0.00, 0.01, 0.01)+Global::domainInfo.getCellLength()*vec3<double>(1,0,0);
  vec3<double> gMax1 = gMin1 + Global::domainInfo.getExtension()*vec3<double>(0.05, 0.98, 0.98);

  vec3<double> gMin2 = vec3<double>(gMax1.x, gMin1.y, gMin1.z);
  vec3<double> gMax2 = gMin2 + Global::domainInfo.getExtension()*vec3<double>(0.04, 0.98, 0.98);

  vec3<double> gMin3 = Global::domainInfo.getMinBounds() + Global::domainInfo.getExtension()*vec3<double>(0.495, 0.01, 0.01);
  vec3<double> gMax3 = gMin3 + Global::domainInfo.getExtension()*vec3<double>(0.01, 0.98, 0.98);

  auto inCube1 = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= gMin1;
    vec3<bool> M = r <  gMax1;
    return m.x && m.y && m.z && M.x && M.y && M.z;
  };

  auto inCube2 = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= gMin2;
    vec3<bool> M = r <  gMax2;
    bool inMainCube = m.x && m.y && m.z && M.x && M.y && M.z;
    vec3<double> gCenter2 = gMin2*vec3<double>(0.0, 0.5, 0.5) + gMax2*vec3<double>(1.0, 0.5, 0.5);
    bool def = inMainCube && (norm2(r-gCenter2)>0.25*auxSq(gMax2.x-gMin2.x));
    return def;
  };

  auto inCube3 = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= gMin3;
    vec3<bool> M = r <  gMax3;
    return m.x && m.y && m.z && M.x && M.y && M.z;
  };

#endif
  

    // Some special case where a sinusoidal default is created in the upper x-surface.
#ifdef SIN_DEFECT

  vec3<double> xyzMin = Global::domainInfo.getMinBounds();
  vec3<double> xyzMax = Global::domainInfo.getMaxBounds();

  double yCenter = (xyzMax.y + xyzMin.y)*0.5;
  double yWidth = 15.;
  double xDepth = 5.;
  
  auto outOfDefault = [&] (const vec3<double>& r) -> bool {

    double fy = xDepth * cos( 0.5*Constant::pi/yWidth * (r.y-yCenter) );
    
    return r.x < xyzMax.x - fy;
  };

#endif


  // Tests for all the cases
  auto theTest = [&] (const vec3<double>& r) -> bool {
#if defined INIT_CROP
    return inCube(r) && inBounds(r);
#elif defined BIG_CRASH
    return (inCube1(r) || inCube2(r) || inCube3(r)) && inBounds(r);
#elif defined SIN_DEFECT
    return outOfDefault(r) && inBounds(r);
#else
    return inBounds(r);
#endif
  };


  for (uint64_t c=begin; c<end; ++c) {
    for (int atom=0; atom<lattice->numberOfAtoms; ++atom) {
    
      vec3<double> r = center + lattice->parameters[0]*(index_convert(c) + lattice->base[atom]);

      const uint64_t id = c * (uint64_t) lattice->numberOfAtoms + atom;

      if (theTest(r)) {

	MPI__Mesoparticle p;

	p.id = id;

	// Check whether the particle is in a wall and set type accordingly
	bool typeSet = false;
	for (uint dim = 0; dim < 3 && !typeSet; dim++) {
	  switch (BC[dim]) {
	  case Input::WALL:
	    if (r[dim] < wallMinBounds[dim]) {
	      p.ti = Global::reference.findTypeIndex(wallLowerTypes[dim]);
	      typeSet = true;
	    }
	    if (r[dim] > wallMaxBounds[dim]) {
	      p.ti = Global::reference.findTypeIndex(wallUpperTypes[dim]);
	      typeSet = true;
	    }
	    break;

	  case Input::WALL_FREE:
	    if (r[dim] < wallMinBounds[dim]) {
	      p.ti = Global::reference.findTypeIndex(wallLowerTypes[dim]);
	      typeSet = true;
	    }
	    break;
	    
	  case Input::FREE_WALL:
	    if (r[dim] > wallMaxBounds[dim]) {
	      p.ti = Global::reference.findTypeIndex(wallUpperTypes[dim]);
	      typeSet = true;
	    }
	    break;

	  case Input::PERIODIC:
	  case Input::FREE:
	  default:
	    p.ti = lattice->types[atom];
	    break;
	  }
	}

	
	// Get the velocity for this atom
	vec3<double> v = Global::reference.getVelocity(p.ti);
	if (Global::reference.getInvMass(p.ti) != 0) {
#ifndef PUSH_X
	  v += sigma[lattice->types[atom]] * vec3<double>(Saru::randN(id, -1), Saru::randN(id, -2), Saru::randN(id, -3));
#else
	  v += vec3<double>(0., 0., +8.0) + sigma[lattice->types[atom]] * vec3<double>(Saru::randN(id, -1), Saru::randN(id, -2), Saru::randN(id, -3));
#endif

#ifdef BIG_CRASH
	  if (inCube2(r)) v += sigma[lattice->types[atom]] * vec3<double>(-10.0, 0.0, 0.0);
#endif
	}

	p.r  = r;
	p.v  = v;
	p.ei = 0.;
	p.progress = 0.;

	// Add this atom to the subset
	if (Global::reference.find(p.ti)!=nullptr)
	  particles.emplace_back(p);
	
      }
      
    }
  }

}



/// @brief Initialize the smooth particles in a subset of the lattice
/// @param [in] begin First cell of the subset
/// @param [in] end Last cell of the subset
/// @param [out] particles Initialized particles
/// @param [in] minBounds Lower bounds of the system
/// @param [in] maxBounds Upper bounds of the system
template <> 
void Configuration<Particle>::buildSubset<MPI__SmoothParticle>(const uint64_t begin, const uint64_t end, std::vector<MPI__SmoothParticle>& particles, const vec3<double>& minBounds, const vec3<double>& maxBounds) const {

  vec3<double> center = Global::domainInfo.getMinBounds() + 0.5*Global::domainInfo.getExtension() - 0.5*NCELLS*lattice->parameters[0];

  auto inBounds = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= minBounds;
    vec3<bool> M = r <  maxBounds;
    return m.x && m.y && m.z && M.x && M.y && M.z;
  };

  // Get the boundary conditions
  auto BC = Global::domainInfo.getBoundaryConditions();

  const uint64_t nyz = (uint64_t) NCELLS.y * (uint64_t) NCELLS.z;
  const uint64_t ny  = (uint64_t) NCELLS.y;
  const uint64_t nz  = (uint64_t) NCELLS.z;

  auto index_convert = [&] (const uint64_t& c) -> vec3<double> {
    uint64_t x = c/nyz;
    uint64_t y = (c/nz) % ny;
    uint64_t z = c%nz;
    return vec3<double>( (double) x, (double) y, (double) z);
  };

  Array<double> sigma(Global::reference.getNumberOfTypes());

  for (uint i=0; i<sigma.size(); ++i)
    sigma[i] = auxSqrt( Stamp_Constant::boltzmann * initialTemperature / Global::reference.find(i)->getMass());

    // Some special case where only a part of the system is initialized (70% in each dimension)
#ifdef INIT_CROP

  vec3<double> gMin = Global::domainInfo.getMinBounds() + Global::domainInfo.getExtension()*vec3<double>(0.10, 0.10, 0.10);
  vec3<double> gMax = gMin + Global::domainInfo.getExtension()*vec3<double>(0.8, 0.8, 0.8);

  auto inCube = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= gMin;
    vec3<bool> M = r <  gMax;
    return m.x && m.y && m.z && M.x && M.y && M.z;
  };

#endif
  

  // Some special case where only a part of the system is initialized (approximatively 5% of x and 97% of y and z)
#ifdef BIG_CRASH

  vec3<double> gMin1 = Global::domainInfo.getMinBounds()+Global::domainInfo.getExtension()*vec3<double>(0.00, 0.01, 0.01)+Global::domainInfo.getCellLength()*vec3<double>(1,0,0);
  vec3<double> gMax1 = gMin1 + Global::domainInfo.getExtension()*vec3<double>(0.05, 0.98, 0.98);

  vec3<double> gMin2 = vec3<double>(gMax1.x, gMin1.y, gMin1.z);
  vec3<double> gMax2 = gMin2 + Global::domainInfo.getExtension()*vec3<double>(0.04, 0.98, 0.98);

  vec3<double> gMin3 = Global::domainInfo.getMinBounds() + Global::domainInfo.getExtension()*vec3<double>(0.495, 0.01, 0.01);
  vec3<double> gMax3 = gMin3 + Global::domainInfo.getExtension()*vec3<double>(0.01, 0.98, 0.98);

  auto inCube1 = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= gMin1;
    vec3<bool> M = r <  gMax1;
    return m.x && m.y && m.z && M.x && M.y && M.z;
  };

  auto inCube2 = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= gMin2;
    vec3<bool> M = r <  gMax2;
    bool inMainCube = m.x && m.y && m.z && M.x && M.y && M.z;
    vec3<double> gCenter2 = gMin2*vec3<double>(0.0, 0.5, 0.5) + gMax2*vec3<double>(1.0, 0.5, 0.5);
    bool def = inMainCube && (norm2(r-gCenter2)>0.25*auxSq(gMax2.x-gMin2.x));
    return def;
  };

  auto inCube3 = [&] (const vec3<double>& r) -> bool {
    vec3<bool> m = r >= gMin3;
    vec3<bool> M = r <  gMax3;
    return m.x && m.y && m.z && M.x && M.y && M.z;
  };

#endif


    // Some special case where a sinusoidal default is created in the upper x-surface.
#ifdef SIN_DEFECT

  vec3<double> xyzMin = Global::domainInfo.getMinBounds();
  vec3<double> xyzMax = Global::domainInfo.getMaxBounds();

  double yCenter = (xyzMax.y + xyzMin.y)*0.5;
  double yWidth = 154.62;
  double xDepth = 50.;
  
  auto outOfDefault = [&] (const vec3<double>& r) -> bool {

    double fy = 0.5*xDepth * (1 + cos( 2*Constant::pi/yWidth * (r.y-yCenter) ));
    
    return r.x < xyzMax.x - fy;
  };

#endif


  // Tests for all the cases
  auto theTest = [&] (const vec3<double>& r) -> bool {
#if defined INIT_CROP
    return inCube(r) && inBounds(r);
#elif defined BIG_CRASH
    return (inCube1(r) || inCube2(r) || inCube3(r)) && inBounds(r);
#elif defined SIN_DEFECT
    return outOfDefault(r) && inBounds(r);
#else
    return inBounds(r);
#endif
  };

  // Loop on the cells of the subset
  for (uint64_t c=begin; c<end; ++c) {

    // Loop on the atoms of a cell
    for (int atom=0; atom<lattice->numberOfAtoms; ++atom) {

      // Get a position for this atom
      vec3<double> r = center + lattice->parameters[0]*(index_convert(c) + lattice->base[atom]);

      // Get an index for this atom
      const uint64_t id = c * (uint64_t) lattice->numberOfAtoms + atom;

      // If the atom is in the system
      if (theTest(r)) {

	MPI__SmoothParticle p;
	
	p.id = id;

	// Check whether the particle is in a wall and set type accordingly
	bool typeSet = false;
	for (uint dim = 0; dim < 3 && !typeSet; dim++) {
	  switch (BC[dim]) {
	  case Input::WALL:
	    if (r[dim] < wallMinBounds[dim]) {
	      p.ti = Global::reference.findTypeIndex(wallLowerTypes[dim]);
	      typeSet = true;
	    }
	    if (r[dim] > wallMaxBounds[dim]) {
	      p.ti = Global::reference.findTypeIndex(wallUpperTypes[dim]);
	      typeSet = true;
	    }
	    break;

	  case Input::WALL_FREE:
	    if (r[dim] < wallMinBounds[dim]) {
	      p.ti = Global::reference.findTypeIndex(wallLowerTypes[dim]);
	      typeSet = true;
	    }
	    break;
	    
	  case Input::FREE_WALL:
	    if (r[dim] > wallMaxBounds[dim]) {
	      p.ti = Global::reference.findTypeIndex(wallUpperTypes[dim]);
	      typeSet = true;
	    }
	    break;

	  case Input::PERIODIC:
	  case Input::FREE:
	  default:
	    p.ti = lattice->types[atom];
	    break;
	  }
	}

	// Get the velocity for this atom
	vec3<double> v = Global::reference.getVelocity(p.ti);
	if (Global::reference.getInvMass(p.ti) != 0) {
#ifndef PUSH_X
	  v += sigma[lattice->types[atom]] * vec3<double>(Saru::randN(id, -1), Saru::randN(id, -2), Saru::randN(id, -3));
#else
	  v += vec3<double>(0., 0., +8.0) + sigma[lattice->types[atom]] * vec3<double>(Saru::randN(id, -1), Saru::randN(id, -2), Saru::randN(id, -3));
#endif

#ifdef BIG_CRASH
	  if (inCube2(r)) v += sigma[lattice->types[atom]] * vec3<double>(-10.0, 0.0, 0.0);
#endif
	}

	p.r  = r;
	p.v  = v;
	p.ei = 0.;
	p.progress = 0;
	
	// Add this particle to the subset
	// Add this atom to the subset
	if (Global::reference.find(p.ti)!=nullptr)
	  particles.emplace_back(p);
	p.ei = 0;
      }

    }
  }

}
