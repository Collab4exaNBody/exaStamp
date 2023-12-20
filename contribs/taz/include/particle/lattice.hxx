/// @file
/// @brief Templated functions for the lattice

// >>>
// #define INIT_CROP
// #define PUSH_X
// #define BIG_CRASH
// #define SIN_DEFECT
// <<<

/// @brief Initialize the particles in a subset of the lattice
/// @param [in] begin First cell of the subset
/// @param [in] end Last cell of the subset
/// @param [out] particles Initialized particles
/// @param [in] minBounds Lower bounds of the system
/// @param [in] maxBounds Upper bounds of the system
template <class Exchange_P>
void Configuration<Particle>::buildSubset(const uint64_t begin, const uint64_t end, std::vector<Exchange_P>& particles, const vec3<double>& minBounds, const vec3<double>& maxBounds) const {

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

	Exchange_P p;
	
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
	
	// Add this atom to the subset
	if (Global::reference.find(p.ti)!=nullptr)
	  particles.emplace_back(p);

      }
      
    }
  }
  
}

template <> 
void Configuration<Particle>::buildSubset<MPI__Mesoparticle>(const uint64_t begin, const uint64_t end, std::vector<MPI__Mesoparticle>& particles, const vec3<double>& minBounds, const vec3<double>& maxBounds) const;

template <> 
void Configuration<Particle>::buildSubset<MPI__SmoothParticle>(const uint64_t begin, const uint64_t end, std::vector<MPI__SmoothParticle>& particles, const vec3<double>& minBounds, const vec3<double>& maxBounds) const;
