#pragma once

//#define USEBUFFER
//#define ANALYSE

/// @brief Compute force for specified potential, between specified particle types, on specified traversal (for a pair potential). Don't forget to change the linked cell version
///
/// Cast a Pair potential into its real type to call an optimized version
/// @param [in] potential Pair potential
/// @param [in] typeIndexA First particle type
/// @param [in] typeIndexB Second particle type
/// @param [in] T Working traversal
TMPLSG inline void TMPL_AMRGrid::computeForceVerlet(PairPotential* potential, uint8_t typeIndexA, uint8_t typeIndexB, ::Traversal T) {


  // Get potential type and cast it
  const auto& type = potential->getType();
  switch (type) {

  case IPotential::Type::LJ :
    __computeForcePairVerlet(static_cast<LennardJonesPotential*>(potential), T, typeIndexA, typeIndexB);
    break;

  default :
    __computeForcePairVerlet(potential, T, typeIndexA, typeIndexB);
    break;

  }

}

/// @brief Compute force for specified potential, between specified particle types, on specified traversal (for an EAM potential). Don't forget to change the linked cell version
///
/// @brief Cast an EAM potential into its real type to call an optimized version
/// @param [in] potential EAM potential
/// @param [in] typeIndexA First particle 
/// @param [in] typeIndexB Second particle type
/// @param [in] T Working traversal
TMPLSG inline void TMPL_AMRGrid::computeForceVerlet(EAMPotential* potential, uint8_t typeIndexA, uint8_t typeIndexB, ::Traversal T) {

  // Get potential type and cast it
  const auto& type = potential->getType();
  switch (type) {

  case IPotential::Type::SUTTON_CHEN :
    __computeForceEAMVerlet(static_cast<SuttonChenPotential*>(potential), T, typeIndexA, typeIndexB);
    break;

  case IPotential::Type::EAM_VNIITF :
    //__computeForceEAMVerlet(static_cast<EamVniitfPotential*>(potential), T, typeIndexA, typeIndexB);
    break;

  case IPotential::Type::MEAM_ :
    __computeForceMEAMVerlet(static_cast<MeamPotential*>(potential), T, typeIndexA, typeIndexB);
    break;

  default :
    std::cout << " no potential " << std::endl; abort(); //__computeForceEAMVerlet(potential, T, typeIndexA, typeIndexB);
    break;

  }

}


/// @brief Compute force for pair potential with Verlet lists. Don't forget to change the linked cell version
/// @tparam Pot_t Type of the potential
/// @param [in] pot Specified pair potential
/// @param [in,out] cells Working cells
/// @param [in] typeIndexA First particle type
/// @param [in] typeIndexB Second particle type
TMPLSG template <class Pot_t>
inline void TMPL_AMRGrid::__computeForcePairVerlet(Pot_t* pot, ::Traversal T, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  static const uint8_t A(typeIndexA);
  static const uint8_t B(typeIndexB); 

  auto ljOctreeBlocsVerlet = [&] (uint i)->void {
    octrees[i]->computeForcePair_perso_per_block_mutex(static_cast<LennardJonesPotential*>(pot), A, B); 
  };

  auto ljOctreeVerlet = [&] (uint i)->void {
    octrees[i]->computeForcePairVerlet_mutex(static_cast<LennardJonesPotential*>(pot), A, B); 
  };

  auto ljDep = [&] (uint idx_cell)->void { 
    if(Global::reference.isBlockVerlet()) particles[idx_cell]. computeForcePair_perso_per_block(static_cast<LennardJonesPotential*>(pot), typeIndexA, typeIndexB) ;
    else particles[idx_cell]. computeForcePair_perso(static_cast<LennardJonesPotential*>(pot), typeIndexA, typeIndexB) ;
  };

  int nbThread; 
# pragma omp parallel shared(nbThread)
  {
    nbThread=omp_get_num_threads();
  }

  const int nbWave = 8;

  assert(nbThread != 0);

#if _OPENMP >= 201511	

  launchGraph(ljDep);

#else // _OPENMP

  if(Global::reference.isBlockVerlet())
    {
  #pragma omp parallel for schedule(dynamic)
    for(int i =0; i<octrees.size(); i++)   
      ljOctreeBlocsVerlet(i);
    }
  else
    {
      #pragma omp parallel for schedule(dynamic)
      for(int i =0; i<octrees.size(); i++)
        ljOctreeVerlet(i);
    }
#endif // _OPENMP
}



/// @brief Reset EAM data (rho and embedding terms)
/// @param [in,out] cells Working cells
TMPLSG inline void TMPL_AMRGrid::__resetEAMData(const Array<uint>& cells) {

  // Parallel loop on cells
  parallel_region(0, cells.size(), [&](const uint begin, const uint end) {

      for (uint i=begin; i<end; ++i)
	{
	  particles[cells[i]].resetEAMData();
	  particles[cells[i]].checkEAMData();
	}

    }, force_partitioner);    

}


/// @brief Compute embedding terms for EAM Potential
/// @tparam Pot_t Type of the potential
/// @param [in] pot Specified EAM potential
/// @param [in,out] cells Working cells
/// @param [in] typeIndexA First particle type
/// @param [in] typeIndexB Second particle type
TMPLSG template <class Pot_t>
inline void TMPL_AMRGrid::__computeForceEmb(Pot_t* pot, const Array<uint>& cells, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  // Parallel loop on octrees
  parallel_region(0, octrees.size(), [&](const uint begin, const uint end) {

      for (uint i=begin; i<end; ++i) 
	octrees[i]->computeForceEmb(pot, typeIndexA, typeIndexB);

    }, force_partitioner);   
}


/// @brief Compute force for EAM potential. Don't forget to change the linked cell version
/// @tparam Pot_t Type of the potential
/// @param [in] potential Specified EAM potential
/// @param [in] T Working traversal
/// @param [in] typeIndexA First particle type
/// @param [in] typeIndexB Second particle type
TMPLSG template <class Pot_t>
void TMPL_AMRGrid::__computeForceEAMVerlet(Pot_t* potential, ::Traversal T, uint8_t typeIndexA, uint8_t typeIndexB) {

  // Get level to call right traversal
  static const uint8_t LEVEL = (uint8_t) auxMin(info.getGhostThickness()-1, (int) potential->getGhostThickness());

  static const uint8_t A(typeIndexA);
  static const uint8_t B(typeIndexB); 

  // Get potential kernels
  auto rho = [&] (uint idx_cell)->void {
    if(Global::reference.isBlockVerlet()) particles[idx_cell]. computeForceRhoVerlet_without_buffer_per_block (potential, A, B);
    else particles[idx_cell]. computeForceRhoVerlet_without_buffer (potential, A, B, this->domain->symmetrize());
  };


  auto finalFunc = [&] (uint idx_cell)->void {
    if(Global::reference.isBlockVerlet()) particles[idx_cell]. computeForceFinalVerlet_without_buffer_per_block (potential, A, B);
    else particles[idx_cell]. computeForceFinalVerlet_without_buffer (potential, A, B, this->domain->symmetrize()); 
  };

  // Get the working cells
  auto& cells = this->getTraversal(T, LEVEL);

  __resetEAMData(cells);

  // Compute rho terms
  launchGraph(rho);   

  // Compute embedding terms
  __computeForceEmb<Pot_t>(potential, cells, typeIndexA, typeIndexB);

  // Exchange embedding terms
  fillEmbSelf();
  sendEmb();
  collectEmb();

  // Compute the forces
  launchGraph(finalFunc);   
}

/// @brief Compute force for MEAM potential. Don't forget to change the linked cell version
/// @tparam Pot_t Type of the potential
/// @param [in] potential Specified MEAM potential
/// @param [in] T Working traversal
/// @param [in] typeIndexA First particle type
/// @param [in] typeIndexB Second particle type
TMPLSG template <class Pot_t>
void TMPL_AMRGrid::__computeForceMEAMVerlet(Pot_t* potential, ::Traversal T, uint8_t typeIndexA, uint8_t typeIndexB) {

  static const uint8_t A(typeIndexA);
  static const uint8_t B(typeIndexB); 


  constexpr bool mutex_on = true;
  constexpr bool mutex_off = false;

  auto meam = [&] (uint idx_cell)->void {
    particles[idx_cell].computeForcesAndPotentialMEAMVerlet_AMR<mutex_off>(
	potential, 
	A, 
	B, 
	idx_cell, 
	this->ghostLayer.data()
	); 
  };


  auto meamOctree = [&] (uint idx)->void {
    particles[idx].computeForcesAndPotentialMEAMVerlet_AMR<mutex_on>(
	potential, 
	A, 
	B, 
	idx, 
	this->ghostLayer.data() 
	);
  };

#pragma omp parallel for schedule(dynamic)
	for(int i =0; i<octreesMEAM.size(); i++) 
	  meamOctree(octreesMEAM[i]); 

}

