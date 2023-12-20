#pragma once

/// @brief Reset forces of all the particles of the grid
TMPLSG inline void TMPL_AMRGrid::resetForce() {

  #pragma omp parallel for
  for( int i = 0 ; i < octrees.size() ; i++)
  {
    size_t numberOfParticles = octrees[i]->getNumberOfParticles();
    simdAMR::kernels::reset(octrees[i]->getForceX(), numberOfParticles);
    simdAMR::kernels::reset(octrees[i]->getForceY(), numberOfParticles);
    simdAMR::kernels::reset(octrees[i]->getForceZ(), numberOfParticles);
  } 
}

/// @brief Reset energies of all the particles of the grid
TMPLSG inline void TMPL_AMRGrid::resetEnergy() {


  #pragma omp parallel for
  for( int i = 0 ; i < octrees.size() ; i++)
  {
    size_t numberOfParticles = octrees[i]->getNumberOfParticles();
    simdAMR::kernels::reset(octrees[i]->getPotentialEnergy(), numberOfParticles);
  } 
}


/// @brief Reset forces and energies of all the particles of the grid
///
///
TMPLSG inline void TMPL_AMRGrid::resetForceAndEnergy() {

  #pragma omp parallel for
  for( int i = 0 ; i < octrees.size() ; i++)
  {
    size_t numberOfParticles = octrees[i]->getNumberOfParticles();
    simdAMR::kernels::reset(octrees[i]->getForceX(), numberOfParticles);
    simdAMR::kernels::reset(octrees[i]->getForceY(), numberOfParticles);
    simdAMR::kernels::reset(octrees[i]->getForceZ(), numberOfParticles);
    simdAMR::kernels::reset(octrees[i]->getPotentialEnergy(), numberOfParticles);
  }  
}


