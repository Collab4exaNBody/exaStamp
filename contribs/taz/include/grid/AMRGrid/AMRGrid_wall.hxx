#pragma once

/// @brief Apply a repulsive force on each particle near a specified wall boundary
/// @param [in] dim Dimension of the boundary : 0 for x, 1 for y, 2 for z
/// @param [in] direction Direction of the boundary : upper bound (1) or lower bound (-1)
TMPLSG void TMPL_AMRGrid::applyWallCondition(const int dim, const int direction) {

  const auto bound = (direction<0) ? info.getCellInf()[dim] : info.getCellSup()[dim];
  const auto limit = (direction<0) ? 0 : Global::domainInfo.getNumberOfCellsPerDim()[dim];
  

  if (bound!=limit) return;

  // Work on edge cells only (as other cannot be next to a wall)
  const auto& cells = this->getTraversal(TraversalManager::EDGE_ONE);

  #pragma omp parallel for 
	for(size_t o=0; o<cells.size(); ++o) {
	
		// Check if cell is on specified boundary
		if ( (direction==-1 && this->onMinBounds(cells[o], dim)) ||
				(direction==+1 && this->onMaxBounds(cells[o], dim)) )
				
	  {
	       
      const double globMinBounds = Global::domainInfo.getMinBounds()  [dim];
      const double globMaxBounds = Global::domainInfo.getMaxBounds()  [dim];
      const double maxRcut       = Global::reference.getMaxRcut();
      
      // Access the good component of the positions and forces
      double* r_ = nullptr;
      double* f_ = nullptr;
      
      double* ep = particles[cells[o]].getPotentialEnergy();
      uint8_t* ti = particles[cells[o]].getType();
      
      switch (dim) {
        case 0: r_ = particles[cells[o]].getPositionX(); f_ = particles[cells[o]].getForceX(); break;
        case 1: r_ = particles[cells[o]].getPositionY(); f_ = particles[cells[o]].getForceY(); break;
        case 2: r_ = particles[cells[o]].getPositionZ(); f_ = particles[cells[o]].getForceZ(); break;
        default:                                break;
      }
      
      // Get the position of the wall
      const double rWall  = (direction<0) ? globMinBounds : globMaxBounds;
      
      // Compute and add the force and potential energy for all the particles
      for (uint i=0; i<particles[cells[o]].getNumberOfParticles(); ++i) {
        if (Global::reference.getInvMass(ti[i]) != 0) {
          double r     = auxAbs(r_[i]-rWall);
          
          if(r>maxRcut) continue;
          
          double ir    = 1./r;
          double ratio = (1-maxRcut*ir);
          double tmp   = ratio*ratio;
          tmp = tmp*tmp;
          tmp = tmp*tmp*tmp;
          
          ep[i] += tmp;
          f_[i] += direction*12.0 * maxRcut*ir*ir * tmp/ratio;

        }
      }
	  }
	}


}


/// @brief Check if there is escaped particles at free boundaries
/// @param [out] outOfFreeBounds Presence of escaped particles for each free boundary
TMPLSG void TMPL_AMRGrid::checkFreeBoundaries(Array<int>& outOfFreeBounds) {

  const vec3<bool> testInf = Global::domainInfo.isFreeInf();
  const vec3<bool> testSup = Global::domainInfo.isFreeSup();

  const auto& boundInf = info.getCellInf();
  const auto& boundSup = info.getCellSup();
  const auto& limitSup = Global::domainInfo.getNumberOfCellsPerDim();

  // Work on edge cells only (as other cannot be next to a system boundary)
  const auto& cells = this->getTraversal(TraversalManager::EDGE_ONE);

  // Parallelize the work by distributing the cells between the threads
  parallel_region(0, cells.size(), cells.size()/thread_grain,
  		[&](const uint begin, const uint end) {

  	for (uint i=begin; i<end; ++i) {
  		for (uint8_t dim=0; dim<VEC3_NDIMS; ++dim) {

  			if ( testInf[dim] && (outOfFreeBounds[0+dim]==0) && (boundInf[dim]==0) && this->onMinBounds(cells[i], dim) ) {
  				if (particles[cells[i]].getNumberOfParticles()>0) {

  					this->lock(0+dim);
  					++outOfFreeBounds[0+dim];
  					this->unlock(0+dim);
  				}
  			}
  			if ( testSup[dim] && (outOfFreeBounds[3+dim]==0) && (boundSup[dim]==limitSup[dim]) && this->onMaxBounds(cells[i], dim) ) {
  				if (particles[cells[i]].getNumberOfParticles()>0) {

  					this->lock(3+dim);
  					++outOfFreeBounds[3+dim];
  					this->unlock(3+dim);
  				}
  			}

  		}
  	}

  });

}

/// @brief Stop the wall particles
TMPLSG void TMPL_AMRGrid::stopWalls() 
{

  #pragma omp parallel for
  for(size_t i=0; i<octrees.size(); ++i)
  {
    double* ptr_vx = octrees[i]->getVelocityX();
    double* ptr_vy = octrees[i]->getVelocityY();
    double* ptr_vz = octrees[i]->getVelocityZ();
    auto ptr_ti = octrees[i]->getType();
    
    size_t numberOfParticles = octrees[i]->getNumberOfParticles();
    
    for (size_t j = 0; j<numberOfParticles; ++j)
    {
      if (Global::reference.isWallType(ptr_ti[j]))
      {
        ptr_vx[j] = 0.;
        ptr_vy[j] = 0.;
        ptr_vz[j] = 0.;
      }
    }
  }  
}

