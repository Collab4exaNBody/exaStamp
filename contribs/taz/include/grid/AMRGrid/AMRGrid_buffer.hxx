#pragma once


/// @brief Fill a buffer with all the particles of the node, case of a ParticleOutput buffer
/// @brief [in] Buffer to fill
TMPLSG void TMPL_AMRGrid::fillBuffer(ParticleOutput* buffer) {
    
  const auto& cells = this->getTraversal(TraversalManager::REAL);

  Array<uint> start(cells.size(), 0);

  // set start for each cell
  for (uint i=1, size=cells.size(); i<size; ++i) {
    start[i] = start[i-1] + particles[ cells[i-1] ].getNumberOfParticles();
  }

  parallel_region(0, cells.size(), [&](const uint begin, const uint end) {

    for(uint i=begin; i<end; ++i)
	    particles[ cells[i] ].fillBuffer(buffer, start[i]);

    });
  
}

/// @brief Fill a buffer with all the particles of the node, case of a ParticleInSitu buffer
/// @brief [in] Buffer to fill
TMPLSG void TMPL_AMRGrid::fillBuffer(ParticleInSitu* buffer) {
    
  const auto& cells = this->getTraversal(TraversalManager::REAL);

  Array<uint> start(cells.size(), 0);

  // set start for each cell
  for (uint i=1, size=cells.size(); i<size; ++i) {
    start[i] = start[i-1] + particles[ cells[i-1] ].getNumberOfParticles();
  }

  // fill the particles attributes
  parallel_region(0, cells.size(), [&](const uint begin, const uint end) {

      for(uint i=begin; i<end; ++i)
	particles[ cells[i] ].fillBuffer(buffer, start[i]);

    });

}

/// @brief Fill a buffer with the ghost information of the node, case of a ParticleInSitu buffer
/// @brief [in] Buffer to fill
TMPLSG void TMPL_AMRGrid::fillGhostBuffer(ParticleInSitu* buffer) {

  // Get cells information
  const auto& cells = this->getTraversal(TraversalManager::REAL);
  Array<uint> start(cells.size(), 0);
  for (uint i=1, size=cells.size(); i<size; ++i) {
    start[i] = start[i-1] + particles[ cells[i-1] ].getNumberOfParticles();
  }
  
  // Construction of ghost data
  const auto& ghostCells = this->getTraversal(TraversalManager::GHOST);
  Array<uint> ghostStart(ghostCells.size(), 0);
  uint size = ghostCells.size();
  for (uint i=1; i<size; ++i) {
    ghostStart[i] = ghostStart[i-1] + particles[ ghostCells[i-1] ].getNumberOfParticles();
  }

  // Get number of ghost particles for allocation
  int nbGhostParticles = ghostStart[size-1] + particles[ ghostCells[size-1] ].getNumberOfParticles();
  buffer->setNbGhost(nbGhostParticles);

  if (buffer->ghostAllocated)
    {

      // fill the cells data for the particles
      for (uint i=0, size=cells.size(); i<size; ++i)
        {
          buffer->particles_start_size[ cells[i] ][0] = start[i]; // fix for q param analytics
          buffer->particles_start_size[ cells[i] ][1] = particles[ cells[i] ].getNumberOfParticles();
          for (uint j=0; j<Neighbor::num_neighbors; ++j)
            buffer->particles_neigh[ cells[i] ][j] = this->neighbors[ cells[i] ][j];
        }

      // fill the ghost positions
      parallel_region(0, ghostCells.size(), [&](const uint begin, const uint end) {

          for(uint i=begin; i<end; ++i)
            particles[ ghostCells[i] ].fillGhostBuffer(buffer, ghostStart[i]);

        });

      // fill the ghost cells data
      for (uint i=0, size=ghostCells.size(); i<size; ++i)
        {
          buffer->ghost_particles_start_size[ ghostCells[i] ][0] = ghostStart[i] + buffer->nPart;
          buffer->ghost_particles_start_size[ ghostCells[i] ][1] = particles[ ghostCells[i] ].getNumberOfParticles();
        }
      
    }
  
}


/// @brief Fill a buffer with all the particles of the node, case of a legacy Stamp dump buffer
/// @brief [in] Buffer to fill
/// @tparam DumpStruct Dump structure (legacy or legacy for DPDE)
/// @param [out] buffer Particles data to fill
// Needed to fill the legacy Stamp V3 dump
TMPLSG template <class DumpStruct>
void TMPL_AMRGrid::fillBuffer(Array<DumpStruct>& buffer) {

  const auto& cells = this->getTraversal(TraversalManager::REAL);

  Array<uint> start(cells.size(), 0);

  // set start for each cell
  for (uint i=1, size=cells.size(); i<size; ++i) {
    start[i] = start[i-1] + particles[ cells[i-1] ].getNumberOfParticles();
  }
  
  parallel_region(0, cells.size(), [&](const uint begin, const uint end) {

      for(uint i=begin; i<end; ++i)
	particles[ cells[i] ].fillBufferAMR(buffer, start[i]);

    });
  
}


/// @brief Fill a buffer with all the particles of the node, case of a Hercule dump buffer
/// @brief [in] Buffer to fill
/// @tparam DumpStruct Dump structure (legacy or legacy for DPDE)
/// @param [out] buffer Particles data to fill
// Needed to fill the hercule dump
TMPLSG template <class DumpStruct>
void TMPL_AMRGrid::fillBuffer(DumpStruct& buffer) {

  const auto& cells = this->getTraversal(TraversalManager::REAL);

  Array<uint> start(cells.size(), 0);

  // set start for each cell
  for (uint i=1, size=cells.size(); i<size; ++i) {
    start[i] = start[i-1] + particles[ cells[i-1] ].getNumberOfParticles();
  }

  parallel_region(0, cells.size(), [&](const uint begin, const uint end) {

      for(uint i=begin; i<end; ++i){
	particles[ cells[i] ].fillBuffer(buffer, start[i]);
      }

    });
  
}


/// @brief Fill a buffer with all the cells of the node, case of a CellOutput buffer
/// @brief [in] Buffer to fill
TMPLSG void TMPL_AMRGrid::fillBuffer(CellOutput* buffer) {

  const auto& cells = this->getTraversal(TraversalManager::REAL);

  parallel_region(0, cells.size(), [&](const uint begin, const uint end) {

      for(uint i=begin; i<end; ++i) {
      	vec3<int> cellCoord = this->coords[ cells[i] ];
      	uint idx = cellCoord.x + cellCoord.y*Global::domainInfo.getNumberOfCellsPerDim().x + cellCoord.z*Global::domainInfo.getNumberOfCellsPerDim().x*Global::domainInfo.getNumberOfCellsPerDim().y;

      	// Write cell index
      	buffer->id[i] = idx;
      	// Convert and write cell positions
      	buffer->r[i] = Global::domainInfo.getMinBounds() + cellCoord*Global::domainInfo.getCellLength();
      	buffer->r[i].x = convert(buffer->r[i].x,Stamp_Units::length,SI_Units_base::meter);
      	buffer->r[i].y = convert(buffer->r[i].y,Stamp_Units::length,SI_Units_base::meter);
      	buffer->r[i].z = convert(buffer->r[i].z,Stamp_Units::length,SI_Units_base::meter);
      	particles[ cells[i] ].fillBuffer(buffer, i);	
      }

    });

}

/// @brief Resize local buffers
///
///
TMPLSG inline void TMPL_AMRGrid::checkLocalBuffers() {
 
	// Parallelize the work by distributing the cells between the threads
  #pragma omp parallel for
	for(uint i=0; i<leafcells.size(); ++i)
	{
		leafcells[i]->checkLocals(
		  leafcells[i]->getNumberOfParticles()
		);
	}
}


/// @brief Resize local buffers
///
///
TMPLSG inline void TMPL_AMRGrid::checkLocalBuffersGHOST() {
 
	const auto& cells = this->getTraversal(TraversalManager::GHOST);

  #pragma omp parallel for
	for(uint i=0; i<cells.size(); ++i){
		particles[cells[i]].checkLocals(
		  particles[cells[i]].getNumberOfParticles()
		);
  }
}
