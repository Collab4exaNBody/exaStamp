#pragma once

/// @brief Get the particles and initialize those in the grid
/// Check init type and call a specialized function
/// @param [in] particleInit Configuration of the particles
/// @param [in] begin Start of the set
/// @param [in] end End of the set
TMPLSG void TMPL_AMRGrid::initParticles(Configuration<Particle>& particleInit, uint64_t begin, uint64_t end) {

	switch (particleInit.getInitType()) {

		case InitType::INIT_DEFAULT :
		initParticlesHard(particleInit, begin, end);
		break;

		case InitType::INIT_STAMP_LEGACY_DUMP :
		initParticlesFile(particleInit, begin, end);
		break;

		case InitType::INIT_HERCULE_DUMP :
		#ifdef __use_lib_hercule
		initHerculeParticlesFile(particleInit, begin, end);
		#endif

		break;
		case InitType::INIT_STAMPV4 :
		std::cerr << "Error ! Particles initialization from Configuration<Particule> should not occur in a molecular initialization." << std::endl;
	}
}



/// @brief Initialize the velocity of a particle
///
/// This should be completed with shock initialization and other strange velocities initilizations
/// @param [in,out] p Pointer to the particle
/// @param [in] sigma Sigma parameters for all the particles
TMPLSG void TMPL_AMRGrid::initVelocities(P& p, Array<double>& sigma) {
	p.velocity()=Global::reference.getVelocity(p.type());
	if (Global::reference.getInvMass(p.type()) != 0) 
		p.velocity() += sigma[p.type()] * vec3<double>(Saru::randN(p.index(), -1), Saru::randN(p.index(), -2), Saru::randN(p.index(), -3));
}


/// @brief Get a chunk of particles from the configuration and initialize those that are in the grid (lattice case)
/// @param [in] particleInit Configuration of the particles
/// @param [in] cell_begin Start of the set
/// @param [in] cell_end End of the set
TMPLSG void TMPL_AMRGrid::initParticlesHard(Configuration<Particle>& particleInit, uint64_t cell_begin, uint64_t cell_end) {

	vec3<int> NCELLS = particleInit.NCELLS;
	double A = particleInit.lattice->parameters[0];

	// Get the boundary conditions
	auto BC = Global::domainInfo.getBoundaryConditions();

	const uint64_t nyz = (uint64_t) NCELLS.y * (uint64_t) NCELLS.z;
	const uint64_t ny  = (uint64_t) NCELLS.y;
	const uint64_t nz  = (uint64_t) NCELLS.z;

		// Function to convert the index of a cell into a coordinate
	auto index_convert = [&] (const uint64_t& c) -> vec3<double> 
	{
		uint64_t x = c/nyz;
		uint64_t y = (c/nz) % ny;
		uint64_t z = c%nz;
		return vec3<double>( (double) x, (double) y, (double) z);
	};

	vec3<double> center = Global::domainInfo.getMinBounds() + 0.5*Global::domainInfo.getExtension() - 0.5*NCELLS*A;


	vec3<double> minAtom = center + A*index_convert(cell_begin);
	vec3<double> maxAtom = center + A*(index_convert(cell_end))+A;  
	uint64_t xbeg = cell_begin/nyz;
	uint64_t xend = cell_end/nyz;
	if(minAtom.x <= info.getLimitSup().x && maxAtom.x >= info.getLimitInf().x) {}
	else return;

	const vec3<double> invCellLength = 1./Global::domainInfo.getCellLength();
	const vec3<double>& minBounds = Global::domainInfo.getMinBounds();
	spin_mutex mutex_np;

	#ifdef __use_lib_omp
	uint64_t grain = std::max ( uint64_t(double(cell_end-cell_begin)/double(256*4)),uint64_t(150));
	#endif

	// Parallelize the work by distributing the chunk between the threads
	parallel_region(cell_begin, cell_end,

		#ifdef __use_lib_omp
		grain,
		#endif
		[&](const uint64_t begin, const uint64_t end) 
		{

			// get some particles to init into tmp
			std::vector<P> tmp;
			particleInit.buildSubset(begin, end, tmp, info.getLimitInf(), info.getLimitSup());

			// remove some particles
			for (auto it=tmp.rbegin(); it!=tmp.rend(); ++it) {
				if (!info.inGrid((*it).r)) {
					*it = tmp.back();
					tmp.pop_back();
				}
			}

			// add remaining particles
			addParticles(tmp.data(), tmp.size(), false, [&] (const P& p) -> vec3<int> {return auxFloor<int>( invCellLength * (p.r-minBounds) );});

			// update num particles
			mutex_np.lock();
			this->numberOfParticles += (uint64_t) tmp.size();
			mutex_np.unlock();

		}
	);
}


/// @brief Get a chunk of particles from the configuration and initialize those that are in the grid
///(case of an initialization from a Stamp file)
/// @param [in] particleInit Configuration of the particles
/// @param [in] begin Start of the set
/// @param [in] end End of the set
TMPLSG void TMPL_AMRGrid::initParticlesFile(Configuration<Particle>& particleInit, uint64_t begin,  uint64_t end) {


	const vec3<double> invCellLength = 1./Global::domainInfo.getCellLength();
	const vec3<double>& minBounds     = Global::domainInfo.getMinBounds();

	// array to import particles
	Array<MPI__Particle> tmp(end-begin);

	// read file
	if (this->getDomainIndex()==0) {    
		InputOutputManager::Offset offset = sizeof(LegacyHeaderIOStruct) + begin * sizeof(LegacyParticleIOStruct);
		readLegacyParticles(particleInit.getFileId(), offset, tmp.size(), tmp.data());
	}

	this->domain->getCommManager()->broadcast(tmp);

	size_t tmpSize = tmp.size();

	#ifdef __use_lib_omp
	uint64_t grain = std::max<uint64_t> ( uint64_t(double(tmpSize)/double(256*4)),uint64_t(150));
	#endif

	parallel_region(0, tmpSize,
	#ifdef __use_lib_omp
	grain,
	#endif
	[&](const size_t begin, const size_t end) 
		{
			uint64_t countAddParticles=0;

			for(size_t it = begin ; it < end ; it++)
			{
				int index = info.inGridAndGetIndex( invCellLength * (tmp[it].r-minBounds) );
				if(index!=-1)
				{
					assert(tmp[it].r.x - minBounds.x >= (particles[index].position.x-1)*Global::domainInfo.getCellLength().x);
					assert(tmp[it].r.y - minBounds.y >= (particles[index].position.y-1)*Global::domainInfo.getCellLength().y);
					assert(tmp[it].r.z - minBounds.z >= (particles[index].position.z-1)*Global::domainInfo.getCellLength().z);

					assert(tmp[it].r.x - minBounds.x <= (particles[index].position.x)*Global::domainInfo.getCellLength().x);
					assert(tmp[it].r.y - minBounds.y <= (particles[index].position.y)*Global::domainInfo.getCellLength().y);
					assert(tmp[it].r.z - minBounds.z <= (particles[index].position.z)*Global::domainInfo.getCellLength().z);
					this->lock(index);
					particles[index].add(tmp[it]);
					this->unlock(index);
					countAddParticles++;
				}
			}

			this->lock(0);
			this->numberOfParticles += countAddParticles;
			this->unlock(0); 
		}
	); 
}

#ifdef __use_lib_hercule
/// @brief Get a chunk of particles from the configuration and initialize those that are in the grid
///(case of an initialization from a Hercule file)
/// @param [in] particleInit Configuration of the particles
/// @param [in] numSSDomBegin Start of the set
/// @param [in] numSSDomEnd End of the set
TMPLSG void TMPL_AMRGrid::initHerculeParticlesFile(Configuration<Particle>& particleInit, uint64_t numSSDomBegin, uint64_t numSSDomEnd) {

  const auto invCellLength = 1./Global::domainInfo.getCellLength();
  const auto& minBounds = Global::domainInfo.getMinBounds();

  spin_mutex mutex_np;

  // open ssdom
  HerculeParticleIODumpStruct ioDumpReadStruct;

  for(u_int numSSDom = numSSDomBegin ; numSSDom<numSSDomEnd;numSSDom++){
    this->m_inputOutputManager->herculeDumpR_readParticles(numSSDom,ioDumpReadStruct);

    {

      // Get some particles to init into tmp
      //std::vector<MPI__Particle> tmp;
      Array<MPI__Particle> tmp(ioDumpReadStruct.nbParticle);
      uint tmpSize=0;
      for(int i=0;i<ioDumpReadStruct.nbParticle;i++){
      	vec3<double> r(ioDumpReadStruct.particles_coordinates[0][i],ioDumpReadStruct.particles_coordinates[1][i],ioDumpReadStruct.particles_coordinates[2][i]);

      	// Add the particles of the grid in a temporary array
      	if (info.inGrid(r)) {
      		MPI__Particle p;
      		p.r.x = ioDumpReadStruct.particles_coordinates[0][i];
      		p.r.y = ioDumpReadStruct.particles_coordinates[1][i];
      		p.r.z = ioDumpReadStruct.particles_coordinates[2][i];
      		p.id =ioDumpReadStruct.particles_iD[i]-1;
      		p.v.x=ioDumpReadStruct.particles_velocity[0][i];
      		p.v.y=ioDumpReadStruct.particles_velocity[1][i];
      		p.v.z=ioDumpReadStruct.particles_velocity[2][i];
      		tmp[tmpSize++]=p;
      	}

      }

      // Add the particles to the grid
      addParticles(tmp.data(), tmpSize, false, [&] (const MPI__Particle& p) -> vec3<int> {return auxFloor<int>( invCellLength * (p.r-minBounds) );});

      // Update the number of particles
      mutex_np.lock();
      this->numberOfParticles += (uint64_t) tmpSize;
      mutex_np.unlock();
    }
  }


}
#endif
