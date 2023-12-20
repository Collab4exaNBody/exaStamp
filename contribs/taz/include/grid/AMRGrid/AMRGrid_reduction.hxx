#pragma once

/// @brief Compute total energy of the particles of the grid
///
///
TMPLSG inline void TMPL_AMRGrid::computeEnergy() {

	this->energyKin = 0.;
	this->energyPot = 0.;

	// Parallelize the work by distributing the cells between the threads
	double tmpKin(0.), tmpPot(0.);

	#pragma omp parallel for reduction(+:tmpKin, tmpPot)
	for(int o = 0; o < octrees.size(); ++o)
	{
    double eKin = 0.;
    double ePot = 0.;
 
    uint8_t* ptr_ti = octrees[o]->getType();
    double* ptr_vx = octrees[o]->getVelocityX();
    double* ptr_vy = octrees[o]->getVelocityY();
    double* ptr_vz = octrees[o]->getVelocityZ();	
    double* ptr_ep = octrees[o]->getPotentialEnergy();
    
    // Sum on all the particles
    #pragma omp simd reduction(+:eKin,ePot)
    for (uint i=0; i< octrees[o]->getNumberOfParticles(); ++i) 
    {
      // avoid NaN
      assert(Global::reference.getInvMass(ptr_ti[i]) != 0) ; 
      eKin += (ptr_vx[i]*ptr_vx[i] + ptr_vy[i]*ptr_vy[i] + ptr_vz[i]*ptr_vz[i]) * Global::reference.getMass(ptr_ti[i]);
      ePot += ptr_ep[i];
    }
    
	
		tmpKin += 0.5*eKin;
		tmpPot += ePot;
	}

	this->energyKin += tmpKin;
	this->energyPot += tmpPot;
}



/// @brief Compute workload of the particles of the grid
///
///
TMPLSG inline void TMPL_AMRGrid::computeWorkload() {

	this->m_workload = 0.;

	const auto& cells = this->getTraversal(TraversalManager::REAL);
	
	double w(0.);
	
	#pragma omp parallel for reduction(+:w)
	for(int i = 0 ; i < cells.size() ; i++)
	{
	  w += particles[cells[i]].computeWorkloadAMR();
	}

  this->m_workload += w;
}


/// @brief Compute total momentum of the particles of the grid
/// @return Total momentum
TMPLSG inline vec3<double> TMPL_AMRGrid::computeTotalMomentum() {

  
  double totalMomentum_x(0.), totalMomentum_y(0.), totalMomentum_z(0.);
  
  #pragma omp parallel for reduction(+:totalMomentum_x,totalMomentum_y,totalMomentum_z)
	for(int o = 0; o < octrees.size(); ++o)
	{
	  vec3<double> localMomentum = 0.;
    uint8_t* ptr_ti = octrees[o]->getType();
    double* ptr_vx = octrees[o]->getVelocityX();
    double* ptr_vy = octrees[o]->getVelocityY();
    double* ptr_vz = octrees[o]->getVelocityZ();		

    for (uint i=0; i< octrees[o]->getNumberOfParticles(); ++i) 
    {
      if (Global::reference.getInvMass(ptr_ti[i]) != 0)
        localMomentum += Global::reference.getMass(ptr_ti[i]) * vec3<double>(ptr_vx[i], ptr_vy[i], ptr_vz[i]);
    }
    
    totalMomentum_x += localMomentum.x;
    totalMomentum_y += localMomentum.y;
    totalMomentum_z += localMomentum.z;        
	}
   return vec3<double>(totalMomentum_x, totalMomentum_y, totalMomentum_z);
}


/// @brief Compute total mass of the particles of the grid
/// @return Total mass
TMPLSG inline double TMPL_AMRGrid::computeTotalMass() {

  double totalMass(0.);
  
  #pragma omp parallel for reduction(+:totalMass)
	for(size_t o = 0; o < octrees.size(); ++o)
	{
	  double localMass = 0.;
    uint8_t* ptr_ti = octrees[o]->getType();

    for (size_t i=0 ; i< octrees[o]->getNumberOfParticles(); ++i) 
    {
      if (Global::reference.getInvMass(ptr_ti[i]) != 0)
      {
        localMass += Global::reference.getMass(ptr_ti[i]);
      }
    }
    totalMass += localMass;
	}
  
  return totalMass;
}

/// @brief Get the total kinetic energy in the center of momentum frame
/// @param [in] vshift Global velocity of the system
/// @return Total kinetic energy in the center of momentum frame
TMPLSG inline double TMPL_AMRGrid::computeShiftedKineticEnergy(const vec3<double>& vshift) {

	// Parallelize the work by distributing the cells between the threads
	double totalShiftedKE(0.);
	
	#pragma omp parallel for reduction(+:totalShiftedKE)
	for(int i = 0; i < octrees.size(); ++i)
	{
    double localShiftedKE(0.);
 
    uint8_t* ptr_ti = octrees[i]->getType();
    double* ptr_vx = octrees[i]->getVelocityX();
    double* ptr_vy = octrees[i]->getVelocityY();
    double* ptr_vz = octrees[i]->getVelocityZ();	
    double* ptr_ep = octrees[i]->getPotentialEnergy();
    
    // Sum on all the particles
    #pragma omp simd reduction(+:localShiftedKE)
    for (size_t j=0; j< octrees[i]->getNumberOfParticles(); ++j) 
    {
      // avoid NaN
      assert(Global::reference.getInvMass(ptr_ti[j]) != 0) ; 
      localShiftedKE += norm2(vec3<double>(ptr_vx[j]-vshift.x, ptr_vy[j]-vshift.y, ptr_vz[j]-vshift.z)) * Global::reference.getMass(ptr_ti[j]);
    }
		totalShiftedKE += localShiftedKE;
	}
  
  totalShiftedKE *= 0.5;

  return totalShiftedKE;
}


/// @brief Get the pressure tensor
/// @param [in] vshift Global velocity of the system
/// @return Total pressure tensor
TMPLSG inline mat3<double> TMPL_AMRGrid::computePressure(const vec3<double>& vshift) {


	// Parallelize the work by distributing the cells between the threads
	mat3<double> totalPressure(0.);
	
	#pragma omp parallel for
	for(size_t i = 0; i < octrees.size(); ++i)
	{
    mat3<double> localPressure(0.);
 
    uint8_t* ptr_ti = octrees[i]->getType();
    double* ptr_vx = octrees[i]->getVelocityX();
    double* ptr_vy = octrees[i]->getVelocityY();
    double* ptr_vz = octrees[i]->getVelocityZ();	
    double* ptr_ep = octrees[i]->getPotentialEnergy();
    
    localPressure = 0.5 * octrees[i]->pressureTensor;
    
    // Sum on all the particles
    for (size_t j=0; j< octrees[i]->getNumberOfParticles(); ++j) 
    {
      // avoid NaN
      assert(Global::reference.getInvMass(ptr_ti[j]) != 0) ; 
      
      vec3<double> velocity = vec3<double>(ptr_vx[j], ptr_vy[j], ptr_vz[j]) - vshift;
      localPressure += tensor(velocity,velocity) * Global::reference.getMass(ptr_ti[j]) ;
    }
    
    this->lock(0);
		totalPressure += localPressure;
		this->unlock(0);
	}
  return totalPressure;
}

/// @brief Determines whether an atom has moved more than 1/2 of the verlet radius
/// @return boolean value
TMPLSG bool TMPL_AMRGrid::checkVerlet() {

  // check verlet list 
  double rVerlet =  Global::reference.getrVerlet()*0.5;
  rVerlet *= rVerlet;

  int value = 0;
  
  #pragma omp parallel for reduction(max:value) 
  for(size_t i = 0; i < octrees.size() ; ++i)
  {
    const int _bool = octrees[i]->checkVerlet(rVerlet);
    if(value < _bool)
      value = _bool;
  }  
  
  if(value >= 1)
    return true;
  else
    return false;
   
}
