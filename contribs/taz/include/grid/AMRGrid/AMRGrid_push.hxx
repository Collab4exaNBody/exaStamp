#pragma once

/// @brief Update particles positions using to first order approximation of the movement equations (r += time * v)
/// @param [in] time Time step
TMPLSG inline void TMPL_AMRGrid::pushPositions1stOrder(double time) {

	constexpr size_t order=1;

  #pragma omp parallel for
  for( int i = 0 ; i < octrees.size() ; i++)
			pushNOrderHelperAssert<order>::pushNOrderAssert(time,
				size_t (octrees[i]->size),
				octrees[i]->getVelocityX(), octrees[i]->getPositionX(), //X
				octrees[i]->getVelocityY(), octrees[i]->getPositionY(), //Y
				octrees[i]->getVelocityZ(), octrees[i]->getPositionZ()  //Z
			);
}


/// @brief Update particles positions using to second order approximation of the movement equations (r += time * v + 0.5*time*time * f)
/// @param [in] time Time step
TMPLSG inline void TMPL_AMRGrid::pushPositions2ndOrder(double time) {

	constexpr size_t order=2;
	
  #pragma omp parallel for
  for( int i = 0 ; i < octrees.size() ; i++)
			pushNOrderHelperAssert<order>::pushNOrderAssert(
				time,
				size_t(octrees[i]->size),
				octrees[i]->getForceX(), octrees[i]->getVelocityX(), octrees[i]->getPositionX(),//X
				octrees[i]->getForceY(), octrees[i]->getVelocityY(), octrees[i]->getPositionY(),//Y
				octrees[i]->getForceZ(), octrees[i]->getVelocityZ(), octrees[i]->getPositionZ() //Z
			);
}


/// @brief Update particles velocities using to first order approximation of the movement equations (v += time * f)
/// @param [in] time Time step
TMPLSG inline void TMPL_AMRGrid::pushVelocities1stOrder(double time) {

	constexpr size_t order=1;

  #pragma omp parallel for
  for( int i = 0 ; i < octrees.size() ; i++)
				pushNOrderHelperAssert<order>::pushNOrderAssert(time,
				size_t (octrees[i]->size),
				octrees[i]->fx.data(), octrees[i]->vx.data(),//getVelocityX(), //X
				octrees[i]->fy.data(), octrees[i]->vy.data(),//getVelocityY(), //Y
				octrees[i]->fz.data(), octrees[i]->vz.data()//getVelocityZ()  //Z
			);
}

/// @brief Update particles velocities by applying Langevin dissipation
/// @param [in] time Time step
/// @param [in] gamma Gamma parameter of the Langevin scheme
/// @param [in] beta Beta parameter of the Langevin scheme
TMPLSG inline void TMPL_AMRGrid::pushDissipationLangevin(double time, double gamma, double beta) {

	static uint iter = 0;
	
	// Calculate arguments of the computation
  const double minus_gamma_dt = -1.00 * gamma * time;
	
	#pragma omp parallel for schedule(runtime)
	for(int i = 0 ; i < octrees.size() ; i++)
	{
	    uint8_t* ptr_ti = octrees[i]->getType();
	    uint64_t* ptr_id = octrees[i]->getId();
	    double* ptr_vx = octrees[i]->getVelocityX();
	    double* ptr_vy = octrees[i]->getVelocityY();
	    double* ptr_vz = octrees[i]->getVelocityZ();  
	    
	    size_t numberOfParticles = octrees[i]->getNumberOfParticles();
	    
	       // Fill an array with inverse masses
	    std::vector<double> invMass(numberOfParticles);
	    for (uint i=0; i<numberOfParticles; ++i) 
	      invMass[i] = Global::reference.getInvMass(ptr_ti[i]);

	    // Fill two arrays with random doubles
	    std::vector<double> rnd1   (numberOfParticles);
	    std::vector<double> rnd2   (numberOfParticles);

	    for (uint i=0; i<numberOfParticles; ++i) {
	      rnd1[i] = (double) Saru::saru((ptr_id[i]<<2),   Saru::int32_max, 3*(iter+Global::seed)+0);
	      rnd2[i] = (double) Saru::saru((ptr_id[i]<<2)+1, Saru::int32_max, 3*(iter+Global::seed)+0);
	    }
	    // Call vectorized computation for the x component
	    simdAMR::kernels::pushFluctuationLangevin(minus_gamma_dt, ptr_vx, beta, invMass.data(), rnd1.data(), rnd2.data(), numberOfParticles);

	    // Get new random doubles
	    for (uint i=0; i<numberOfParticles; ++i) {
	      rnd1[i] = (double) Saru::saru((ptr_id[i]<<2),   Saru::int32_max, 3*(iter+Global::seed)+1);
	      rnd2[i] = (double) Saru::saru((ptr_id[i]<<2)+1, Saru::int32_max, 3*(iter+Global::seed)+1);
	    }
	    // Call vectorized computation for the y component
	    simdAMR::kernels::pushFluctuationLangevin(minus_gamma_dt, ptr_vy, beta, invMass.data(), rnd1.data(), rnd2.data(), numberOfParticles);

	    // Get new random doubles
	    for (uint i=0; i<numberOfParticles; ++i) {
	      rnd1[i] = (double) Saru::saru((ptr_id[i]<<2),   Saru::int32_max, 3*(iter+Global::seed)+2);
	      rnd2[i] = (double) Saru::saru((ptr_id[i]<<2)+1, Saru::int32_max, 3*(iter+Global::seed)+2);
	    }
	    // Call vectorized computation for the z component
	    simdAMR::kernels::pushFluctuationLangevin(minus_gamma_dt, ptr_vz, beta, invMass.data(), rnd1.data(), rnd2.data(), numberOfParticles);
 
	    #pragma omp atomic	
	    ++iter;   
	}
	


}
