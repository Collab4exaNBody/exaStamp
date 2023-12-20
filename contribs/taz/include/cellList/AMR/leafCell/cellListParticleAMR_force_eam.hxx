#pragma once

/// @brief Function to compute rho terms for a EAM potential with SIMD optimized version --> Modify computeForceRho
/// between particles of two specified types
/// @tparam Pot_t Type of the potential
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
/// @param [in] cells All the cells in the grid
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] correcter Tool to correct distances
/// @param [in] symmetrize Indicates if the force computation is symmetrized
template <class Pot_t>
void leafCell::computeForceRhoVerlet_without_buffer(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, const bool symmetrize) {

  NeighborList_base::nbr_id* start;
  uint n;  
  double rcut2 = Global::reference.getRcut2(typeIndexA, typeIndexB);

  // Loop on the particles of the cell
  for (uint i=0; i<size; ++i) {

  	// Get the complementary type to the type of the particle
    uint8_t otherType;
    if      (getType(i) == typeIndexA) otherType = typeIndexB;
    else if (getType(i) == typeIndexB) otherType = typeIndexA;
    else continue;
    
    neighborListGetNeighbors(i, otherType, start, n);

    int shifti = shift+i;

    double rhoCum=0;

    const double rxi = granny->rx[shifti];
    const double ryi = granny->ry[shifti];
    const double rzi = granny->rz[shifti];

    leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(shifti);

    // Loop on the neighbors
    for (size_t j=0; j<n; ++j) {

      auto nbr  = neighborCellsList[std::get<0>(start[j])];
      int idx   = std::get<1>(start[j]) + nbr->shift;

      assert(nbr->getIsLeaf() == true );

      // compute dist (xi, xj)
      double rxij = rxi - nbr->granny->rx[idx];    
      double ryij = ryi - nbr->granny->ry[idx];    
      double rzij = rzi - nbr->granny->rz[idx];  

      const double r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      
      assert(r2 != 0);

      if(r2>rcut2) continue;

      double rhotmp=0;

      double r = std::sqrt(r2);
      double ratio = pot->p.a0/r;
      rhotmp  = std::pow(ratio, pot->p.m);
      rhotmp -= pot->rhoCut;


      rhoCum += rhotmp;

      if(symmetrize)
      {
	if(!nbr->getIsGhost())
	{
		auto tmpId = std::get<1>(start[j]);
	      	nbr->rhoAMR(tmpId) += rhotmp;
	}
      }
    }

    this->rhoAMR(i) += rhoCum;    
  }
	    
}


/// @brief Function to compute rho terms for a EAM potential with SIMD optimized version --> Modify computeForceRho
/// between particles of two specified types
/// @tparam Pot_t Type of the potential
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
/// @param [in] cells All the cells in the grid
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] correcter Tool to correct distances
/// @param [in] symmetrize Indicates if the force computation is symmetrized
template <class Pot_t>
void leafCell::computeForceRhoVerlet_without_buffer_per_block(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  NeighborList_base::nbr_id* start;
  uint n;  
  double rcut2 = Global::reference.getRcut2(typeIndexA, typeIndexB);

  double   rxij, ryij, rzij;
  ssize_t idx;

  // Loop on the particles of the cell
  for (ssize_t i=0; i<this->size; ++i) {

  	// Get the complementary type to the type of the particle
    uint8_t otherType;
    if      (getType(i) == typeIndexA) otherType = typeIndexB;
    else if (getType(i) == typeIndexB) otherType = typeIndexA;
    else continue;
    
    neighborListGetNeighbors(i, otherType, start, n);

    int shifti = shift+i;

    double rhoCum=0;

    const double rxi = granny->rx[shifti];
    const double ryi = granny->ry[shifti];
    const double rzi = granny->rz[shifti];

    leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(shifti);

    // Loop on the neighbors
    for (size_t j=0; j<n; ++j) {

      auto* nbr  = neighborCellsList[std::get<0>(start[j])];
      ssize_t idxs   = std::get<1>(start[j]) + nbr->shift;
      ssize_t idxNoShift = std::get<1>(start[j]);
      ssize_t end = nbr->shift + nbr->size;//idxs+SIZEVECTOR < nbr->getNumberOfParticles() ? idxs+SIZEVECTOR :  nbr->getNumberOfParticles() ;

      double *rho_ptr = &(nbr->rhoAMR(std::get<1>(start[j])));

      #pragma omp simd  reduction(+:rhoCum)
      #pragma vector aligned
      for(ssize_t s = 0; s < SIZEVECTOR; s++)
      {
        idx = idxs + s;

	if(idx <end)
	{

		rxij = rxi - nbr->granny->rx[idx];    
		ryij = ryi - nbr->granny->ry[idx];    
		rzij = rzi - nbr->granny->rz[idx];  


		const double r2 = rxij*rxij + ryij*ryij + rzij*rzij;

		
		if(r2<rcut2 )
		{
		      double dRho = 0.;
		      double rhotmp=0;
		      //assert(r2!=0);
		      double r = std::sqrt(r2);
		      double ratio = pot->p.a0/r;
		      rhotmp  = std::pow(ratio, pot->p.m);
		      rhotmp -= pot->rhoCut;

		      rhoCum += rhotmp;

			if(!nbr->getIsGhost())
			{
			      	rho_ptr[s] += rhotmp;
			}
		 }
	       }
	}
    }
    this->rhoAMR(i) += rhoCum;    

  }
	    
}

// @brief Function to compute forces for a EAM potential with SIMD optimized version. Don't forget to modify Verlet version.
/// between particles of two specified types
/// @tparam Pot_t Type of the potential
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
/// @param [in] cells All the cells in the grid
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] correcter Tool to correct distances
/// @param [in] symmetrize Indicates if the force computation is symmetrized
template <class Pot_t>
void leafCell::computeForceFinalVerlet_without_buffer(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, const bool symmetrize) {

  NeighborList_base::nbr_id* start;
  uint n;  
  double rcut2 = Global::reference.getRcut2(typeIndexA, typeIndexB);

  // Loop on the particles of the cell
  for (uint i=0; i<size; ++i) {

  	// Get the complementary type to the type of the particle
    uint8_t otherType;
    if      (getType(i) == typeIndexA) otherType = typeIndexB;
    else if (getType(i) == typeIndexB) otherType = typeIndexA;
    else continue;

    const double pim = Global::reference.getInvMass(getType(i));
    const double nim = -Global::reference.getInvMass(otherType);
    
    neighborListGetNeighbors(i, otherType, start, n);

    int shifti = shift+i;

    const double rxi = granny->rx[shifti];
    const double ryi = granny->ry[shifti];
    const double rzi = granny->rz[shifti];

    leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(shifti);

    double dex(0.),dey(0.),dez(0.), dep(0.);

    double emb=this->embAMR(i);

    // Loop on the neighbors
    for (size_t j=0; j<n; ++j) {

      auto nbr  = neighborCellsList[std::get<0>(start[j])];
      int idx   = std::get<1>(start[j]) + nbr->shift;
      int idxNoShift = std::get<1>(start[j]);

      assert(nbr->getIsLeaf() == true );

      // compute dist (xi, xj)
      double rxij = rxi - nbr->granny->rx[idx];    
      double ryij = ryi - nbr->granny->ry[idx];    
      double rzij = rzi - nbr->granny->rz[idx];  

      const double r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      
      assert(r2 != 0);

      if(r2>rcut2) continue;

      double phi  = 0.; 
      double dPhi = 0.;
      double rho  = 0.; 
      double dRho = 0.;

      const double inv_r = 1./std::sqrt(r2);

      const double ratio = pot->p.a0*inv_r;
      rho = std::pow(ratio, pot->p.m);
      dRho = -1 * pot->p.m * rho * inv_r;

      phi  = pot->p.epsilon * std::pow(ratio, pot->p.n);
      dPhi = -1 * pot->p.n * phi * inv_r;
      phi -= pot->phiCut;


      // Get the force and energy
      double de = -1. *(dRho * (emb + nbr->embAMR(idxNoShift)) + dPhi)*inv_r;



      dex += de * rxij;
      dey += de * ryij;
      dez += de * rzij;
      dep += .5 * phi;

      if(symmetrize)
      {
	      de *= nim;
	      nbr->granny->fx[idx] += de * rxij;
	      nbr->granny->fy[idx] += de * ryij;
	      nbr->granny->fz[idx] += de * rzij;
	      nbr->granny->ep[idx] += .5 * phi;	
      }
    }

    granny->fx[shifti] += dex *pim;
    granny->fy[shifti] += dey *pim;
    granny->fz[shifti] += dez *pim;
    granny->ep[shifti] += dep;

  }
	    
}


// @brief Function to compute forces for a EAM potential with SIMD optimized version. Don't forget to modify Verlet version.
/// between particles of two specified types
/// @tparam Pot_t Type of the potential
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
/// @param [in] cells All the cells in the grid
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] correcter Tool to correct distances
/// @param [in] symmetrize Indicates if the force computation is symmetrized
template <class Pot_t>
void leafCell::computeForceFinalVerlet_without_buffer_per_block(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  NeighborList_base::nbr_id* start;
  uint n;  
  double rcut2 = Global::reference.getRcut2(typeIndexA, typeIndexB);
  size_t end,idx;

  double   rxij, ryij, rzij;

  // Loop on the particles of the cell
  for (ssize_t i=0; i<size; ++i) {

  	// Get the complementary type to the type of the particle
    uint8_t otherType;
    if      (getType(i) == typeIndexA) otherType = typeIndexB;
    else if (getType(i) == typeIndexB) otherType = typeIndexA;
    else continue;

    const double pim = Global::reference.getInvMass(getType(i));
    const double nim = -Global::reference.getInvMass(otherType);
    
    neighborListGetNeighbors(i, otherType, start, n);

    int shifti = shift+i;

    const double rxi = granny->rx[shifti];
    const double ryi = granny->ry[shifti];
    const double rzi = granny->rz[shifti];

    leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(shifti);

    double dex(0.),dey(0.),dez(0.), dep(0.);

    double emb=this->embAMR(i);

    // Loop on the neighbors
    for (ssize_t j=0; j<n; ++j) {

      auto* nbr  = neighborCellsList[std::get<0>(start[j])];
      size_t idxs   = std::get<1>(start[j]) + nbr->shift;

      int idxNoShift = std::get<1>(start[j]);
      
      end = nbr->shift + nbr->getNumberOfParticles();//idxs+SIZEVECTOR < nbr->getNumberOfParticles() ? idxs+SIZEVECTOR :  nbr->getNumberOfParticles() ;

      double *ptrEmbAmr = &(nbr->embAMR(idxNoShift));      

      #pragma omp simd reduction(+:dep,dex,dey,dez)
      #pragma vector aligned
      for(ssize_t s = 0; s < SIZEVECTOR; s++)
      {
        idx = idxs + s;

        rxij = rxi - nbr->granny->rx[idx];    
        ryij = ryi - nbr->granny->ry[idx];    
        rzij = rzi - nbr->granny->rz[idx];  


        const double r2 = rxij*rxij + ryij*ryij + rzij*rzij;

        
        if(r2<rcut2 && idx <end )
        {
	      const double inv_r    = 1./std::sqrt(r2);
	      double phi  = 0.; 
	      double dPhi = 0.;
	      double dRho = 0.;
	      
	      // Apply the potential the get the rho and phi terms

	      const double ratio = pot->p.a0*inv_r;
              const double rho = std::pow(ratio, pot->p.m);
	      dRho = -1 * pot->p.m * rho * inv_r;

              phi  = pot->p.epsilon * std::pow(ratio, pot->p.n);
              dPhi = -1 * pot->p.n * phi * inv_r;
              phi -= pot->phiCut;


	      // Get the force and energy
	      double de = -1. *(dRho * (emb + ptrEmbAmr[s]) + dPhi)*inv_r;

	      dex += de * rxij;
	      dey += de * ryij;
	      dez += de * rzij;
	      dep += .5 * phi;

	      de *= nim;
	      nbr->granny->fx[idx] += de * rxij;
	      nbr->granny->fy[idx] += de * ryij;
	      nbr->granny->fz[idx] += de * rzij;
	      nbr->granny->ep[idx] += .5 * phi;	
	 }
       }

    }
    granny->fx[shifti] += dex *pim;
    granny->fy[shifti] += dey *pim;
    granny->fz[shifti] += dez *pim;
    granny->ep[shifti] += dep;
  }
}
