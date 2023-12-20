#pragma once



/// @brief Function to compute forces for a pair potential with SIMD optimized version. Don't forget to modify the linked cell version.
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
template < class Pot_t>
void leafCell::computeForcePair_buffer(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  // No embedding term in a pair potential
  static constexpr bool fill_emb = false;

  // Get the SIMD optimized operator to compute the force
  typename Pot_t::simd_opt_t fastOp;
  pot->setParameters(fastOp);
  
  // Get vectorization buffer
  auto& fb = self_type::getVectBuffer(this->getMaxNumberOfNeighbors());	

  // Loop on the particles of the cell
  for (uint i=0; i<this->size; ++i) {

  	// Get the complementary type to the type of the particle
    uint8_t otherType;
    if      (getType(i) == typeIndexA) otherType = typeIndexB;
    else if (getType(i) == typeIndexB) otherType = typeIndexA;
    else continue;


    // Fill the buffer with atoms included of Verlet lists of that type
    uint nbrSize = fillForceBufferVerlet<fill_emb>(i, otherType);

    // Apply the potential to get the forces
    fastOp(fb.dfx(), fb.dfy(), fb.dfz(), fb.den(), fb.drx(), fb.dry(), fb.drz(), nbrSize);

    // Update particles values
    writeForce(i, otherType, true, nbrSize);
  }
	    
}



/// @brief Function to compute forces for a pair potential with SIMD optimized version
/// between particles of two specified types
/// @tparam Pot_t Type of the potential
/// @tparam symmetrize Indicates if the force computation is symmetrized
/// @tparam Verlet Indicates if the verlet lists are used
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
template <class Pot_t>
void leafCell::computeForcePair_perso(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  uint8_t *ptr_ti(getType());
  const auto& param = pot->getParameters();
                 
  const double   epsilon   ( param.epsilon ),
                 sigma6    ( param.sigma*param.sigma*param.sigma*param.sigma*param.sigma*param.sigma ),
                _2epsilon  (  2. * epsilon ),
                _24epsilon ( 24. * epsilon ),
                 ecut      ( 0.5 * pot->getEcut() );
                

  double   rxij, ryij, rzij,
            dex,  dey,  dez, 
              e, tmpe,tmpde, rcut2;
  
  uint8_t otherType;  
  NeighborList_base::nbr_id* start;  
  uint n;

  rcut2 = Global::reference.getRcut2(typeIndexA, typeIndexB);

  // Loop on the particles of the cell
  for (int i=0; i<size ; ++i) {

    // Get the complementary type to the type of the particle
    if      (ptr_ti[i] == typeIndexA) otherType = typeIndexB;
    else if (ptr_ti[i] == typeIndexB) otherType = typeIndexA;
    else continue;
    
    const double pim = Global::reference.getInvMass(ptr_ti[i]);
    const double nim = -Global::reference.getInvMass(otherType);

    neighborListGetNeighbors(i, otherType, start, n);

    e=0.; tmpe=0.; dex=0.; dey=0.; dez=0.;

    int shifti = shift+i;

    const double rxi = granny->rx[shifti];
    const double ryi = granny->ry[shifti];
    const double rzi = granny->rz[shifti];

    leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(shifti);
   
    for(int j=0 ; j<n ; ++j)
    {  
      auto nbr  = neighborCellsList[std::get<0>(start[j])];
      int idx   = std::get<1>(start[j]) + nbr->shift;

      assert(nbr->getIsLeaf() == true );

      // compute dist (xi, xj)
      rxij = rxi - nbr->granny->rx[idx];    
      ryij = ryi - nbr->granny->ry[idx];    
      rzij = rzi - nbr->granny->rz[idx];    

      const double r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      
      assert(r2 != 0);

      if(r2>rcut2) continue;
       
      const double ir2     = 1.0/r2;
      const double ratio6  = sigma6*ir2*ir2*ir2; 
      const double ratio12 = ratio6*ratio6;
      const double diff = ratio12 -ratio6;

      tmpe  = _2epsilon  * (diff)   - ecut;
      tmpde = _24epsilon * (diff+ratio12) * ir2;
      
      e   += tmpe;
      dex += rxij * tmpde;
      dey += ryij * tmpde;
      dez += rzij * tmpde;
      
      // third newton's law
	tmpde  *= nim;
	nbr->granny->fx[idx] += rxij*tmpde;
	nbr->granny->fy[idx] += ryij*tmpde;
	nbr->granny->fz[idx] += rzij*tmpde;
	nbr->granny->ep[idx] += tmpe;
    }

    granny->fx[shifti] += pim*dex;
    granny->fy[shifti] += pim*dey;
    granny->fz[shifti] += pim*dez;
    granny->ep[shifti] += e;
  }
	    
}


/// @brief Function to compute forces for a pair potential with SIMD optimized version
/// between particles of two specified types
/// @tparam Pot_t Type of the potential
/// @tparam symmetrize Indicates if the force computation is symmetrized
/// @tparam Verlet Indicates if the verlet lists are used
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
template < class Pot_t>
void leafCell::computeForcePairVerlet_mutex(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  uint8_t *ptr_ti(getType());
  const auto& param = pot->getParameters();
                 
  const double   epsilon   ( param.epsilon ),
                 sigma6    ( param.sigma*param.sigma*param.sigma*param.sigma*param.sigma*param.sigma ),
                _2epsilon  (  2. * epsilon ),
                _24epsilon ( 24. * epsilon ),
                 ecut      ( 0.5 * pot->getEcut() );
                

  double   rxij, ryij, rzij,
            dex,  dey,  dez, 
              e, tmpe,tmpde, rcut2;
  
  uint8_t otherType;  
  NeighborList_base::nbr_id* start;  
  uint n;

  rcut2 = Global::reference.getRcut2(typeIndexA, typeIndexB);

  // Loop on the particles of the cell
  for (int i=0; i<size ; ++i) {

    // Get the complementary type to the type of the particle
    if      (ptr_ti[i] == typeIndexA) otherType = typeIndexB;
    else if (ptr_ti[i] == typeIndexB) otherType = typeIndexA;
    else continue;
    
    const double pim = Global::reference.getInvMass(ptr_ti[i]);
    const double nim = -Global::reference.getInvMass(otherType);

    neighborListGetNeighbors(i, otherType, start, n);

    e=0.; tmpe=0.; dex=0.; dey=0.; dez=0.;

    int shifti = shift+i;

    const double rxi = granny->rx[shifti];
    const double ryi = granny->ry[shifti];
    const double rzi = granny->rz[shifti];

    leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(shifti);

    for(int j=0 ; j<n ; ++j)
    {  
      auto nbr  = neighborCellsList[std::get<0>(start[j])];
      int idx   = std::get<1>(start[j]) + nbr->shift;

      assert(nbr->getIsLeaf() == true );

      // compute dist (xi, xj)
      rxij = rxi - nbr->granny->rx[idx];    
      ryij = ryi - nbr->granny->ry[idx];    
      rzij = rzi - nbr->granny->rz[idx];    

      const double r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      
      assert(r2 != 0);

      if(r2>rcut2) continue;
       
      const double ir2     = 1.0/r2;
      const double ratio6  = sigma6*ir2*ir2*ir2; 
      const double ratio12 = ratio6*ratio6;
      const double diff = ratio12 -ratio6;

      tmpe  = _2epsilon  * (diff)   - ecut;
      tmpde = _24epsilon * (diff+ratio12) * ir2;
      
      e   += tmpe;
      dex += rxij * tmpde;
      dey += ryij * tmpde;
      dez += rzij * tmpde;
      
      // third newton's law
      tmpde  *= nim;
 
      if(!nbr->granny->getIsGhost())
      {
        nbr->granny->lock(idx);
        nbr->granny->fx[idx] += rxij*tmpde;
        nbr->granny->fy[idx] += ryij*tmpde;
        nbr->granny->fz[idx] += rzij*tmpde;
        nbr->granny->ep[idx] += tmpe;
        nbr->granny->unlock(idx);
      }
    }
    granny->lock(shifti);
    granny->fx[shifti] += pim*dex;
    granny->fy[shifti] += pim*dey;
    granny->fz[shifti] += pim*dez;
    granny->ep[shifti] += e;
    granny->unlock(shifti);
  }
	    
}



/// @brief Function to compute forces for a pair potential with SIMD optimized version
/// between particles of two specified types
/// @tparam Pot_t Type of the potential
/// @tparam symmetrize Indicates if the force computation is symmetrized
/// @tparam Verlet Indicates if the verlet lists are used
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
template < class Pot_t>
void leafCell::computeForcePair_perso_per_block(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  uint8_t *ptr_ti(getType());
  const auto&  param = pot->getParameters();
                
  const double   epsilon   ( param.epsilon ),
                 sigma6    ( param.sigma*param.sigma*param.sigma*param.sigma*param.sigma*param.sigma ),
                _2epsilon  (   2 * epsilon ),
                _24epsilon (  24 * epsilon ),
                 ecut      ( 0.5 * pot->getEcut() );
  
  double   rxij, ryij, rzij,
            rxi,  ryi,  rzi,
            dex,  dey,  dez,
              e, tmpe,tmpde,   
                      rcut2;
  
  uint8_t otherType;  
  NeighborList_base::nbr_id* start;  
  uint n, idx, end;
  
  rcut2 = Global::reference.getRcut2(typeIndexA, typeIndexB);
  
  // Loop on the particles of the cell
  for (size_t i=0; i<size ; ++i) {

    // Get the complementary type to the type of the particle
    if      (ptr_ti[i] == typeIndexA) otherType = typeIndexB;
    else if (ptr_ti[i] == typeIndexB) otherType = typeIndexA;
    else continue;
    
    double pim = Global::reference.getInvMass(ptr_ti[i]);
    double nim = -Global::reference.getInvMass(otherType);
    
    neighborListGetNeighbors(i, otherType, start, n);
    
    e=0.; tmpe=0.; tmpde=0.; dex=0.; dey=0.; dez=0.;

    int shifti = shift+i;

    const double rxi = granny->rx[shifti];
    const double ryi = granny->ry[shifti];
    const double rzi = granny->rz[shifti];
     
    leafCell ** neighborCellsList = neighborListGetLeafCellNeighbors(shifti);

    for(size_t j=0 ; j<n ; ++j)
    {  
      auto* nbr  = neighborCellsList[std::get<0>(start[j])];
      size_t idxs   = std::get<1>(start[j]) + nbr->shift;
      
      end = nbr->shift + nbr->getNumberOfParticles();

	assert(idxs<=end);

	// moche mais je laisse comme ça pour l'instant
      double *cell_ptr_rx = &(nbr->granny->rx[idxs]);   
      double *cell_ptr_ry = &(nbr->granny->ry[idxs]);  
      double *cell_ptr_rz = &(nbr->granny->rz[idxs]);  

      double *cell_ptr_fx = &(nbr->granny->fx[idxs]);   
      double *cell_ptr_fy = &(nbr->granny->fy[idxs]);  
      double *cell_ptr_fz = &(nbr->granny->fz[idxs]);  

      double *cell_ptr_ep = &(nbr->granny->ep[idxs]); 
      
      #pragma omp simd reduction(+:e,dex,dey,dez)
      #pragma vector aligned
      for(size_t s = 0; s < SIZEVECTOR; s++)
      {

        idx = idxs + s;

	if(idx <end )	{

        rxij = rxi - cell_ptr_rx[s];    
        ryij = ryi - cell_ptr_ry[s];    
        rzij = rzi - cell_ptr_rz[s];   
        
        const double r2 = rxij*rxij + ryij*ryij + rzij*rzij;
        
        if(r2<rcut2)
        {
          const double ir2     = 1.0/r2;
          const double ratio6  = sigma6*ir2*ir2*ir2; 
          const double ratio12 = ratio6*ratio6;
          const double diff = ratio12 -ratio6;

          tmpe  = _2epsilon  * (diff)   - ecut;
          tmpde = _24epsilon * (diff+ratio12) * ir2;
        
          e   += tmpe;
          dex += rxij * tmpde;
          dey += ryij * tmpde;
          dez += rzij * tmpde;
        
	  tmpde  *= nim;
	  cell_ptr_fx[s] += rxij*tmpde;
	  cell_ptr_fy[s] += ryij*tmpde;
	  cell_ptr_fz[s] += rzij*tmpde;
	  cell_ptr_ep[s] += tmpe;
        }
      }
	}
    }
    granny->fx[shifti] += pim*dex;
    granny->fy[shifti] += pim*dey;
    granny->fz[shifti] += pim*dez;
    granny->ep[shifti] += e;
  }
	    
}

/// @brief Function to compute forces for a pair potential with SIMD optimized version
/// between particles of two specified types
/// @tparam Pot_t Type of the potential
/// @tparam symmetrize Indicates if the force computation is symmetrized
/// @tparam Verlet Indicates if the verlet lists are used
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
template < class Pot_t>
void leafCell::computeForcePair_perso_per_block_mutex(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  uint8_t *ptr_ti(getType());
  const auto&  param = pot->getParameters();
                
  const double   epsilon   ( param.epsilon ),
                 sigma6    ( param.sigma*param.sigma*param.sigma*param.sigma*param.sigma*param.sigma ),
                _2epsilon  (   2 * epsilon ),
                _24epsilon (  24 * epsilon ),
                 ecut      ( 0.5 * pot->getEcut() );
  
  double   rxij, ryij, rzij,
            rxi,  ryi,  rzi,
            dex,  dey,  dez,
              e, tmpe,tmpde,   
                      rcut2;
  
  uint8_t otherType;  
  NeighborList_base::nbr_id* start;  
  uint n, idx, end;
  
  rcut2 = Global::reference.getRcut2(typeIndexA, typeIndexB);
  
  // Loop on the particles of the cell
  for (size_t i=0; i<size ; ++i) {

    // Get the complementary type to the type of the particle
    if      (ptr_ti[i] == typeIndexA) otherType = typeIndexB;
    else if (ptr_ti[i] == typeIndexB) otherType = typeIndexA;
    else continue;
    
    double pim = Global::reference.getInvMass(ptr_ti[i]);
    double nim = -Global::reference.getInvMass(otherType);
    
    neighborListGetNeighbors(i, otherType, start, n);
    
    e=0.; tmpe=0.; tmpde=0.; dex=0.; dey=0.; dez=0.;

    int shifti = shift+i;

    const double rxi = granny->rx[shifti];
    const double ryi = granny->ry[shifti];
    const double rzi = granny->rz[shifti];
    
    leafCell ** neighborCellsList = neighborListGetLeafCellNeighbors(shifti);

    for(size_t j=0 ; j<n ; ++j)
    {  
      auto* nbr  = neighborCellsList[std::get<0>(start[j])];
      size_t idxs   = std::get<1>(start[j]) + nbr->shift;
      
      end = nbr->shift + nbr->getNumberOfParticles();//idxs+SIZEVECTOR < nbr->getNumberOfParticles() ? idxs+SIZEVECTOR :  nbr->getNumberOfParticles() ;

	assert(idxs<=end);

	// moche mais je laisse comme ça pour l'instant
      double *cell_ptr_rx = &(nbr->granny->rx[idxs]);   
      double *cell_ptr_ry = &(nbr->granny->ry[idxs]);  
      double *cell_ptr_rz = &(nbr->granny->rz[idxs]);  

      double *cell_ptr_fx = &(nbr->granny->fx[idxs]);   
      double *cell_ptr_fy = &(nbr->granny->fy[idxs]);  
      double *cell_ptr_fz = &(nbr->granny->fz[idxs]);  

      double *cell_ptr_ep = &(nbr->granny->ep[idxs]);  

	bool ghost = nbr->granny->getIsGhost();
      
     for(size_t i = 0; i < SIZEVECTOR; i++)
      {
	 idx = idxs + i;
	if(idx<end && !ghost) nbr->granny->lock(idx);
      }



      #pragma omp simd reduction(+:e,dex,dey,dez)
      #pragma vector aligned
      for(size_t s = 0; s < SIZEVECTOR; s++)
      {

        idx = idxs + s;

	if(idx <end )	{

        rxij = rxi - cell_ptr_rx[s];    
        ryij = ryi - cell_ptr_ry[s];    
        rzij = rzi - cell_ptr_rz[s];   
        
        const double r2 = rxij*rxij + ryij*ryij + rzij*rzij;
        
        if(r2<rcut2)
        {
          const double ir2     = 1.0/r2;
          const double ratio6  = sigma6*ir2*ir2*ir2; 
          const double ratio12 = ratio6*ratio6;
          const double diff = ratio12 -ratio6;

          tmpe  = _2epsilon  * (diff)   - ecut;
          tmpde = _24epsilon * (diff+ratio12) * ir2;
        
          e   += tmpe;
          dex += rxij * tmpde;
          dey += ryij * tmpde;
          dez += rzij * tmpde;
        
	  tmpde  *= nim;
	  cell_ptr_fx[s] += rxij*tmpde;
	  cell_ptr_fy[s] += ryij*tmpde;
	  cell_ptr_fz[s] += rzij*tmpde;
	  cell_ptr_ep[s] += tmpe;
        }
      }
	}
     for(size_t i = 0; i < SIZEVECTOR; i++)
      {
	 idx = idxs + i;
	if(idx<end && !ghost) nbr->granny->unlock(idx);
      }
    }
    

    granny->lock(shifti);
    granny->fx[shifti] += pim*dex;
    granny->fy[shifti] += pim*dey;
    granny->fz[shifti] += pim*dez;
    granny->ep[shifti] += e;
    granny->unlock(shifti);
  }
	    
}

