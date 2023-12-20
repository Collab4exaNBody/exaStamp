
/// @brief Bring back computed forces from the vectorization buffer to the particle (Verlet lists)
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] cells Pointer to the cells of the grid
/// @param [in] i Index of the particle
/// @param [in] typeIndex Type of the neighbors to consider
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] symmetrize Indicates if the force computation is symmetrized
/// @param [in] nbrSize number of element in the buffer
void leafCell::writeForce(const uint i, const uint8_t typeIndex, const bool symmetrize, uint nbrSize) {
  // Inverse mass for particle and neighbors
  double pim = Global::reference.getInvMass(getType(i));
  double nim = Global::reference.getInvMass(typeIndex);

  // Get neighbors data
  NeighborList_base::nbr_id* start;

  uint nbrSize2;
  neighborListGetNeighbors(i, typeIndex, start, nbrSize2);
  
  assert(nbrSize2 >= nbrSize);

  // Buffer to read
  auto& fb = self_type::getVectBuffer();
  
  // Sum the contributions for the particle
  double local_fe[4] = {0.};
  mat3<double> virial(0.);
  for (uint j=0; j<nbrSize; ++j) {
    local_fe[0] += fb.dfx(j);
    local_fe[1] += fb.dfy(j);
    local_fe[2] += fb.dfz(j);
    local_fe[3] += fb.den(j);
    // Virial for particle i
    virial += tensor(vec3<double>(fb.dfx(j),fb.dfy(j),fb.dfz(j)),vec3<double>(fb.drx(j),fb.dry(j),fb.drz(j)));
  }


  // Update force and energy on the particle
  /**/ incForceX(i,pim*local_fe[0]);
  /**/ incForceY(i,pim*local_fe[1]);
  /**/ incForceZ(i,pim*local_fe[2]);
  /**/ incPontentialEnergy(i,local_fe[3]);
  // Update pressure tensor
  this->pressureTensor += virial;

  // If symmetrized force computation
  if (symmetrize) {

    leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(i+this->shift);
  
    // Loop on neighbor
    uint indx;
    mat3<double> tmp_tensor(0.);
    for (uint k=0; k<nbrSize; ++k) {

      indx = fb.indxBuffer(k);
      const auto& nbrCellIndex = std::get<NeighborList_base::CELL>(start[indx]);

      // Get neighbor cell and index
      auto& nbrCell  = neighborCellsList[nbrCellIndex];
      auto& nbrIndex = std::get<NeighborList_base::INDX>(start[indx]);
      	
      // If the neighbor is not a ghost
      if(!nbrCell->getIsGhost()) {

        // Update force and energy on the neighbor
      	
      	double* ptr_fx = nbrCell->getForceX();
       	double* ptr_fy = nbrCell->getForceY();
      	double* ptr_fz = nbrCell->getForceZ();
      	double* ptr_ep = nbrCell->getPotentialEnergy();

      	ptr_fx[nbrIndex] -= nim*fb.dfx(k);
      	ptr_fy[nbrIndex] -= nim*fb.dfy(k);
      	ptr_fz[nbrIndex] -= nim*fb.dfz(k);
      	ptr_ep[nbrIndex] += fb.den(k);

	      /**/nbrCell->pressureTensor +=  tensor(vec3<double>(fb.dfx(k),fb.dfy(k),fb.dfz(k)),vec3<double>(fb.drx(k),fb.dry(k),fb.drz(k)));

      }
    }
  }
}


/// @brief Bring back computed rho terms from the vectorization buffer to the EAM storage (Verlet lists)
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] cells Pointer to the cells of the grid
/// @param [in] i Index of the particle
/// @param [in] typeIndex Type of the neighbors to consider
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] symmetrize Indicates if the force computation is symmetrized
/// @param [in] nbrSize number of element in the buffer
void leafCell::writeRho(const uint i, const uint8_t typeIndex, const bool symmetrize, uint nbrSize) {
  // Get neighbor data
  NeighborList_base::nbr_id* start;

  uint nbrSize2;
  neighborListGetNeighbors(i, typeIndex, start, nbrSize2);
  
  assert(nbrSize2 >= nbrSize);

  // Buffer to fill
  auto& fb = self_type::getVectBuffer();
  // Sum the contributions for the particle
  double local_rho = 0.;

  for (uint j=0; j<nbrSize; ++j) local_rho += fb.rho(j);
  
  // Update rho term on the particle
  this->rhoAMR(i) += local_rho; 


  leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(i+this->shift);

  // If symmetrized force computation
  if (symmetrize) {
    // Loop on neighbors
    for (uint k=0; k<nbrSize; ++k) {
      
    	// Get neighbor cell and index
    	self_type* nbrCell = neighborCellsList[std::get<NeighborList_base::CELL>(start[fb.indxBuffer(k)])];
    	uint nbrIndex      = std::get<NeighborList_base::INDX>(start[fb.indxBuffer(k)]);
    	
      assert(nbrIndex>=0);
      
    	// Update rho term on the neighbor
    	if(!nbrCell->getIsGhost())
    	  nbrCell->rhoAMR(nbrIndex) += fb.rho(k);

    }
  }
}

/// @file
/// @brief MEAM force writings for CellList


/// @brief Bring back computed forces of pair interaction from the vectorization buffer to the particle
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] cells Pointer to the cells of the grid
/// @param [in] cell number of the Cell
/// @param [in] i Index of the particle
/// @param [in] typeIndex Type of the neighbors to consider
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] size Number of neighbors
template<bool mutex> void leafCell::writeForceMEAM( const uint cell, const uint i, const uint8_t typeIndex, const uint8_t* ghostLayer, const uint sizeOfInteraction) {

  // Inverse mass for particle and neighbors
  double pim = Global::reference.getInvMass(getType(i));
  double nim = Global::reference.getInvMass(typeIndex);

  // Get neighbors data
  NeighborList_base::nbr_id* start;
  uint size1;
  neighborListGetNeighbors(i, typeIndex, start, size1);

  // Buffer to read
  auto& fb = self_type::getVectBufferWithScreenterm();

  // Sum the contributions for the particle
  double local_fe[4] = {0.};
  mat3<double> virial(0.);
  for (uint j=0; j<sizeOfInteraction; ++j) {
  local_fe[0] += fb.dfx(j);
  local_fe[1] += fb.dfy(j);
  local_fe[2] += fb.dfz(j);
  local_fe[3] += fb.den(j);
  // Virial for particle i
  virial += tensor(vec3<double>(fb.dfx(j),fb.dfy(j),fb.dfz(j)),vec3<double>(fb.drx(j),fb.dry(j),fb.drz(j)));
  }
  
  if(mutex)
  {
    if (!granny->getIsGhost())
    {
      uint positionAtom = i + shift;
      granny->lock(positionAtom);    
      
      granny->ep[positionAtom] += local_fe[3];
      granny->fx[positionAtom] -= pim*local_fe[0];
      granny->fy[positionAtom] -= pim*local_fe[1];
      granny->fz[positionAtom] -= pim*local_fe[2];
      this->pressureTensor += virial;

      granny->unlock(positionAtom);
    }
  }
  else
  {
    if (ghostLayer[cell]==0)
    {
      // Update derivate of rho 0,1,2,3 on the particle    
      /**/ incForceX(i,-pim*local_fe[0]);
      /**/ incForceY(i,-pim*local_fe[1]);
      /**/ incForceZ(i,-pim*local_fe[2]);
      /**/ incPontentialEnergy(i,local_fe[3]);
       this->pressureTensor += virial;
    }
  }


  // Loop on neighbor
  for (uint k=0; k<sizeOfInteraction; ++k)
  {

	leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(i+this->shift);
	uint indx = fb.indxBuffer(k);
	const auto& nbrCellIndex = std::get<NeighborList_base::CELL>(start[indx]);

	// Get neighbor cell and index
	self_type* nbrCell  = neighborCellsList[nbrCellIndex];
	uint nbrIndex = std::get<NeighborList_base::INDX>(start[indx]);

        // Update force and energy on the neighbor
      	
      	double* ptr_fx = nbrCell->getForceX();
       	double* ptr_fy = nbrCell->getForceY();
      	double* ptr_fz = nbrCell->getForceZ();


	if(mutex)
	{
		if (nbrCell->getIsGhost()) continue;

		nbrCell->granny->lock(nbrIndex+nbrCell->shift);

		ptr_fx[nbrIndex] += nim*fb.dfSx(k);
		ptr_fy[nbrIndex] += nim*fb.dfSy(k);
		ptr_fz[nbrIndex] += nim*fb.dfSz(k);

		nbrCell->granny->unlock(nbrIndex+nbrCell->shift);
	}
	else
	{
		// If the neighbor is not a ghost
		if (!nbrCell->granny->getIsGhost())
		{
			ptr_fx[nbrIndex] += nim*fb.dfSx(k);
			ptr_fy[nbrIndex] += nim*fb.dfSy(k);
			ptr_fz[nbrIndex] += nim*fb.dfSz(k);
		}
	}
  }
}


/// @brief Bring back computed forces (rho part) from the vectorization buffer to the MEAM storage
/// @tparam Pot_t Class of potential
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] pot Potential
/// @param [in] cells Pointer to the cells of the grid
/// @param [in] cell number of the Cell
/// @param [in] i Index of the particle
/// @param [in] typeIndex Type of the neighbors to consider
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] size Number of neighbors
/// @param [in] rho0 rho 0
/// @param [in] rhod0 Derivative of rho 0
/// @param [in] rho1 rho 1
/// @param [in] rhod1 Derivative of rho 1
/// @param [in] rho2 rho 2
/// @param [in] rhod2 Derivative of rho 2
/// @param [in] rho3 rho 3
/// @param [in] rhod3 Derivative of rho 3
template <bool mutex, class Pot_t>
void leafCell::writeRho_0123_MEAM(Pot_t* pot, const uint cell, const uint i, const uint8_t typeIndex, const uint8_t* ghostLayer, const uint sizeOfInteraction, double rho0, vec3<double> &rhod0, double rho1, vec3<double> &rhod1, double rho2, vec3<double> &rhod2, double rho3, vec3<double> &rhod3)
{

  // Inverse mass for particle
  double pim = Global::reference.getInvMass(getType(i));
  double tmpRho0,tmpRho1,tmpRho2,tmpRho3,tmpG,tmpdG;
  double rho,fEmbed,df;

  // Get neighbors data
  NeighborList_base::nbr_id* start;
  uint size1;
  neighborListGetNeighbors(i, typeIndex, start, size1);

  // Buffer to read
  auto& fb = self_type::getVectBufferWithScreenterm();

  // Apply the potential to get the rho terms
  pot->rho(rho, rho0,rho1,rho2,rho3);

  // compute F(rho/Z)
  pot->fEmbed(rho, fEmbed, df);

  pot->getCoeffDerivateRho(rho0,rho1,rho2,rho3, tmpRho0, tmpRho1, tmpRho2, tmpRho3, tmpG, tmpdG );

  double Coeff=pim*df/ pot->p.Z; // Z=12

  if(mutex)
  {
    if (!granny->getIsGhost())
    {
      uint positionAtom = i + shift;
      granny->lock(positionAtom);    
      
      granny->ep[positionAtom] += fEmbed;
      granny->fx[positionAtom] -= Coeff*(tmpG*rhod0.x+rho0*tmpdG*(tmpRho0*rhod0.x+tmpRho1*rhod1.x+tmpRho2*rhod2.x+tmpRho3*rhod3.x));
      granny->fy[positionAtom] -= Coeff*(tmpG*rhod0.y+rho0*tmpdG*(tmpRho0*rhod0.y+tmpRho1*rhod1.y+tmpRho2*rhod2.y+tmpRho3*rhod3.y));
      granny->fz[positionAtom] -= Coeff*(tmpG*rhod0.z+rho0*tmpdG*(tmpRho0*rhod0.z+tmpRho1*rhod1.z+tmpRho2*rhod2.z+tmpRho3*rhod3.z));

      granny->unlock(positionAtom);
    }
  }
  else
  {
    if (!granny->getIsGhost())
    {

      uint positionAtom = i + shift;

      // Update derivate of rho 0,1,2,3 on the particle    
      granny->ep[positionAtom] += fEmbed;
      granny->fx[positionAtom] -= Coeff*(tmpG*rhod0.x+rho0*tmpdG*(tmpRho0*rhod0.x+tmpRho1*rhod1.x+tmpRho2*rhod2.x+tmpRho3*rhod3.x));
      granny->fy[positionAtom] -= Coeff*(tmpG*rhod0.y+rho0*tmpdG*(tmpRho0*rhod0.y+tmpRho1*rhod1.y+tmpRho2*rhod2.y+tmpRho3*rhod3.y));
      granny->fz[positionAtom] -= Coeff*(tmpG*rhod0.z+rho0*tmpdG*(tmpRho0*rhod0.z+tmpRho1*rhod1.z+tmpRho2*rhod2.z+tmpRho3*rhod3.z));
    }
  }

  // Loop on neighbor
  for (uint k=0; k<sizeOfInteraction; ++k)
  {

	leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(i+this->shift);
	uint indx = fb.indxBuffer(k);
	const auto& nbrCellIndex = std::get<NeighborList_base::CELL>(start[indx]);

	// Get neighbor cell and index
	self_type* nbrCell  = neighborCellsList[nbrCellIndex];
	uint nbrIndex = std::get<NeighborList_base::INDX>(start[indx]);
      	
      	double* ptr_fx = nbrCell->getForceX();
       	double* ptr_fy = nbrCell->getForceY();
      	double* ptr_fz = nbrCell->getForceZ();

	if(mutex)
	{
		if (nbrCell->granny->getIsGhost()) continue;

		nbrCell->granny->lock(nbrIndex+nbrCell->shift);

		ptr_fx[nbrIndex] += Coeff*(rho0*tmpdG*(tmpRho0*fb.drho0x(k)+tmpRho1*fb.drho1x(k)+tmpRho2*fb.drho2x(k)+tmpRho3*fb.drho3x(k))+tmpG*fb.drho0x(k));
		ptr_fy[nbrIndex] += Coeff*(rho0*tmpdG*(tmpRho0*fb.drho0y(k)+tmpRho1*fb.drho1y(k)+tmpRho2*fb.drho2y(k)+tmpRho3*fb.drho3y(k))+tmpG*fb.drho0y(k));
		ptr_fz[nbrIndex] += Coeff*(rho0*tmpdG*(tmpRho0*fb.drho0z(k)+tmpRho1*fb.drho1z(k)+tmpRho2*fb.drho2z(k)+tmpRho3*fb.drho3z(k))+tmpG*fb.drho0z(k));

		nbrCell->granny->unlock(nbrIndex+nbrCell->shift);
	}
	else
	{

		// If the neighbor is real
		if (!nbrCell->granny->getIsGhost())
		{
			ptr_fx[nbrIndex] += Coeff*(rho0*tmpdG*(tmpRho0*fb.drho0x(k)+tmpRho1*fb.drho1x(k)+tmpRho2*fb.drho2x(k)+tmpRho3*fb.drho3x(k))+tmpG*fb.drho0x(k));
			ptr_fy[nbrIndex] += Coeff*(rho0*tmpdG*(tmpRho0*fb.drho0y(k)+tmpRho1*fb.drho1y(k)+tmpRho2*fb.drho2y(k)+tmpRho3*fb.drho3y(k))+tmpG*fb.drho0y(k));
			ptr_fz[nbrIndex] += Coeff*(rho0*tmpdG*(tmpRho0*fb.drho0z(k)+tmpRho1*fb.drho1z(k)+tmpRho2*fb.drho2z(k)+tmpRho3*fb.drho3z(k))+tmpG*fb.drho0z(k));
		}
	}
  }
}
