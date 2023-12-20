#pragma once


/// @brief Fill the vectorization buffer of a particle with its neighbors of a specified type
/// @tparam fill_emb Indicates if the buffer must be filled with embedding terms too
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] cells All the cells in the grid
/// @param [in] i Index of the particle
/// @param [in] typeIndex Type of the neighbors to consider
/// @param [in]  correcter correct distance on edge cells.
/// @return Number of neighbors
template <bool fill_emb> uint leafCell::fillForceBufferVerlet(const uint i, const uint8_t typeIndex) {

  // Get neighbors data
  NeighborList_base::nbr_id* start;
  uint n;
  neighborListGetNeighbors(i, typeIndex, start, n);

  // Buffer to fill
  auto& fb = self_type::getVectBuffer();
  
  const auto& ptr_rx = getPositionX();
  const auto& ptr_ry = getPositionY();
  const auto& ptr_rz = getPositionZ();
  
  const double xi = ptr_rx[i];
  const double yi = ptr_ry[i];
  const double zi = ptr_rz[i];  

  leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(this->shift+i);
  
  // Loop on neighbors
  for (uint j=0; j<n; ++j) 
  {
    // Get neighbor cell and index
    const auto& nbr   = neighborCellsList[std::get<NeighborList_base::CELL>(start[j])];
    const auto& idx   = std::get<NeighborList_base::INDX>(start[j]) + nbr->shift;

    assert ( idx + nbr->shift < nbr->granny->size );
 
    if (fill_emb)
      assert( nbr->granny->size == nbr->granny->m_eamStorage.size() );
      
    // Fill the buffer with distances
    fb.drx(j) = xi - nbr->granny->rx[idx];
    fb.dry(j) = yi - nbr->granny->ry[idx];
    fb.drz(j) = zi - nbr->granny->rz[idx];
    fb.indxBuffer(j)=j;  

    // Fill the buffer with embedding terms if required
    if (fill_emb)
     fb.emb(j) = nbr->embAMR(idx-nbr->shift);
  }

  // Atom interactions are computed by typeAtoms.
  double rcut2 = Global::reference.getRcut2(getType(i), typeIndex);
  
  // Find atoms too far away
  simd::kernels::Verlet_corrector(fb.drx(), fb.dry(), fb.drz(), rcut2, fb.sij(), n); // fb.sij = tmp3

  for(int j=n-1; j>=0; --j) 
  {
    // Overwrite atoms too far away
    if( fb.sij(j) != 1 )
    {
      n--;
      fb.drx(j) = fb.drx(n);
      fb.dry(j) = fb.dry(n);
      fb.drz(j) = fb.drz(n);
      fb.indxBuffer(j)=fb.indxBuffer(n);     

      if (fill_emb)
        fb.emb(j) = fb.emb(n);
    }
  }

  // Return the real number of neighbors
  return n;

}


/// @brief Fill the vectorization buffer of a particle with its neighbors of a specified type
/// @tparam fill_emb Indicates if the buffer must be filled with embedding terms too
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] cells All the cells in the grid
/// @param [in] i Index of the particle
/// @param [in] typeIndex Type of the neighbors to consider
/// @param [in]  correcter correct distance on edge cells.
/// @return Number of neighbors
uint leafCell::fillForceBufferVerletMEAM(const uint i, const uint8_t typeIndex) {

  // Get neighbors data
  NeighborList_base::nbr_id* start;
  uint n;
  neighborListGetNeighbors(i, typeIndex, start, n);

  // Buffer to fill
  auto& fb = self_type::getVectBufferWithScreenterm();
  fb.resizeMeam(n);
  
  const auto& ptr_rx = getPositionX();
  const auto& ptr_ry = getPositionY();
  const auto& ptr_rz = getPositionZ();
  
  const double xi = ptr_rx[i];
  const double yi = ptr_ry[i];
  const double zi = ptr_rz[i];  

  
  leafCell** neighborCellsList = neighborListGetLeafCellNeighbors(this->shift+i);

  // Loop on neighbors
  for (size_t j=0; j<n; ++j) 
  {
    // Get neighbor cell and index
    const auto& nbr   = neighborCellsList[std::get<NeighborList_base::CELL>(start[j])];
    const auto& idx   = std::get<NeighborList_base::INDX>(start[j])+nbr->shift;

    assert ( idx + nbr->shift < nbr->granny->size );
      
    // Fill the buffer with distances
    fb.drx(j) = xi - nbr->granny->rx[idx];
    fb.dry(j) = yi - nbr->granny->ry[idx];
    fb.drz(j) = zi - nbr->granny->rz[idx];
    fb.indxBuffer(j)=j;  
  }

  // Atom interactions are computed by typeAtoms.
  double rcut2 = Global::reference.getRcut2(getType(i), typeIndex);
  
  // Find atoms too far away
  simd::kernels::Verlet_corrector(fb.drx(), fb.dry(), fb.drz(), rcut2, fb.sij(), n); // fb.sij = tmp3


  for(int j=n-1; j>=0; --j) 
  {
    // Overwrite atoms too far away
    if( fb.sij(j) != 1 )
    {
      n--;
      fb.drx(j) = fb.drx(n);
      fb.dry(j) = fb.dry(n);
      fb.drz(j) = fb.drz(n);
      fb.indxBuffer(j)=fb.indxBuffer(n);     
    }
  }

  // Return the real number of neighbors
  return n;

}

/// @brief Resize buffer if interaction was screened
/// @param [in] i Index of the particle
/// @param [in] typeIndex Type of the neighbors to consider
/// @return Number of neighbors
inline uint leafCell::resizeForceBuffer(uint sizeBuffer) {

  // Buffer to fill
  auto& fb = self_type::getVectBufferWithScreenterm();

  for(ssize_t j=sizeBuffer-1; j>=0; --j) 
  {
    if(fb.S(j) == 0)
    {
      sizeBuffer--;
      fb.drx(j)        = fb.drx(sizeBuffer);
      fb.dry(j)        = fb.dry(sizeBuffer);
      fb.drz(j)        = fb.drz(sizeBuffer);
      fb.S(j)          = fb.S(sizeBuffer);
      fb.indxBuffer(j) = fb.indxBuffer(sizeBuffer);
    }
  }
  return sizeBuffer;
} 
