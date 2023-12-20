#pragma once


/// @brief Fill a buffer with all the particles of the cell, case of a ParticleOutput buffer
/// @param [in,out] buffer Buffer to fill
/// @param [in] start Where to write this cell in the buffer
void Octree::fillBuffer(ParticleOutput* buffer, const uint start) {

  const auto ptr_id = getId();
  const auto ptr_ti = getType();
  const auto ptr_rx = getPositionX();
  const auto ptr_ry = getPositionY();
  const auto ptr_rz = getPositionZ();
  const auto ptr_vx = getVelocityX();
  const auto ptr_vy = getVelocityY();
  const auto ptr_vz = getVelocityZ();
          
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // crade mais c'est bien le MPI Comm World qui est utilisé

  // Set pointers to the buffer
  auto vi = buffer->id.data()   + start;
  auto vt = buffer->type.data() + start;
  auto vr = buffer->r.data()    + start;
  auto vv = buffer->v.data()    + start;

  /// @brief For each particle
  for (uint i=0; i<this->size; ++i) {

    // Write the index
    vi[i] = ptr_id[i]+1;

    // Write the type if needed
    if (buffer->writeType) 
      vt[i] = ptr_ti[i];

    // Write the position
    vr[i].x = ptr_rx[i];    
    vr[i].y = ptr_ry[i];    
    vr[i].z = ptr_rz[i];    

    // Write the velocity if needed
    if (buffer->writeVelocity) {
      vv[i].x = neighborList.numberOfNeighbors(i);   //WARNING 
      vv[i].y = level;    
      vv[i].z = rank;    
    }

  }


}

/// @brief Fill a buffer with all the particles of the cell, case of a ParticleInSitu buffer
/// @param [in,out] buffer Buffer to fill
/// @param [in] start Where to write this cell in the buffer
void Octree::fillBuffer(ParticleInSitu* buffer, const uint start) {  


  auto ptr_id = getId();
  auto ptr_ti = getType();
  auto ptr_rx = getPositionX();
  auto ptr_ry = getPositionY();
  auto ptr_rz = getPositionZ();
  auto ptr_vx = getVelocityX();
  auto ptr_vy = getVelocityY();
  auto ptr_vz = getVelocityZ();

  // Set pointers to the buffer
  auto vi = buffer->id   + start;
  auto vt = buffer->type + start;
  auto vrx = buffer->rx    + start;
  auto vry = buffer->ry    + start;
  auto vrz = buffer->rz    + start;
  auto vvx = buffer->vx    + start;
  auto vvy = buffer->vy    + start;
  auto vvz = buffer->vz    + start;

  /// @brief For each particle
  for (uint i=0; i<this->size; ++i) {

    // Write the index
    vi[i] = ptr_id[i];

    // Write the type
    vt[i] = ptr_ti[i];

    // Write the position
    vrx[i] = ptr_rx[i];    
    vry[i] = ptr_ry[i];    
    vrz[i] = ptr_rz[i];    

    // Write the velocity
    vvx[i] = ptr_vx[i];    
    vvy[i] = ptr_vy[i];    
    vvz[i] = ptr_vz[i];    

  }


}

/// @brief Fill a buffer with all the particles of the ghost cell, case of a ParticleInSitu buffer
/// @param [in,out] buffer Buffer to fill
/// @param [in] start Where to write this cell in the buffer
void Octree::fillGhostBuffer(ParticleInSitu* buffer, const uint start) {  


  const auto& ptr_id = getId();
  const auto& ptr_rx = getPositionX();
  const auto& ptr_ry = getPositionY();
  const auto& ptr_rz = getPositionZ();

  // Set pointers to the buffer
  auto vi = buffer->ghost_id   + start;
  auto vrx = buffer->ghost_rx    + start;
  auto vry = buffer->ghost_ry    + start;
  auto vrz = buffer->ghost_rz    + start;

  /// @brief For each particle
  for (uint i=0; i<this->size; ++i) {

    // Write the index
    vi[i] = ptr_id[i];

    // Write the position
    vrx[i] = ptr_rx[i];    
    vry[i] = ptr_ry[i];    
    vrz[i] = ptr_rz[i];
    
  }


}


/// @brief Fill a buffer with all the particles of the cell, case of a Hercule dump buffer
/// @tparam DumpStruct Dump structure (hercule or hercule for DPDE)
/// @param [in,out] _buffer Buffer to fill
/// @param [in] start Where to write this cell in the buffer
// Needed to fill the Hercule dump
template <class DumpStruct> void Octree::fillBuffer(DumpStruct& _buffer, const uint start) {

// not implemented
  
}

/// @brief Fill a buffer with the cell data
/// @param [in,out] buffer Buffer to fill
/// @param [in] start Where to write this cell in the buffer
void Octree::fillBuffer(CellOutput* buffer, const uint start) {

  // Set pointers to the buffer
  auto vp = buffer->p.data();
  auto vt = buffer->t.data();
  auto vd = buffer->d.data();

  // Convert and write the pressure tensor
  vp[start] = this->pressureTensor/(3.*Global::domainInfo.getCellVolume());
  vp[start].m11 = convert(vp[start].m11, Stamp_Units::pressure, SI_Units_base::pascal);
  vp[start].m12 = convert(vp[start].m12, Stamp_Units::pressure, SI_Units_base::pascal);
  vp[start].m13 = convert(vp[start].m13, Stamp_Units::pressure, SI_Units_base::pascal);
  vp[start].m21 = convert(vp[start].m21, Stamp_Units::pressure, SI_Units_base::pascal);
  vp[start].m22 = convert(vp[start].m22, Stamp_Units::pressure, SI_Units_base::pascal);
  vp[start].m23 = convert(vp[start].m23, Stamp_Units::pressure, SI_Units_base::pascal);
  vp[start].m31 = convert(vp[start].m31, Stamp_Units::pressure, SI_Units_base::pascal);
  vp[start].m32 = convert(vp[start].m32, Stamp_Units::pressure, SI_Units_base::pascal);
  vp[start].m33 = convert(vp[start].m33, Stamp_Units::pressure, SI_Units_base::pascal);
  
/*  // Convert and write the temperature
  auto cellMomentum = this->computeTotalMomentum();
  double cellMass = this->computeTotalMass();
  double shiftedKE = this->computeShiftedKineticEnergy(cellMomentum/cellMass);
  vt[start] = convert(shiftedKE/auxMax(uint(1),numberOfParticles), Stamp_Units::energy, SI_Units_base::joule) / (1.5*Constant::boltzmann);
  
  // Convert and write the density
  vd[start] = convert(cellMass/Global::domainInfo.getCellVolume(), Stamp_Units::mass/(Stamp_Units::length*Stamp_Units::length*Stamp_Units::length),SI_Units_base::kilogram/(SI_Units_base::meter*SI_Units_base::meter*SI_Units_base::meter));*/
  
}

