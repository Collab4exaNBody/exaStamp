inline void Octree::ghostVelCopy(const int destDomain, const vec3<int>& destCell, std::vector< std::tuple<int, ExchangeGhostV3> >& leaving) const{

  auto tmp = std::make_tuple(destDomain, ExchangeGhostV3());

  std::get<1>(tmp).cell = destCell;

  auto& id_ = std::get<1>(tmp).id;
  auto& vel = std::get<1>(tmp).v3;

  for (uint i=0; i< this->size; ++i) {
    id_ = id[i];
    vel = vec3<double>(vx[i],vy[i],vz[i]);
    leaving.push_back(tmp);
  }
}


/// @brief Copy the particle velocities into a local ghost cell
/// @param [out] to Ghost cell to fill
inline void Octree::ghostVelCopy(self_type& to) const {
  to.vx = vx;
  to.vy = vy;
  to.vz = vz;
}


/// @brief Set the embedding term of a particle from exchanged data
/// @param [in] exchange Received data
inline void Octree::setEmbAMR(const ExchangeEAM& exchange) {

  // Loop on the particles of the cell
  for (uint i=0; i< this->size; ++i) {

    // If the received embedding term is of this particle
    if (getId(i)==exchange.id) {

      // Set the embedding term
      this->embAMR(i) = exchange.emb;

    }

  }

}

/// @brief Get the destination cell of a particle
/// @param [in] i Index of the particle
/// @param [in] minBounds Lower boundaries of the cell
/// @param [in] maxBounds upper boundaries of the cell
/// @return Destination cell, Neighbor::null if the cell doesn't change cell
inline int8_t Octree::computeNeighborMov(const uint i, const vec3<double>& minBounds, const vec3<double>& maxBounds) const {

  vec3<int> m(0);

  // Get the movement of the particle in the cells
  if      (getPositionX(i) >= maxBounds.x) m.x += 1;
  else if (getPositionX(i) <  minBounds.x) m.x -= 1;
  
  if      (getPositionY(i) >= maxBounds.y) m.y += 1;
  else if (getPositionY(i) <  minBounds.y) m.y -= 1;
  
  if      (getPositionZ(i) >= maxBounds.z) m.z += 1;
  else if (getPositionZ(i) <  minBounds.z) m.z -= 1;
  
  // Get the index of the destination cell
  return Neighbor::getIndex(m);

}


/// @brief Store the particles leaving the cell and remove them
/// @param [out] list Array to store the leaving particles
/// @param [in] nbrs Cell neighbors
/// @param [in] minBounds Cell lower boundaries
/// @param [in] maxBounds Cell upper boundaries
void Octree::internalReorganization(std::list< std::tuple<int, MPI__Particle> >& list, const Array<int>& nbrs, const vec3<double>& minBounds, const vec3<double>& maxBounds) {

	// For each particle
  for (uint i=0; i<this->size; ++i) {

  	// Get the movement of the particle
    auto mov = computeNeighborMov(i, minBounds, maxBounds);

    // If leaving
    if (mov != Neighbor::null) {

    	// Create associated MPI particle
      MPI__Particle tmp;
      getExchange(tmp, i);

      // Add to leaving array and mark to remove
      list.emplace_back(nbrs[mov], std::move(tmp));
      this->mark(i);

    }

  }

  // Remove each particle whose index is in the marked list
  for (uint i=this->marks.size()-1; i!=(uint)-1; --i) {
    remove(this->marks[i]);
  }
  
  this->clear_marked();

}

/// @brief Convert a particle of the cell into an MPI type
/// @param [in] i Index of the particle
/// @param [out] p Particle in MPI type
inline void Octree::getExchange(MPI__Particle& p, const uint i) const {
  p.id  = getId(i);
  p.ti  = getType(i);
  p.r.x = getPositionX(i);
  p.r.y = getPositionY(i);
  p.r.z = getPositionZ(i);
  p.v.x = getVelocityX(i);
  p.v.y = getVelocityY(i);
  p.v.z = getVelocityZ(i);

}




