/// @brief Function to compute embedding terms for a EAM potential. Don't forget to modify Verlet version.
/// between particles of two specified types
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
void Octree::computeForceEmb(EAMPotential* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {
  
  // Loop on the particles of the cell
  double *e = this->m_eamStorage.emb();
  double *r = this->m_eamStorage.rho();
  double fEmbed;

  assert(this->m_eamStorage.size() == this->size);
  
 // #pragma omp simd
 // //#pragma vector aligned
  for (uint i=0; i<this->size; i++) {

    // Compute embedding term only if the particle is of one of the considered types
    if (ti[i]!=typeIndexA || ti[i]!=typeIndexB) continue;

    fEmbed = 0.;

    // Compute embedding term from rho
    pot->fEmbed(r[i], fEmbed, e[i]);
    ep[i] += fEmbed;

  }

}

/// @brief Compute the potential energy on the system
/// @return Total potential energy in the cell
inline double Octree::computePotentialEnergyAMR() const {

  double ePot = 0.;

  // Sum on all the particles
  for (uint i=0; i< this->size; ++i) {
    
    // avoid NaN
    assert(Global::reference.getInvMass(ti[i]) != 0) ; 
    ePot +=  ep[i];
  }
  
  return ePot;

}
