/// @brief Shortcut for the arguments of a MEAM force computation (rho0)
#define _rho0_  fb.drho0x(),fb.drho0y(),fb.drho0z()
/// @brief Shortcut for the arguments of a MEAM force computation (rho1)
#define _rho1_   rho1,fb.drho1x(),fb.drho1y(),fb.drho1z()
/// @brief Shortcut for the arguments of a MEAM force computation (rho2)
#define _rho2_   rho2,fb.drho2x(),fb.drho2y(),fb.drho2z()
/// @brief Shortcut for the arguments of a MEAM force computation (rho3)
#define _rho3_   rho3,fb.drho3x(),fb.drho3y(),fb.drho3z()

/// @brief Function to compute forces for a MEAM potential with SIMD optimized version. Don't forget to modify linked cell version.
/// between particles of two specified types
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
/// @param [in] cell the numeration of cell
/// @param [in] cells All the cells in the grid
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] correcter Tool to correct distances
template < bool mutex, class Pot_t>
void leafCell::computeForcesAndPotentialMEAMVerlet_AMR(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, const int cell, const uint8_t * ghostLayer)
{
  double rho0,rho1,rho2,rho3;

  // Get vectorization buffer
  auto& fb = self_type::getVectBufferWithScreenterm();

  // Loop on the particles of the cell
  for (uint i=0; i<size; ++i) 
  {
    // Get the complementary type to the type of the particle
    uint8_t otherType;
    if    (getType(i) == typeIndexA) otherType = typeIndexB;
    else if (getType(i) == typeIndexB) otherType = typeIndexA;
    else continue;

    vec3<double> rhod0,rhod1,rhod2,rhod3;

    // Fill the buffer with neighbors of that type
    uint nbrSizeLast = fillForceBufferVerletMEAM(i, otherType);

    if(nbrSizeLast>0) {

    	//compute Screening term
    	pot->screeningFunction(fb.drx(), fb.dry(), fb.drz(), fb.S(), nbrSizeLast);

    	//resize buffer if pair interactions are screened
    	uint nbrSize = resizeForceBuffer(nbrSizeLast);

    	// Apply the potential to get the forces

    	// -- Pair term --
    	//compute phi_i
    	pot->phi(fb.dfx(), fb.dfy(), fb.dfz(), fb.dfSx(), fb.dfSy(), fb.dfSz(), fb.den(), fb.drx(), fb.dry(), fb.drz(), fb.S(), nbrSize);

    	//write forces and energy of pair term
    	writeForceMEAM<mutex>( cell, i, otherType, ghostLayer, nbrSize);

    	// -- Electronique term --
    	// compute rho0
    	rho0=pot->compute_rho_0(fb.S(), fb.drx(), fb.dry(), fb.drz(), nbrSize);

    	//to avoid issue because if rho0=0 => rho=0 and log(rho) * rho != 0
    	rho0=std::max(rho0, 1.e-30);

    	//compute first part of rho0,rho 2, rho 3 and rho1 (pour l'instant)
    	pot->derivedRhoS(fb.S(),fb.drx(), fb.dry(), fb.drz(), _rho0_, rhod0, _rho1_, rhod1, _rho2_, rhod2, _rho3_, rhod3, nbrSize);

    	//write energy of F(rho/Z) end first part of derivate of F(rho/Z)
    	writeRho_0123_MEAM<mutex>(pot, cell, i, otherType, ghostLayer, nbrSize, rho0, rhod0, rho1, rhod1, rho2, rhod2, rho3, rhod3);
    }
  }
}

#undef _rho0_
#undef _rho1_
#undef _rho2_
#undef _rho3_

