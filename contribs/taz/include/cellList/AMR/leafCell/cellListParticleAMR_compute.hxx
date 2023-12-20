// ---------------   Verlet Version  ------------------- //

/// @brief Function to compute forces for a pair potential with SIMD optimized version
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
template < class CList>
void leafCell::computeForcePairVerlet(PairPotential* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter, const bool symmetrize) {

	    
}

/// @brief Function to compute forces for a EAM potential with SIMD optimized version
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
template < class CList>
void leafCell::computeForceRhoVerlet(EAMPotential* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter, const bool symmetrize) {


	    
}

/// @brief Function to compute forces for a EAM potential with SIMD optimized version
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
template < class CList>
void leafCell::computeForceFinalVerlet(EAMPotential* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter, const bool symmetrize) {


	    
}


/// @brief Function to compute forces for a pair potential with SIMD optimized version
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
template < class Pot_t, class CList, class>
void leafCell::computeForcePairVerlet(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter, const bool symmetrize) {

  if(isleaf == 1)
  #ifndef _perBLOCK
    leafCell::computeForcePair_perso(pot, typeIndexA, typeIndexB) ;
  #else  
    leafCell::computeForcePair_perso_per_block<Pot_t>(pot, typeIndexA, typeIndexB);
  #endif
      
  else 
   if(symmetrize)
    for(uint8_t i = 0 ; i<8 ; ++i)
     daughter[i]->   computeForcePairVerlet(pot, typeIndexA, typeIndexB) ;
}



template < class Pot_t>
void leafCell::computeForcePairVerlet(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB) {

  if(isleaf == 1)
  #ifndef _perBLOCK
    leafCell::computeForcePair_perso(pot, typeIndexA, typeIndexB) ;
   #else  
    computeForcePair_perso_per_block<Pot_t>(pot, typeIndexA, typeIndexB);
  #endif
  else 
    for(uint8_t i = 0 ; i<8 ; ++i)
     daughter[i]->   computeForcePairVerlet(pot, typeIndexA, typeIndexB) ;
	    
}

/// ONLY SYMETRIZED
template < class Pot_t, class CList, class>
void leafCell::computeForcePairVerlet_mutex(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter, const bool symmetrize) {

  if(isleaf == 1)
      leafCell::computeForcePairVerlet_mutex(pot, typeIndexA, typeIndexB) ;
  else 
    for(uint8_t i = 0 ; i<8 ; ++i)
     daughter[i]->   computeForcePairVerlet_mutex(pot, typeIndexA, typeIndexB, cells, ghostLayer, correcter, symmetrize) ;

}

/// @brief Function to compute forces for a EAM potential with SIMD optimized version
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
template < class Pot_t, class CList, class>
void leafCell::computeForceRhoVerlet(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter, const bool symmetrize) {

  if(isleaf == 1)
     computeForceRhoVerlet_2(pot, typeIndexA, typeIndexB, symmetrize) ;
  else 
    for(uint8_t i = 0 ; i<8 ; ++i)
     daughter[i]->   computeForceRhoVerlet(pot, typeIndexA, typeIndexB, cells, ghostLayer, correcter, symmetrize) ;
	    
}

/// @brief Function to compute forces for a EAM potential with SIMD optimized version
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
template < class Pot_t, class CList, class>
void leafCell::computeForceFinalVerlet(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, CList* cells, const uint8_t* ghostLayer, const Correcter& correcter, const bool symmetrize) {

  if(isleaf == 1)
     computeForceFinalVerlet_2(pot, typeIndexA, typeIndexB, symmetrize) ;
  else 
    for(uint8_t i = 0 ; i<8 ; ++i)
     daughter[i]->   computeForceFinalVerlet(pot, typeIndexA, typeIndexB, cells, ghostLayer, correcter, symmetrize) ;
	    
}


/// @brief Function to compute forces for a MEAM potential with SIMD optimized version. Don't forget to modify Verlet version.
/// between particles of two specified types
/// @tparam CList Class of child cellList doing the force computation
/// @param [in] pot Potential
/// @param [in] typeIndexA First type of particle
/// @param [in] typeIndexB Second type of particle
/// @param [in] cell the numeration of cell
/// @param [in] cells All the cells in the grid
/// @param [in] ghostLayer Ghost layer for each cell (cf Grid)
/// @param [in] correcter Tool to correct distances
template <bool mutex, class Pot_t, class CList>
void leafCell::computeForcesAndPotentialMEAMVerlet(Pot_t* pot, const uint8_t typeIndexA, const uint8_t typeIndexB, const int cell, CList* cells, const uint8_t * ghostLayer, const Correcter& correcter)
{
  if(isleaf == 1)
     leafCell::computeForcesAndPotentialMEAMVerlet_AMR<mutex>(pot, typeIndexA, typeIndexB, cell, ghostLayer);
  else 
    for(uint8_t i = 0 ; i<8 ; ++i)
     daughter[i]->   computeForcesAndPotentialMEAMVerlet<mutex>(pot, typeIndexA, typeIndexB, cell, cells, ghostLayer, correcter);

}
