#pragma once

#include <grid/AMRGrid/externAMRGrid.hpp>


/// @file
/// @brief Implementation of the neighbors search related methods in class AMRGrid

/// @brief Compute the neighbor list on all the cells (except last layer of ghost)
///
/// Call a specialized function depending on the ghost thickness
TMPLSG inline void TMPL_AMRGrid::makeNeighborLists() {

	if(Global::reference.isMEAM()) __makeNeighborLists_Traversal_MEAM_AMR(TraversalManager::ALL);
	else __makeNeighborLists_gt1_AMR() ;


	if(Global::reference.isMEAM()) cWave.resortAllTraversal( this->tManager, particles.data(), info, numberOfCellsPerDim, false);
	else cWave.resort(TraversalManager::REAL, this->tManager, particles.data(), info, numberOfCellsPerDim, false);
}


/// @brief Compute the neighbor list on the inside cells
///
/// Call a specialized function depending on the ghost thickness
TMPLSG inline void TMPL_AMRGrid::makeNeighborListsInside() {

	__makeNeighborLists_Traversal_AMR(TraversalManager::INSIDE_GTHICK);

	cWaveEdge.resort( TraversalManager::EDGE_ONE, this->tManager, particles.data(), info, numberOfCellsPerDim,  false);
}


/// @brief Compute the neighbor list on edge cells
///
/// Call a specialized function depending on the ghost thickness
TMPLSG inline void TMPL_AMRGrid::makeNeighborListsOnEdges() {
  
	__makeNeighborLists_Traversal_AMR(TraversalManager::EDGE_ONE);

	cWaveEdge.resort( TraversalManager::EDGE_ONE, this->tManager, particles.data(), info, numberOfCellsPerDim, false);
  
}


/// @brief Clear neighbor lists
///
///
TMPLSG void TMPL_AMRGrid::clearNeighborList() {

	#pragma omp parallel for schedule(runtime)
  for(size_t i = 0 ; i < octrees.size() ; i++)
			octrees[i]->clearNeighborList();
	

	if(Global::reference.isMEAM()) {
		const auto& cells = this->getTraversal(TraversalManager::GHOST);

		// Parallel loop on traversal cells
		parallel_region(0, cells.size(),
		[&](const uint begin, const uint end) {
			for (uint i=begin; i!=end; ++i)
				particles[cells[i]].clearNeighborList() ;
		});
	}
}


/// @brief Sort the neighbor lists so that so that particles
/// of the same type are contiguous in the lists
/// @param [in] traversal Identify the cells to sort
TMPLSG void TMPL_AMRGrid::sortNeighbors(::Traversal traversal) {


  #pragma omp parallel for schedule(runtime)
  for(size_t i = 0 ; i < octrees.size() ; i++)
			octrees[i]->sortNeighbors();

}



#define Cell leafCell

/// @brief Compute the neighbor list specified traversal for a ghost thickness of 1
/// @tparam elem_chunk Base chunk for cell arrays
/// @param [in] t Identify the cells to work on
TMPLSG void TMPL_AMRGrid::__makeNeighborLists_gt1_AMR() {

	  #pragma omp parallel for schedule(runtime)
	  for(size_t i = 0 ; i < octrees.size() ; i++)
	    octrees[i]->fillVerlet();


	if(Global::reference.isBlockVerlet())
	{
		#pragma omp parallel for schedule(runtime)
		for (size_t i=0; i<leafcells.size(); ++i)
		{
			if(leafcells[i]->size != 0)
				__makeNeighborLists_gt1_AMR<true>(leafcells[i]) ;
		} 
	}
	else
	{
		
		#pragma omp parallel for schedule(runtime)
		for (size_t i=0; i<leafcells.size(); ++i)
		{
			if(leafcells[i]->size != 0)
				__makeNeighborLists_gt1_AMR<false>(leafcells[i]) ;
		} 
	}
	const auto& cells = this->getTraversal(TraversalManager::REAL);


	#pragma omp parallel for schedule(runtime)
	for(size_t i = 0 ; i < cells.size() ; i++)
		particles[cells[i]].checkMtx();
}

TMPLSG  template<typename F>
void TMPL_AMRGrid::__RecApply(F function, Cell *C) {

	if(C->size == 0 ) return ;

	if(C->isleaf == 1) function(C) ;
	else 
		for(int i = 0; i<8; i++) __RecApply(function,C->getDaughterCell(i));

}


/// @brief Compute the neighbor list specified traversal for a ghost thickness of 1
/// @tparam elem_chunk Base chunk for cell arrays
/// @param [in] t Identify the cells to work on
TMPLSG void TMPL_AMRGrid::__makeNeighborLists_Traversal_AMR(TraversalManager::CellTraversal t) {


	const auto& cells = this->getTraversal(t);

	if(Global::reference.isBlockVerlet())
	{
		auto fun = [&] (Cell *C)->void {__makeNeighborLists_gt1_AMR<true>(C);};
		#pragma omp parallel for schedule(runtime)
		for(int i = 0 ; i < cells.size() ; i++)
		{
		  particles[cells[i]].fillVerlet();
			__RecApply(fun, &particles[cells[i]]);
			particles[cells[i]].checkMtx();
		}
	}
	else 
	{
		auto fun = [&] (Cell *C)->void {__makeNeighborLists_gt1_AMR<false>(C);};
		#pragma omp parallel for schedule(runtime)
		for(int i = 0 ; i < cells.size() ; i++)
		{
		  particles[cells[i]].fillVerlet();
			__RecApply(fun, &particles[cells[i]]);
			particles[cells[i]].checkMtx();
		}
	}
}


/// @brief Compute the neighbor list specified traversal for a ghost thickness of 1
/// @tparam elem_chunk Base chunk for cell arrays
/// @param [in] t Identify the cells to work on
TMPLSG void TMPL_AMRGrid::__makeNeighborLists_Traversal_MEAM_AMR(TraversalManager::CellTraversal t) {

	checkLocalBuffers();
	checkLocalBuffersGHOST(); //ALL

	
	#pragma omp parallel for schedule(runtime) // Pas besoin de stocker pour les ghosts octrees
	for(size_t i = 0 ; i < octreesMEAM.size() ; i++)
		particles[octreesMEAM[i]].fillVerlet();

	#pragma omp parallel for schedule(runtime)
	for(int i = 0 ; i < leafcellsMEAM.size() ; i++)
	{
		if(leafcellsMEAM[i]->size != 0)
		__makeNeighborLists_MEAM_AMR(leafcellsMEAM[i]);
	}

	#pragma omp parallel for schedule(runtime)
	for(int i = 0 ; i < octreesMEAM.size() ; i++)
	{
		particles[octreesMEAM[i]].checkMtx();
	}
}

/// @brief Compute the neighbor list specified traversal for a ghost thickness of 1
/// @tparam elem_chunk Base chunk for cell arrays
/// @param [in] C Identify the cells to work on
TMPLSG template< bool perBlock>
void TMPL_AMRGrid::__makeNeighborLists_gt1_AMR(Cell* C) {


	vec3<int> dcell;  /* Determine the direction of the neighbor cell. Used to compute the masks */
	double iter;      /* To store intermediary result */  
	int other_lvl;    /* Level of refinement of the neighbor cell */

	/* Compute an index to determine if the neighbour cell have yet done (symmetrical case) */
	/* This index is unique for one couple of position and level                            */
	const vec3<int> tmp = vec3<int>(1,numberOfCellsPerDim.x+2,(numberOfCellsPerDim.x+2)*(numberOfCellsPerDim.y+2));
	int tmp1=1<<C->getLevel();
	int cell = dot(C->getPosition(),tmp)+info.getTotalNumberOfCells()*(1<<3*C->getLevel()); 
	const bool symmetrize = this->domain->symmetrize();

	/* detect if atoms are too far away to interact with at least one atom of the neighbour cell */
	bool projMask[C->getNumberOfParticles()];

	/* Retrieve the list of the neighbor cells and their number */
	const auto& List = C->getNeighborCells();
	const int  size = C->getNeighborCellsSize();

	/* Retrieve information about the cell C */
	const auto& coord = C->getPosition() - tmp1 +1;
	const auto& self_position = C->getPosition();
	const int self_lvl = C->getLevel();

	/* We begin by fill the atom neighbour lists with interaction with atoms included the same cell */
	/* Note that the fist parameter of the makeNeighboorsAMR* consist of to store the index of the neighbor */
	/* cell where the neighbor atom is stored. size-1 = himself */
	C->makeNeighborsAMR_self(size-1, symmetrize, perBlock);

	/* Loop on all neighbor cells */
	for (size_t nbr=0; nbr<size-1; ++nbr) {

		auto& neighbor = *List[nbr];
		other_lvl = neighbor.getLevel();

		/* Some debug check */
		assert(List[nbr] != C);
		assert( !((self_lvl == other_lvl) && (C->position == neighbor.position)) );

		/* Only used for symmetrize case */
		const int nbrCell = dot(neighbor.getPosition(),tmp)+info.getTotalNumberOfCells()*(1<<(3*other_lvl));

		/* Determine if the neighbor cell is a ghost or not */
		const bool isReal = !neighbor.getIsGhost();
		const bool do_symmetrize = ((symmetrize && cell<nbrCell) || !symmetrize);

		iter  = self_lvl >= other_lvl ? std::pow(0.5, self_lvl-other_lvl ) : std::pow(0.5, other_lvl - self_lvl );
		vec3<double> tmp3 = neighbor.getPosition()*iter;
		vec3<double> tmp4 = self_position*iter;   
		dcell = self_lvl >= other_lvl ? neighbor.getPosition()-  auxFloor(tmp4): auxFloor(tmp3)-self_position;

		// Set projection mask
		C->getProjAMR(coord, dcell, projMask);

		
		if (isReal) // All real case
			// Find the neighbors in neighbor cells
			C->makeNeighborsAMR (nbr, neighbor, projMask, do_symmetrize, perBlock);
		else  	// Edge case
			// Find the neighbors in neighbor cells
			C->makeNeighborsAMR (nbr, neighbor, projMask, true, perBlock);

	} // End loop on neighbor cells
}


/// @brief Compute the neighbor list specified traversal for a ghost thickness of 1
/// @tparam elem_chunk Base chunk for cell arrays
/// @param [in] t Identify the cells to work on
TMPLSG
void TMPL_AMRGrid::__makeNeighborLists_MEAM_AMR(Cell* C) {
		
	vec3<int> tmp = vec3<int>(1,numberOfCellsPerDim.x+2,(numberOfCellsPerDim.x+2)*(numberOfCellsPerDim.y+2));

	// Find the neighbors in present cell
	int tmp1=1<<C->getLevel();
	//int cell = dot(C->getPosition(),tmp)+info.getTotalNumberOfCells()*tmp1;
	const bool symmetrize = this->domain->symmetrize();

	/* detect if atoms are too far away to interact with at least one atom of the neighbour cell */
	bool projMask[C->getNumberOfParticles()];

	const auto& List = C->getNeighborCells();
	const size_t  size = C->getNeighborCellsSize();
	const vec3<int> coord = C->getPosition() - tmp1 +1;
	const vec3<int> self_position =C->getPosition() - tmp1;
	vec3<int> dcell;  
	double iter;
	const int self_lvl = C->getLevel();
	assert( size > 0 );
	C->makeNeighborsAMR_self (size-1, symmetrize, false);


	// Loop on all neighbors
	for (uint nbr=0; nbr<(size-1); ++nbr)
	{
		assert( List[nbr] != C );

		auto& neighbor = *List[nbr];
		int other_lvl = neighbor.getLevel();
		vec3<int> other_position = neighbor.position - (1<<other_lvl);

		iter  = self_lvl >= other_lvl ? std::pow(0.5, self_lvl - other_lvl ) 
			: std::pow(0.5, other_lvl - self_lvl );

		vec3<double> tmp3 = other_position*iter;
		vec3<double> tmp4 = self_position*iter;   

		dcell = self_lvl >= other_lvl ? other_position - vec3<int>(floor(tmp4.x),floor(tmp4.y),floor(tmp4.z))
			: vec3<int>(floor(tmp3.x),floor(tmp3.y),floor(tmp3.z))-self_position;

		// Set projection mask
		C->getProjAMR(coord, dcell, projMask);

		// All real case
		// Find the neighbors in neighbor cells
		C->makeNeighborsAMR (nbr, neighbor, projMask, true, false);
	} // End loop on neighbor cells

}


 /*!
  *  \brief return last 3bits for a level
  * 
  *  \param b : locational code
  *  \param e : level of the cell
  *  \return 3bits in 3 integer
  */
TMPLSG inline vec3<int> TMPL_AMRGrid::get_Vec3_Val(Cell *C)
{
	return C->getPosition() - 2*C->getMotherCell()->getPosition();
}
  
/*!
*  \brief We check if all field of a is inferior to b
* 
*  \param a : vec3<int> wich contains cartesian position
*  \param b : vec3<int> wich contains cartesian position
*  \return bool if all field of a is inferior to b
*/
TMPLSG inline bool TMPL_AMRGrid::testInf(vec3<int> a, vec3<int> b)
{
	if(a.x >= b.x) return false;
	else if(a.y >= b.y) return false;
	else if(a.z >= b.z) return false;
	else return true;
}

/*!
*  \brief We check if all field of a is superior or equal to 0
* 
*  \param a : vec3<int> wich contains cartesian position
*  \return bool if all field of a is superior or equal to 0
*/
TMPLSG inline bool TMPL_AMRGrid::testSup1(vec3<int> a)
{
	if(a.x<1) return false;
	else if(a.y<1) return false;
	else if(a.z<1) return false;
	else return true;
}



TMPLSG void TMPL_AMRGrid::fillneighborCells_List()
{

	parallel_region(0, leafcells.size(),
	[&](const uint begin, const uint end) {

		for (uint i=begin; i!=end; ++i) {
			fillNeighborsCells_list(leafcells[i]) ;
			leafcells[i]->avoidDoublon(numberOfCellsPerDim, info.getTotalNumberOfCells());
			leafcells[i]->addNeighborCell(leafcells[i]);
		}
	});

	if(Global::reference.isMEAM()) {

		const auto& cells = this->getTraversal(TraversalManager::GHOST);

		parallel_region(0, cells.size(),
		[&](const uint begin, const uint end) {

			for (uint i=begin; i!=end; ++i) {
				fillNeighborsCells_list_GHOST(&particles[cells[i]]) ;
				particles[cells[i]].avoidDoublon(numberOfCellsPerDim, info.getTotalNumberOfCells());
				particles[cells[i]].addNeighborCell(&particles[cells[i]]);
			}
		});
	}
} 



TMPLSG void TMPL_AMRGrid::fillNeighborsCells_list_GHOST(Cell *C) 
{
	C->clearNeighborCells();

	vec3<int> offsetCellNeighbor = C->position; //- info.getCellInf();
	vec3<int> tmp;
	for(int i=-1; i<2 ; i++)
		for(int j=-1; j<2 ; j++)
			for(int k=-1; k<2 ; k++)
			{
				if( (i==0) && (j==0) && (k==0) ) continue;

				tmp= offsetCellNeighbor+vec3<int>(i,j,k);  

				int check = info.indexForAMRWithGhost(tmp);

				if(check != -1)
					recursiveCellSearch_level_0( &particles[check], C); 
			}
}

/*!
*  \brief fill array of neighbor cells of the cell C->
*
*  \param C : Cell
*  \return nothing
*/
TMPLSG void TMPL_AMRGrid::fillNeighborsCells_list(Cell *C) 
{
	C->clearNeighborCells();

	if(C->getMotherCell() !=NULL)
	{
		vec3<int> localPosition = get_Vec3_Val(C); // localPosition = local position in the octree
		int localMorton=localPosition.x*4+localPosition.y*2+localPosition.z;
		cellSearch_1(localMorton, C);
		uint32_t biais;

		// Find in the 7 other nodes
		vec3<int> offsetCellNeighbor, posNeighbor;
		vec3<double> temp;
		int CellSize = (1 << (C->getLevel()-1));
		double InvCellSize = 1./CellSize;
		vec3<double> localPosition2 = vec3<int>( (localPosition.x==1 ? 1:-1),(localPosition.y==1 ? 1:-1), (localPosition.z==1 ? 1:-1));

		for(int i=0; i<2 ; i++)
			for(int j=0; j<2 ; j++)
				for(int k=0; k<2 ; k++)
				{
					if( (i==0) && (j==0) && (k==0) ) continue;

					posNeighbor  = C->getMotherCell()->getPosition();
					posNeighbor += localPosition2 *(vec3<int>(i,j,k));                            
					temp         = vec3<double>(posNeighbor) * InvCellSize;                            
					offsetCellNeighbor  = vec3<int>(floor(temp.x),floor(temp.y),floor(temp.z));
					posNeighbor  = posNeighbor-offsetCellNeighbor*CellSize;                                          
					biais        = encode(posNeighbor.x, posNeighbor.y, posNeighbor.z, levelMax-C->getLevel()+1);
					int index = info.indexForAMRWithGhost(offsetCellNeighbor);
					if(index != -1)
						CellSearch(biais, &particles[index], C); 
				}
	}
	else 
	{
		vec3<int> offsetCellNeighbor = C->getPosition();
		vec3<int> temp;
		for(int i=-1; i<2 ; i++)
			for(int j=-1; j<2 ; j++)
				for(int k=-1; k<2 ; k++)
				{
					if( (i==0) && (j==0) && (k==0) ) continue;

					temp=vec3<int>(i,j,k);  
					int index = info.indexForAMRWithGhost(offsetCellNeighbor+temp);
					if(index != -1)
						recursiveCellSearch_level_0( &particles[index], C); 
				}
	}
}





/*!
*  \brief the first step to find neighbor cells in the other nodes.
*  
*  \param neighborRoot : path to access to the mother neighbor cell in the root.
*  \param root : root of the octree wich contains mother neighbor cell. Root is defined by offset in other function
*  \param r2 : square raduiscutoff of the potential.
*  \return nothing
*/
TMPLSG void TMPL_AMRGrid::CellSearch(uint32_t neighborRoot,  Octree* root, Cell* C)
{


	assert(root != C);

	uint32_t e0=3*(levelMax - root->getLevel());
	Cell* ptr_root = static_cast<Cell*> (root);

	for(uint16_t i=0; (i<C->getLevel()-1) && (ptr_root->getIsLeaf() == 0); ++i)
	{
		e0  -= 3;
		ptr_root = ptr_root->getDaughterCell(  ((neighborRoot >> e0) & 1) + (( (neighborRoot >> (e0+1) ) & 1 ) << 1 ) + (( (neighborRoot >> (e0+2) ) & 1 ) << 2) ); 
	}

	if(ptr_root->getIsLeaf()==1)
		C->addNeighborCell(ptr_root); 
	else
		for(uint8_t j=0;j<=7;j++)
			recursiveCellSearch( ptr_root->getDaughterCell(j), C);
}


/*!
*  \brief the second step to find neighbor cells in the other nodes (recursive search).
*  \param root : neighbors Cell but maynot be a leaf cell. //bad name
*  \param C : this is the cell where we want to build neighbor lists.
*  \return nothing
*/
TMPLSG void TMPL_AMRGrid::recursiveCellSearch_level_0(Cell *root, Cell *C)
{
  vec3<int> p1, p2;
  double temp;


	assert(root != C);

  p1 = root->getPosition();
  p2 = C->getPosition();

  //debug
  p1 = p1 * std::pow(2,C->getLevel());
 // p2 = p2 * std::pow(2,C->getLevel());

  if( (fabs(p1.x-p2.x) < 2 ) && ( fabs(p1.y-p2.y) < 2 ) && (fabs(p1.z-p2.z) < 2 ))
    {
    if(root->getIsLeaf()==1)  C->addNeighborCell(root);
    else
    {
      p2=2*p2; // warning if u must debug
      p1=2*p1;
      for(int bitX=0;bitX<2;bitX++)
      {
        temp=p1.x+bitX-p2.x;
        if( ( temp > -2 ) && ( temp < 3))
        {
          for(int bitY=0;bitY<2;bitY++)
          {
            temp = p1.y+bitY-p2.y;
            if(( temp > -2 ) && ( temp < 3) )
            {
              for(int bitZ=0;bitZ<2;bitZ++)
              {
                temp = p1.z+bitZ-p2.z;
                if( ( temp > -2 ) && ( temp < 3) ) 
                {
                  int tmp = bitX*4 + bitY *2 + bitZ;
                  if(root->getDaughterCell(tmp)->getIsLeaf() == 1) C->addNeighborCell(root->getDaughterCell(tmp) ); 
                  else getDifferentLevelNeighborsCellInOtherNode( root->getDaughterCell(tmp), C);
                }
              }
            }
          }
        }
      }
    }
  }

}


/*!
*  \brief the second step to find neighbor cells in the other nodes (recursive search).
*  \param root : neighbors Cell but maynot be a leaf cell. //bad name
*  \param C : this is the cell where we want to build neighbor lists.
*  \return nothing
*/
TMPLSG void TMPL_AMRGrid::recursiveCellSearch(Cell *root, Cell *C)
{
  vec3<int> p1, p2;
  double temp;


	assert(root != C);

  p1 = root->getPosition();
  p2 = C->getPosition();

  if( (fabs(p1.x-p2.x) < 2 ) && ( fabs(p1.y-p2.y) < 2 ) && (fabs(p1.z-p2.z) < 2 ))
  {
    if(root->getIsLeaf()==1)  C->addNeighborCell(root);
    else
    {
      p2=2*p2; // warning if u must debug
      p1=2*p1;
      for(int bitX=0;bitX<2;bitX++)
      {
        temp=p1.x+bitX-p2.x;
        if( ( temp > -2 ) && ( temp < 3))
        {
          for(int bitY=0;bitY<2;bitY++)
          {
            temp = p1.y+bitY-p2.y;
            if(( temp > -2 ) && ( temp < 3) )
            {
              for(int bitZ=0;bitZ<2;bitZ++)
              {
                temp = p1.z+bitZ-p2.z;
                if( ( temp > -2 ) && ( temp < 3) ) 
                {
                  int tmp = bitX*4 + bitY *2 + bitZ;
                  if(root->getDaughterCell(tmp)->getIsLeaf() == 1) C->addNeighborCell(root->getDaughterCell(tmp) ); 
                  else getDifferentLevelNeighborsCellInOtherNode( root->getDaughterCell(tmp), C);
                }
              }
            }
          }
        }
      }
    }
  }

}

/*!
*  \brief the third step to find neighbor cells in the other nodes (recursive search).
*  \param root : neighbors Cell but maynot be a leaf cell. //bad name
*  \param C : this is the cell where we want to build neighbor lists.
*  \param n : this is the cell where we want to build neighbor lists.
*  \return nothing
*/
TMPLSG void TMPL_AMRGrid::getDifferentLevelNeighborsCellInOtherNode(Cell *root, Cell *C)
{
	assert(root != C);

  vec3<int> p1, p2;
  double temp;
  
  p1 = root->getPosition();
  p2 = C->getPosition();

  int n; 

  if(C->getLevel() > root->getLevel()) 
  {
    n = C->getLevel() - root->getLevel();
    p1 = p1*(1<<n);
  } 
  else
  {
    n= root->getLevel()-C->getLevel();
    p2 = p2*(1<<n);
  }

  n*=2;

  for(int bitX=0;bitX<2;bitX++)
  {
    temp=p1.x+bitX-p2.x;
    if( ( temp > -2 ) && ( temp < n+1))
    {
      for(int bitY=0;bitY<2;bitY++)
      {
        temp = p1.y+bitY-p2.y;
        if(( temp > -2 ) && ( temp < n+1) )
        {
          for(int bitZ=0;bitZ<2;bitZ++)
          {
            temp = p1.z+bitZ-p2.z;
            if( ( temp > -2 ) && ( temp < n+1) ) 
            {
              int tmp = bitX*4 + bitY *2 + bitZ;
              if(root->getDaughterCell(tmp)->getIsLeaf()==true) C->addNeighborCell(root->getDaughterCell(tmp) ); 
              else getDifferentLevelNeighborsCellInOtherNode( root->getDaughterCell(tmp), C);
            }
          }
        }
      }
    }
  }
}


/*!
*  \brief the first step to find neighbor cells in the same node is to find each cell of the same level
* 
*  \param noNeighbor : local index in the node of the cell, used to avoid built neighbor list again.
*  \param C : this is the cell wich want neighbor lists.
*  \param r2 : square raduiscutoff of the potential.
*  \return 3bits in 3 integer
*
/*MPLSG void TMPL_AMRGrid::cellSearch_1(int noNeighbor, cell &C, cell &mother)
{
	for(int j=0; j< noNeighbor ;j++)
		recursiveCellSeach_1(mother.getDaughterCell(j), C);

	for(int j= noNeighbor+1 ; j<8  ;j++)
		recursiveCellSeach_1(mother.getDaughterCell(j), C);
}*/


/*!
*  \brief the first step to find neighbor cells in the same node is to find each cell of the same level
* 
*  \param noNeighbor : local index in the node of the cell, used to avoid built neighbor list again.
*  \param C : this is the cell wich want neighbor lists.
*  \param r2 : square raduiscutoff of the potential.
*  \return 3bits in 3 integer
*/
TMPLSG void TMPL_AMRGrid::cellSearch_1(int noNeighbor, Cell *C)
{
	for(int j=0; j< noNeighbor ;j++)
		recursiveCellSeach_1(C->getMotherCell()->getDaughterCell(j), C);

	for(int j= noNeighbor+1 ; j<8  ;j++)
		recursiveCellSeach_1(C->getMotherCell()->getDaughterCell(j), C);
}


/*!
*  \brief the second step to find neighbor cells in the same node is to find cells wich are leaf
* 
*  \param noNeighbor : local index in the node of the cell, used to avoid built neighbor list again.
*  \param C : this is the cell wich want neighbor lists.
*  \param r2 : square raduiscutoff of the potential.
*  \return 3bits in 3 integer
*/
/*TMPLSG void TMPL_AMRGrid::recursiveCellSeach_1(cell &root, cell &C)
{
  vec3<int> p1,p2,ptemp;

  if(root.isLeaf == true) leafCell.addNeighborCell();
  else
  {
    p1 = root->getPosition();
    p2 = C->getPosition();

    int temp;

    p2=2*p2; // warning if u must debug
    p1=2*p1;
    for(int bitX=0;bitX<2;bitX++)
    {
      temp=p1.x+bitX-p2.x;
      if( ( temp > -2 ) && ( temp < 3))
      {
        for(int bitY=0;bitY<2;bitY++)
        {
          temp = p1.y+bitY-p2.y;
          if(( temp > -2 ) && ( temp < 3) )
          {
            for(int bitZ=0;bitZ<2;bitZ++)
            {
              temp = p1.z+bitZ-p2.z;
              if( ( temp > -2 ) && ( temp < 3) ) 
              {
                int tmp = bitX*4 + bitY *2 + bitZ;
                if(root->getDaughterCell(tmp)->getIsLeaf()==1) C->addNeighborCell(root->getDaughterCell(tmp) ); 
                else getDifferentLevelNeighborsCellInOtherNode( root->getDaughterCell(tmp), C);
              }
            }
          }
        }
      }
    }
  }
}*/


/*!
*  \brief the second step to find neighbor cells in the same node is to find cells wich are leaf
* 
*  \param noNeighbor : local index in the node of the cell, used to avoid built neighbor list again.
*  \param C : this is the cell wich want neighbor lists.
*  \param r2 : square raduiscutoff of the potential.
*  \return 3bits in 3 integer
*/
TMPLSG void TMPL_AMRGrid::recursiveCellSeach_1(Cell *root, Cell *C)
{

	assert(root != C);

  vec3<int> p1,p2,ptemp;

  if(root->getIsLeaf()==1) C->addNeighborCell(root);
  else
  {
    p1 = root->getPosition();
    p2 = C->getPosition();

    int temp;

    p2=2*p2; // warning if u must debug
    p1=2*p1;
    for(int bitX=0;bitX<2;bitX++)
    {
      temp=p1.x+bitX-p2.x;
      if( ( temp > -2 ) && ( temp < 3))
      {
        for(int bitY=0;bitY<2;bitY++)
        {
          temp = p1.y+bitY-p2.y;
          if(( temp > -2 ) && ( temp < 3) )
          {
            for(int bitZ=0;bitZ<2;bitZ++)
            {
              temp = p1.z+bitZ-p2.z;
              if( ( temp > -2 ) && ( temp < 3) ) 
              {
                int tmp = bitX*4 + bitY *2 + bitZ;
                if(root->getDaughterCell(tmp)->getIsLeaf()==1) C->addNeighborCell(root->getDaughterCell(tmp) ); 
                else getDifferentLevelNeighborsCellInOtherNode( root->getDaughterCell(tmp), C);
              }
            }
          }
        }
      }
    }
  }
}
 
