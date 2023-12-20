#pragma once


inline void leafCell::defineNewCell(size_t begin, size_t numberOfElements, leafCell* motherCell, size_t posChild ) 
{
	// define infos;
	isleaf   = true;
	shift    = begin;
	size     = numberOfElements; 
	level    = motherCell->level +1; 

	/* poschild =(x={0,1}y={0,1}z={0,1})_bit */
	/* pos_daughter =  pos_mother + vec3_convert posChild */ 
	position = motherCell->position*2 + vec3<int>( (posChild>>2) & 1 , (posChild>>1) & 1 , posChild & 1 );

	/* define family information */
	granny     = motherCell->granny;
	mother     = motherCell;

	for(uint8_t i=0;i<8;i++) 
		daughter[i]=nullptr;
}

inline void leafCell::updateOldCell(size_t begin, size_t numberOfElements) 
{
	/* define information */;
	isleaf   = true;
	shift    = begin;
	size     = numberOfElements; 
}

inline void leafCell::createChild( const size_t refinement_max, size_t *nAtomPerCell)
{
	/* Computation of constants to extract information related to the refinement level of the cell. */
	/* Compute 2^{refinement max}. */
	const size_t totalCellPerDim = 1 << refinement_max;

	/* Compute 2^{refinement max-level}. */
	const size_t nCellPerDimPerLevel = totalCellPerDim >> level;

	/* Compute 8^{refinement max-level+1} */
	/* Corresponding to the number of cells in the floor of this octree. */
	const size_t nCellPerLevel = 1 << (3*(refinement_max-(level+1)));

	/* Computation of the morton index corresponding to the position of the cell in octree. */
	size_t a = encode( -granny->offset.x*totalCellPerDim + position.x*nCellPerDimPerLevel,
		-granny->offset.y*totalCellPerDim + position.y*nCellPerDimPerLevel,
		-granny->offset.z*totalCellPerDim + position.z*nCellPerDimPerLevel,
	0);

	/* Morton index always > 0*/
	assert(a>=0);

	/* Octree = 8 daughter cells */
	for(uint Child = 0 ; Child < 8 ; ++Child)
	{
		/* The daughter cell has never been allocated. */
		if(daughter[Child] == nullptr)
		{
			leafCell *tmp = new leafCell(nAtomPerCell[a], nAtomPerCell[a+nCellPerLevel]-nAtomPerCell[a], this, Child);
			daughter[Child] = tmp;
		}
		else 
			/* The daughter cell has already been allocated, its information is just updated. */
			daughter[Child]->updateOldCell( nAtomPerCell[a], nAtomPerCell[a+nCellPerLevel]-nAtomPerCell[a]);

		a+=nCellPerLevel;
	}

	isleaf=false;
}


template<typename F> inline void leafCell::refine(const size_t refinement_max, size_t* nAtomPerCell, F criterion)
{
	if(isleaf)
	{
		/* Refinement criteria */
		if(refinement_max > level && criterion(this)) 
		{
			/* Allocation of new daughter cells or update information they are already allocated. */
			createChild(refinement_max, nAtomPerCell);

			for(uint Child=0 ; Child<8 ; ++Child) 
			daughter[Child]->refine(refinement_max, nAtomPerCell, criterion);
		}
		else
			/* If it is a leaf cell, an allocation (arbitrary) is provided to store the n adjacent leaf cells. */
		{}	//neighborCells.resize(60); 
	}
	else 
		for(size_t Child=0 ; Child<8 ; ++Child)
			daughter[Child]->refine(refinement_max, nAtomPerCell, criterion);
};


inline void leafCell::adjust(const size_t refinement_max, size_t *nAtomPerCell) 
{


	/* Computation of constants to extract information related to the refinement level of the cell. */
	/* Compute 2^{refinement max}. */
	const size_t totalCellPerDim = 1 << refinement_max;

	/* Compute 2^{refinement max-level}. */
	const size_t nCellPerDimPerLevel = totalCellPerDim >> level;

	/* Compute 8^{refinement max-level} */
	/* Corresponding to the number of cells in the floor of this octree. */
	const size_t nCellPerLevel = 1 << (3*(refinement_max-(level)));

	/* Computation of the morton index corresponding to the position of the cell in octree. */
	size_t a = encode( -granny->offset.x*totalCellPerDim + position.x*nCellPerDimPerLevel,
		-granny->offset.y*totalCellPerDim + position.y*nCellPerDimPerLevel,
		-granny->offset.z*totalCellPerDim + position.z*nCellPerDimPerLevel,
	0);

	shift = nAtomPerCell[a];
	size  = nAtomPerCell[a+nCellPerLevel]-shift;

	/* Check if the number of atoms + shift is not out of range of the storage  */
	assert (shift + size  <=  granny->size);

	if(isleaf==false)
		for(size_t Child=0; Child<8; ++Child)
			daughter[Child]->adjust(refinement_max, nAtomPerCell);
  
};

