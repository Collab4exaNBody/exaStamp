#include <utils/morton.hpp>


inline void Octree::createChild( const size_t refinement_max, size_t *nAtomPerCell)
{
	/* Computation of constants to extract information related to the refinement level of the cell. */
	/* Compute 2^{refinement max}. */
	const size_t totalCellPerDim = 1 << refinement_max;

	/* Compute 2^{refinement max-level}. */
	const size_t nCellPerDimPerLevel = totalCellPerDim >> this->level;

	/* Compute 8^{refinement max-level+1} */
	/* Corresponding to the number of cells in the floor of this octree. */
	const int nCellPerLevel = 1 << (3*(refinement_max-(this->level+1)));

	/* Computation of the morton index corresponding to the position of the cell in octree. */
	size_t a = encode( -this->offset.x*totalCellPerDim + this->position.x*nCellPerDimPerLevel,
		-this->offset.y*totalCellPerDim + this->position.y*nCellPerDimPerLevel,
		-this->offset.z*totalCellPerDim + this->position.z*nCellPerDimPerLevel,
	0);

	/* Octree = 8 daughter cells */
	for(size_t Child = 0 ; Child < 8 ; ++Child)
	{
		/* The daughter cell has never been allocated. */
		if(this->daughter[Child] == nullptr)
		{
			this->daughter[Child] = new leafCell(nAtomPerCell[a],
				nAtomPerCell[a+nCellPerLevel]-nAtomPerCell[a], 
				static_cast<leafCell *>(this), 
			Child);
		}
		else 
		{
			/* The daughter cell has already been allocated, its information is just updated. */
			this->daughter[Child]->updateOldCell( nAtomPerCell[a], nAtomPerCell[a+nCellPerLevel]-nAtomPerCell[a]);
		}  
		a+=nCellPerLevel;
	}

	this->isleaf=false;
}


template<typename F> inline void Octree::refine(const size_t refinement_max, F criterion)
{
	if(refinement_max > this->level && criterion(this))
	{
		/* This is the root cell */
		this->granny = this;

		/* Maximal number of leaf cells (possible) per octree */
		const int numberOfCellMax = 1 << (3*refinement_max);

		/* Number of atoms per leaf cells if the octree was totaly refined. */
		size_t nAtomPerCell[numberOfCellMax+1];

		/* Intermediate storage to count atoms. */
		size_t iterator(0);

		/* First element contains 0 atoms */
		nAtomPerCell[0] = 0;

		for(size_t cell=1; cell<= numberOfCellMax ; ++cell)
		{
			/* All atoms included in the leaf cells (in a octree refined refinement_max) have the same morton index. */
			/* Atoms are previously sorted by their morton index */
			while(iterator < morton.size() && morton[iterator] < cell)
				iterator++;

			nAtomPerCell[cell] = iterator;

			/* The number of atoms included in the leaf cell 5 is obtained by the substraction nAtomPerCell[6]-nAtomPerCell[5]  */ 
		}

		assert(iterator >=0);

		/* Refinement criteria are OK, root cell is refined in 8 cells */
		/* Daughter cells are allocated only if daughter[Child] == nullptr, else they were already allocated  */
		createChild(refinement_max, nAtomPerCell);

		/* The refinement function is applied recursively until all leaf cells are found. */
		for(size_t Child = 0; Child<8 ; ++Child)
			this->daughter[Child]->refine(refinement_max, nAtomPerCell, criterion);

	}
	else {}
}



inline void Octree::adjust(const size_t refinement_max)
{

	assert(this->granny == this);

	/* No else case because the pointers point by default to the first element of the storage in the root cell. */ 
	if(this->isleaf== false)
	{	
		/* Maximal number of leaf cells (possible) per octree */
		const size_t numberOfCellMax = 1 << (3*refinement_max);

		/* Number of atoms per leaf cells if the octree was totaly refined. */
		size_t nAtomPerCell[numberOfCellMax+1];

		/* Intermediate storage to count atoms. */
		size_t iterator(0);

		/* First element contains 0 atoms */
		nAtomPerCell[0] = 0;

		for(size_t cell=1; cell<= numberOfCellMax ; ++cell)
		{
			/* All atoms included in the leaf cells (in a octree refined refinement_max) have the same morton index. */
			/* Atoms are previously sorted by their morton index */
			while(iterator < morton.size() && morton[iterator] < cell)
				iterator++;

			nAtomPerCell[cell] = iterator;

			/* The number of atoms included in the leaf cell 5 is obtained by the substraction nAtomPerCell[6]-nAtomPerCell[5]  */ 
		}

		assert(nAtomPerCell[numberOfCellMax] == this->size); 

		/* The adjust function is applied recursively until all leaf cells are adjusted */
		for(size_t Child = 0; Child<8 ; ++Child)
			this->daughter[Child]->adjust(refinement_max, nAtomPerCell);
  
	}
}


inline void Octree::fill_indx_morton( const vec3<double>& inv_sizeOfCell, const vec3<double> position)
{
	/* For each atom, a morton index, corresponding to a refined cell after n_max refinement, is associated with it. */
	for(size_t i=0; i<this->size; ++i)
	{
		assert(rx[i]>=position.x);
		assert(ry[i]>=position.y);
		assert(rz[i]>=position.z);

		assert(rx[i]< (position.x + Global::domainInfo.getCellLength().x ));
		assert(ry[i]< (position.y + Global::domainInfo.getCellLength().y ));
		assert(rz[i]< (position.z + Global::domainInfo.getCellLength().z ));

		morton[i] = encode( 
			size_t( (rx[i]-position.x)* inv_sizeOfCell.x ),
			size_t( (ry[i]-position.y)* inv_sizeOfCell.y ),
			size_t( (rz[i]-position.z)* inv_sizeOfCell.z ),
		0); 
	}
}
