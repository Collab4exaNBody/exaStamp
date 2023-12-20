#include "utils/morton.hpp"
#include "cellList/AMR/sortTemplate.cpp"

inline void Octree::sort(const size_t refinement_max , const vec3<double>& sizeOfCell, const vec3<double> &inf)
{
	/* Inverse cell size after refinement_max refinement */
	const auto& inv_s =  (1<<refinement_max)/sizeOfCell;
  
	/* Get the morton index for each atom */
	fill_indx_morton( inv_s, (this->position-1)*sizeOfCell+inf);
      
	/* Atoms are sorted by their morton index */
	introsort(this->size, morton.data(),	// array to sort
		id.data(), ti.data(),		//arrays modified
		rx.data(), ry.data(), rz.data(),
		vx.data(), vy.data(), vz.data());
}



