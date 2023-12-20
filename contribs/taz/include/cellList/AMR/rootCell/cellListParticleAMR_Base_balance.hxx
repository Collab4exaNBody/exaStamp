#pragma once


/// @brief Compute the memory used by the cell
/// @return Memory usage
double Octree::computeMemory() const {
  return (double) (sizeof(*this) + id_capacity() * ( sizeof(uint64_t)+sizeof(uint8_t)+10*sizeof(double) ));
}

/// @brief Compute the workload on the cell
/// @return Workload
inline double Octree::computeWorkloadAMR() const {

  double w=0;

  // On calcul la charge de cet octree durant la dernière mise à jour des listes de Verlet
  size_t n = this->neighborList.nbElements();
  
  for(size_t i = 0 ; i < n ; i++)
      w += this->neighborList.numberOfNeighbors(i+shift)*this->neighborList.numberOfNeighbors(i+shift);

  return w;
 
}
