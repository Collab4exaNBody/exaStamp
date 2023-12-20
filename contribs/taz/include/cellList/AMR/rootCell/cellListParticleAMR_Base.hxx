#pragma once



inline void Octree::resize(const uint& new_size) {

  resizeInfos(new_size);
  resizePositions(new_size);
  resizeVelocities(new_size);
  resizeForces(new_size);

  this->size=new_size;

  assert(this->size == id.size());
  assert(this->size == rx.size());
  
}

inline void Octree::resizeInfos(const uint& new_size) {

  id.resize(new_size);
  ti.resize(new_size);
  morton.resize(new_size);

}

inline void Octree::resizePositions(const uint& new_size) {

  rx.resize(new_size);
  ry.resize(new_size);
  rz.resize(new_size);
  
}
inline void Octree::resizeVelocities(const uint& new_size) {

  ep.resize(new_size);
  vx.resize(new_size);
  vy.resize(new_size);
  vz.resize(new_size);

}

inline void Octree::resizeForces(const uint& new_size) {

  fx.resize(new_size);
  fy.resize(new_size);
  fz.resize(new_size);

}


/// @brief Add a particle to the cell
///
/// Specialization for an MPI_ParticleBase type
/// @param [in] q Particle to add
inline void Octree::add(const MPI__ParticleBase& q) {

  this->size++;

  // Add all the existing elements to the corresponding arrays
  id.push_back(q.id);
  ti.push_back(q.ti);
  morton.push_back(0);
  ep.push_back(0);
  
  rx.push_back(q.r.x);
  ry.push_back(q.r.y);
  rz.push_back(q.r.z);
  
  // Initialize other elements
  vx.push_back(0);
  vy.push_back(0);
  vz.push_back(0); 

  fx.push_back(0);
  fy.push_back(0);
  fz.push_back(0);


  assert(this->size == id.size());
  assert(this->size == rx.size());
  // Increment the number of particles of its type
  ++this->m_numPerType[q.ti];

}


/// @brief Add a particle to the cell
///
/// Specialization for an MPI_Particle type
/// @param [in] q Particle to add
inline void Octree::add(const MPI__Particle& q) {

  this->size++;
  
  // Add all the existing elements to the corresponding arrays
  id.push_back(q.id);
  ti.push_back(q.ti);
  morton.push_back(0);
  ep.push_back(0);
  
  rx.push_back(q.r.x);
  ry.push_back(q.r.y);
  rz.push_back(q.r.z);
  
  // Initialize other elements
  vx.push_back(q.v.x);
  vy.push_back(q.v.y);
  vz.push_back(q.v.z); 
  
  fx.push_back(0);
  fy.push_back(0);
  fz.push_back(0);

 
  assert(this->size == id.size());
  assert(this->size == rx.size());

  // Increment the number of particles of its type
  ++this->m_numPerType[q.ti];

}


/// @brief Add a particle to the cell
///
/// Specialization for an MPI_Mesoparticle type
/// @param [in] q Particle to add
inline void Octree::add(const MPI__Mesoparticle& q) {


  this->size++;

  // Add all the existing elements to the corresponding arrays
  id.push_back(q.id);
  ti.push_back(q.ti);
  morton.push_back(0);
  ep.push_back(0);
  
  rx.push_back(q.r.x);
  ry.push_back(q.r.y);
  rz.push_back(q.r.z);
  
  vx.push_back(q.v.x);
  vy.push_back(q.v.y);
  vz.push_back(q.v.z);
  
  // Initialize other elements
  fx.push_back(0);
  fy.push_back(0);
  fz.push_back(0);

  // Increment the number of particles of its type
  ++this->m_numPerType[q.ti];

  assert(this->size == id.size());
  assert(this->size == rx.size());

}


// yet declare
struct passEraseElem { template<typename ...T> passEraseElem(T...) {} };


template<class... T> 
inline void eraseNElems (size_t e, std::vector<T>&... vec)
{
  passEraseElem{( assert(e < vec.size()) ,vec.erase(vec.begin() + e)   
    ,1)...};
}


/// @brief Remove a particle from the cell
/// @param [in] index Index of the particle
inline void Octree::remove(const uint& index) {

  // Decrement the number of particles of its type
  --this->m_numPerType[ ti[index] ];
  
  this->size--;
  
  eraseNElems(
    index, // elem to delete
    morton, id, ti, ep,
    rx, ry, rz, 
    vx, vy, vz,
    fx, fy, fz  
  );
  
  assert(this->size == id.size());
  assert(this->size == rx.size());

}



/// @brief Remove the last particle of the cell
///
///
inline void Octree::remove_last() {

  // Decrement the number of particles of its type
  --this->m_numPerType[ ti.back() ];
  
  this->size--;

  // Remove all elements of the particle
  id.pop_back();
  ti.pop_back();
  morton.pop_back();
  ep.pop_back();

  rx.pop_back();
  ry.pop_back();
  rz.pop_back();

  vx.pop_back();
  vy.pop_back();
  vz.pop_back();

  fx.pop_back();
  fy.pop_back();
  fz.pop_back();

  assert(this->size == id.size());
  assert(this->size == rx.size());

}


/// @brief Remove all the particles
///
///
inline void Octree::clear() {

	// Set all number of particles per type to 0
  for (uint8_t i=0; i<this->m_numPerType.size(); ++i) this->m_numPerType[i]=0;
  
  this->size=0;

  // Clear all arrays
  morton.clear();
  id.clear();
  ti.clear();
  ep.clear();

  rx.clear();
  ry.clear();
  rz.clear();

  vx.clear();
  vy.clear();
  vz.clear();

  fx.clear();
  fy.clear();
  fz.clear();

  assert(this->size == id.size());
  assert(this->size == rx.size());

}



inline void Octree::ghostAdjustAMR(const Correcter& correcter){

  vec3<double> p = (this->position-1) * Global::domainInfo.getCellLength() + Global::domainInfo.getMinBounds();
  correcter.correctP(this->size, p, rx.data(), ry.data(), rz.data());
  
}

/// @brief Debug print of the particles of the cell
/// @param [in,out] flux Print flux
inline void Octree::__debug_print(std::ostream& flux) const {

  using namespace std;

  const uint pr =  6;
  const uint sz = 16;

  // Sort particles by index
  std::vector<uint> sortId(this->size);

  for (uint j=0; j<this->size; ++j) 
    sortId[j] = j;

  std::sort(sortId.begin(), sortId.end(), [&](uint a, uint b) -> bool { return id[a]<id[b]; });

  // Print each particle
  for (uint i=0; i<this->size; ++i) {

    const uint& j = sortId[i];
    auto  nbrSize = this->neighborList.numberOfNeighbors(j);

    flux << setw(8) << 0 
	 << " | "
	 << setw(8) << id[j] 
	 << " | " 
	 << (int) ti[j] 
	 << " | " 
	 << scientific << setprecision(pr) << setw(sz) << ep[j] 
	 << " | "
	 << scientific << setprecision(pr) << setw(sz) << rx[j] << " "
	 << scientific << setprecision(pr) << setw(sz) << ry[j] << " "
	 << scientific << setprecision(pr) << setw(sz) << rz[j] 
	 << " | "
	 << scientific << setprecision(pr) << setw(sz) << vx[j] << " "
	 << scientific << setprecision(pr) << setw(sz) << vy[j] << " "
	 << scientific << setprecision(pr) << setw(sz) << vz[j] 
	 << " | "
	 << scientific << setprecision(pr) << setw(sz) << fx[j] << " "
	 << scientific << setprecision(pr) << setw(sz) << fy[j] << " "
	 << scientific << setprecision(pr) << setw(sz) << fz[j] 
	 << " | "
	 << setw(8) << nbrSize
	 << ""
	 << std::endl;
    
  }
  
}

/// @brief Get debug data from the cell (not used)
/// @return A tuple of : number of particles, capacity, mean number of neighbors, mean capacity of neighbor list
std::tuple<uint, uint, double, double> Octree::__debug_diag() const {

  uint num_particles = id.size();
  uint capacity = id.capacity();

  uint nbr_size=0, nbr_capa=0;

  for (uint i=0; i<num_particles; ++i)
  {
    nbr_size += this->neighborList.numberOfNeighbors(i);
    nbr_capa += this->neighborList.capacity(i);
  }

  return std::make_tuple(num_particles, capacity, (double)nbr_size/(double)num_particles, (double)nbr_capa/(double)num_particles);
} 



