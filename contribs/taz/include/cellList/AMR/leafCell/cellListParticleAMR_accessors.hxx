#pragma once

//accesors

inline uint8_t* leafCell::getType()           
{ 
  return granny->ti.data()+ shift; 
}

const inline uint32_t* leafCell::getMorton()    const    
{ 
  return granny->morton.data()+ shift; 
}

inline uint64_t* leafCell::getId()            
{ 
  return granny->id.data()+ shift; 
}

inline double* leafCell::getPotentialEnergy() 
{ 
  return granny->ep.data()+ shift; 
}

inline double* leafCell::getForceX()          
{ 
  return granny->fx.data()+ shift; 
}

inline double* leafCell::getForceY()          
{ 
  return granny->fy.data()+ shift; 
}

inline double* leafCell::getForceZ()          
{ 
  return granny->fz.data()+ shift; 
}

inline double* leafCell::getVelocityX()       
{ 
  return granny->vx.data()+ shift; 
}

inline double* leafCell::getVelocityY()       
{ 
  return granny->vy.data()+ shift; 
}

inline double* leafCell::getVelocityZ()       
{ 
  return granny->vz.data()+ shift; 
}

inline double* leafCell::getPositionX()       
{ 
  return granny->rx.data()+ shift; 
}

inline double* leafCell::getPositionY()       
{ 
  return granny->ry.data()+ shift; 
}

inline double* leafCell::getPositionZ()       
{ 
  return granny->rz.data()+ shift; 
}

inline vec3<int> leafCell::getOffset()       
{ 
  return granny->getOffset(); 
}

inline bool leafCell::getIsGhost()       
{ 
  return granny->getIsGhost(); 
}

inline bool leafCell::getIsEdge()       
{ 
  return granny->getIsEdge(); 
}

inline uint8_t leafCell::getType(const int i) const 
{ 
  return granny->ti[i+shift]; 
}

inline uint64_t leafCell::getId(const int i) const 
{ 
  return granny->id[i+shift]; 
}

inline double leafCell::getPotentialEnergy(const int i) const 
{ 
  return granny->ep[i+shift]; 
}

inline double leafCell::getForceX(const int i) const 
{ 
  return granny->fx[i+shift]; 
}

inline double leafCell::getForceY(const int i) const 
{ 
  return granny->fy[i+shift]; 
}

inline double leafCell::getForceZ(const int i) const 
{ 
  return granny->fz[i+shift]; 
}

inline double leafCell::getVelocityX(const int i) const 
{ 
  return granny->vx[i+shift]; 
}

inline double leafCell::getVelocityY(const int i) const 
{ 
  return granny->vy[i+shift]; 
}

inline double leafCell::getVelocityZ(const int i) const 
{ 
  return granny->vz[i+shift]; 
}

inline double leafCell::getPositionX(const int i) const 
{ 
  return granny->rx[i+shift]; 
}

inline double leafCell::getPositionY(const int i) const 
{ 
  return granny->ry[i+shift]; 
}

inline double leafCell::getPositionZ(const int i) const 
{ 
  return granny->rz[i+shift]; 
}

inline vec3<double> leafCell::getRealPosition(vec3<double>& sizeOfCell) 
{
  return position * sizeOfCell;
}

inline leafCell* leafCell::getDaughterCell (size_t i) { 
  assert(i<8);
  return daughter[i]; 
}

inline double& leafCell::embAMR(size_t i) 
{
  return granny->m_eamStorage.emb(i+shift);
}
 
inline double& leafCell::rhoAMR(size_t i) 
{
  return granny->m_eamStorage.rho(i+shift);
} 

inline int leafCell::getIsLeaf () 
{
  return isleaf; 
} 

inline uint leafCell::getNeighborCellsSize ()    
{ 
  return neighborCells.size(); 
}

inline leafCell* leafCell::getMotherCell () 
{ 
  return mother; 
}


inline int leafCell::getLevel () 
{ 
  return level; 
}    

inline uint leafCell::getSize () 
{ 
  return size; 
}   

inline vec3<int> leafCell::getPosition () 
{ 
  return position; 
}

inline vec3<int>& leafCell::getRefPosition () 
{ 
  return position; 
}

inline leafCell** leafCell::getNeighborCells () 
{ 
  return neighborCells.data(); 
}

inline size_t leafCell::getNumberOfParticles()
{
  return size;
}

// mutator

inline void leafCell::incForceX(const uint i, const double value) 
{ 
  granny->fx[i+shift]+=value; 
}

inline void leafCell::incForceY(const uint i, const double value) 
{ 
  granny->fy[i+shift]+=value; 
}

inline void leafCell::incForceZ(const uint i, const double value) 
{ 
  granny->fz[i+shift]+=value; 
}

inline void leafCell::incPontentialEnergy(const uint i, const double value) 
{ 
  granny->ep[i+shift]+=value; 
}


inline void leafCell::setIsLeaf (int i) 
{ 
  isleaf=i; 
}

inline void leafCell::setInfo(vec3<int> pos, int l) 
{
  position= pos; 
  level=l; 
} 

