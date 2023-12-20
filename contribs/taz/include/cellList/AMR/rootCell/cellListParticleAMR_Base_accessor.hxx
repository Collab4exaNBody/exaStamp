//accesors


inline bool Octree::getIsGhost() 
{
  return m_isGhost;
}

inline bool Octree::getIsReal()
{
  return !m_isGhost;
}

inline void Octree::setOffset(vec3<int>& o, bool ghost) 
{
  offset = o; 
  m_isGhost = ghost;
}

inline vec3<int> Octree::getOffset() 
{ 
  return offset;
}

inline void Octree::updateGranny() 
{ 
  this->granny = this;
}

int Octree::id_capacity() const             {return id.capacity();}
inline uint8_t* Octree::getType()           { return ti.data(); }
inline uint64_t* Octree::getId()            { return id.data(); }
inline double* Octree::getPotentialEnergy() { return ep.data(); }
inline double* Octree::getForceX()          { return fx.data(); }
inline double* Octree::getForceY()          { return fy.data(); }
inline double* Octree::getForceZ()          { return fz.data(); }
inline double* Octree::getVelocityX()       { return vx.data(); }
inline double* Octree::getVelocityY()       { return vy.data(); }
inline double* Octree::getVelocityZ()       { return vz.data(); }
inline double* Octree::getPositionX()       { return rx.data(); }
inline double* Octree::getPositionY()       { return ry.data(); }
inline double* Octree::getPositionZ()       { return rz.data(); }


inline uint8_t Octree::getType(const int i) const { return ti[i]; }
inline uint64_t Octree::getId(const int i) const { return id[i]; }
inline double Octree::getPotentialEnergy(const int i) const { return ep[i]; }
inline double Octree::getForceX(const int i) const { return fx[i]; }
inline double Octree::getForceY(const int i) const { return fy[i]; }
inline double Octree::getForceZ(const int i) const { return fz[i]; }
inline double Octree::getVelocityX(const int i) const { return vx[i]; }
inline double Octree::getVelocityY(const int i) const { return vy[i]; }
inline double Octree::getVelocityZ(const int i) const { return vz[i]; }
inline double Octree::getPositionX(const int i) const { return rx[i]; }
inline double Octree::getPositionY(const int i) const { return ry[i]; }
inline double Octree::getPositionZ(const int i) const { return rz[i]; }
