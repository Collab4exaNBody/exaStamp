/** 
 * @file
 * @brief Definition of the shared structure between the molecular dynamics code and the analytics library (for in situ processing)
 * 
 */


#ifndef __PARTICLE_INSITU_HPP_INCLUDED
#define __PARTICLE_INSITU_HPP_INCLUDED


#include <unordered_map>


struct ParticleInSitu {
  
public:

  ///@brief Default constructor
  ParticleInSitu() :
    nPart(0), step(-1), ghostAllocated(false), qAllocated(false)
  {
    for (int i=0; i<3; ++i)
      extension[i] = 0.0;
  }

  /// @brief Constructor with the number of particles
  /// @param[in] n Number of particles
  ParticleInSitu(uint n) :
    nPart(n), step(-1), ghostAllocated(false), qAllocated(false)
  {
    for (int i=0; i<3; ++i)
      extension[i] = 0.0;
    id = new uint[n];
    type = new uint8_t[n];
    rx = new double[n];
    ry = new double[n];
    rz = new double[n];
    vx = new double[n];
    vy = new double[n];
    vz = new double[n];
  }

  /// @brief Destructor
  ~ParticleInSitu()
  {
    delete id;
    delete type;
    delete rx;
    delete ry;
    delete rz;
    delete vx;
    delete vy;
    delete vz;

    if (ghostAllocated)
      {
        delete ghost_id;
        delete ghost_rx;
        delete ghost_ry;
        delete ghost_rz;
        ghostAllocated = false;
      }

    if (qAllocated)
      {
        for (uint i=0; i<this->qLen; ++i)
          delete q[i];
        delete q;
        qAllocated = false;
      }
    
  }

  /// @brief Fix the number of ghosts and allocate arrays
  /// @param[in] _nGhost The number of ghost particles
  void setNbGhost(uint _nGhost)
  {
    this->nGhost = _nGhost;
    if ( this->nGhost > 0 )
      {
        ghost_id = new uint[nGhost];
        ghost_rx = new double[nGhost];
        ghost_ry = new double[nGhost];
        ghost_rz = new double[nGhost];
        ghostAllocated = true;
      }
  }

  /// @brief Set the number of q parameters to be computed and allocate corresponding arrays
  /// @param[in] qlen Number of q parameters to be computed
  void allocateQ(uint qlen)
  {
    this->qLen = qlen;
    q = new double*[qlen+1];
    for(uint i=0; i<qlen+1; ++i)
      q[i] = new double[nPart];
    this->qAllocated = true;
  }

  uint nPart;           ///< Number of particles
  int step;             ///< Iteration the data correspond to

  bool ghostAllocated;  ///< Flag to tell whether the ghost arrays were allocated
  bool qAllocated;      ///< Flag to tell whether the q arrays were allocated
  
  double extension[3];  ///< Domain extension

  uint qLen;            ///< Number of q paramaters to be computed

  uint* id;             ///< Ids of the particles
  uint8_t* type;        ///< Types of the particles
  double* rx;           ///< x-positions of the particles
  double* ry;           ///< y-positions of the particles
  double* rz;           ///< z-positions of the particles
  double* vx;           ///< x-velocities of the particles
  double* vy;           ///< y-velocities of the particles
  double* vz;           ///< z-velocities of the particles
  double** q;           ///< q parameters of the particles

  uint nGhost;          ///< Number of ghost particles
  uint* ghost_id;       ///< Ids of the ghost particles
  double* ghost_rx;     ///< x-positions of the ghost particles
  double* ghost_ry;     ///< y-positions of the ghost particles
  double* ghost_rz;     ///< z-positions of the ghost particles
  std::unordered_map<int,int [26]> particles_neigh;            ///< Ids of the 26 neighbors cells for each cell
  std::unordered_map<int,int [2]> particles_start_size;        ///< Start indices and sizes of the cells in the particles attribute arrays
  std::unordered_map<int,int [2]> ghost_particles_start_size;  ///< Start indices and sizes of the cells in the ghost attribute arrays

};

#endif /* __PARTICLE_INSITU_HPP_INCLUDED */
