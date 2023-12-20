
/// @file 
/// @brief Structures for reading or writing in a legacy dump file
// ECom = to be commented or deleted by Estelle

#ifndef __STAMP_LEGACY_IO_STRUCTURES_HPP_INCLUDED
#define __STAMP_LEGACY_IO_STRUCTURES_HPP_INCLUDED


#include "utils/array/array.hpp"


/// @brief Header structure for a legacy dump file (Stamp V3 MPI-IO format)
class LegacyHeaderIOStruct {

public:

  double xmin; ///< Simulation box lower bound in the x dimension
  double xmax; ///< Simulation box upper bound in the x dimension
  double ymin; ///< Simulation box lower bound in the y dimension
  double ymax; ///< Simulation box upper bound in the y dimension
  double zmin; ///< Simulation box lower bound in the z dimension
  double zmax; ///< Simulation box upper bound in the z dimension
  double time; ///< Physical time elapsed so far
  double  CPUtime; ///< CPU time elapsed so far
  double totalEnergy; ///< Total energy of the system
  double  potentialEnergy; ///< Total potential energy of the system
  double  kineticEnergy; ///< Total kinetic energy of the system
  double  internalEnergy; ///< Total internal energy of the system (not sure this is used in Stamp)
  double  rotationalEnergy; ///< Total rotationalEnergy of the system (related to the not yet implemented quaternions)
  int particlesTotalNumber; ///< Number of particles
  int iterationNumber; ///< Step at dump time
  int domainNumber; ///< Number of domains
  double johnDoe[5]; ///< Some variables not used in ExaStamp

};


/// @brief Particle structure for a legacy dump file (Stamp V3 MPI-IO format)
class LegacyParticleIOStruct {

public:

  // LegacyParticleIOStruct()
  //   : coordinates(),
  //     johnDoe(),
  //     velocity(),
  //     quaternion(),
  //     momentumAngular(),
  //     oriantation(),
  //     particleType(0),
  //     particleID(0) {}

  double coordinates[3] = {0.,0.,0.}; ///< Particle coordinates
  double johnDoe[3] = {0.,0.,0.} ; ///< Some variables that may not even be used in Stamp
  double velocity[3] = {0.,0.,0.}; ///< Particle velocity
  double quaternion[4]= {0.,0.,0.,0.}; ///< Particle quaternion (not implemented in ExaStamp for now)
  double momentumAngular[3]= {0.,0.,0.}; ///< Particle angular momentum
  double oriantation[3] = {0.,0.,0.}; ///< Particle orientation (not implemented in ExaStamp for now)
  int    particleType; ///< Particle type
  int    particleID; ///< Particle index (Stamp indexes begin at 1 while ExaStamp indexes begin at 0)

};


/// @brief Particle structure for a legacy dump file (Stamp V3 MPI-IO format)
class LegacyDPDEParticleIOStruct {

public:

  // LegacyParticleIOStruct()
  //   : coordinates(),
  //     johnDoe(),
  //     velocity(),
  //     internalEnergy(),
  //     progress(),
  //     internalTemperature(),
  //     quaternion(),
  //     momentumAngular(),
  //     oriantation(),
  //     particleType(0),
  //     particleID(0) {}

  double coordinates[3] = {0.,0.,0.}; ///< Particle coordinates
  double johnDoe[3] = {0.,0.,0.} ; ///< Some variables that may not even be used in Stamp
  double velocity[3] = {0.,0.,0.}; ///< Particle velocity
  double internalEnergy = 0.; ///< Particle internal energy
  double internalTemperature = 0.; ///< Particle internal temperature (not used in ExaStamp)
  double progress = 0.; ///< Progress variable (for chemistry)
  double quaternion[4]= {0.,0.,0.,0.}; ///< Particle quaternion (not implemented in ExaStamp for now)
  double momentumAngular[3]= {0.,0.,0.}; ///< Particle angular momentum
  double oriantation[3] = {0.,0.,0.}; ///< Particle orientation (not implemented in ExaStamp for now)
  int    particleType; ///< Particle type
  int    particleID; ///< Particle index (Stamp indexes begin at 1 while ExaStamp indexes begin at 0)

};


/// @brief [ECom] Particles storage structure for a Hercule dump file 
///
/// Store the particles and cells
class HerculeParticleIODumpStruct {
public:

	/// @brief Default constructor
	HerculeParticleIODumpStruct(){}

	/// @brief Initialize the structure
	///
	/// Allocate the arrays of the structure with the right size
	/// @param [in] _nbCell Number of cells
	/// @param [in] _nbParticle Number of particles
	void init(int _nbCell,int _nbParticle){
		nbCell=_nbCell;
		if (nbCell){
			cells_xmin.alloc(nbCell);
			cells_xmax.alloc(nbCell);
			cells_ymin.alloc(nbCell);
			cells_ymax.alloc(nbCell);
			cells_zmin.alloc(nbCell);
			cells_zmax.alloc(nbCell);
			cells_nbParticles.alloc(nbCell);
		}
		nbParticle=_nbParticle;
		if (nbParticle){
			particles_coordinates[0].alloc(nbParticle);
			particles_coordinates[1].alloc(nbParticle);
			particles_coordinates[2].alloc(nbParticle);
			particles_johnDoe[0].alloc(nbParticle);
			particles_johnDoe[1].alloc(nbParticle);
			particles_johnDoe[2].alloc(nbParticle);
			particles_velocity[0].alloc(nbParticle);
			particles_velocity[1].alloc(nbParticle);
			particles_velocity[2].alloc(nbParticle);
			particles_quaternion[0].alloc(nbParticle);
			particles_quaternion[1].alloc(nbParticle);
			particles_quaternion[2].alloc(nbParticle);
			particles_quaternion[3].alloc(nbParticle);
			particles_momentumAngular[0].alloc(nbParticle);
			particles_momentumAngular[1].alloc(nbParticle);
			particles_momentumAngular[2].alloc(nbParticle);
			particles_oriantation[0].alloc(nbParticle);
			particles_oriantation[1].alloc(nbParticle);
			particles_oriantation[2].alloc(nbParticle);
			particles_type.alloc(nbParticle);
			particles_iD.alloc(nbParticle);
		}
	}
	/// @brief [ECom] Destructor (nothing to do)
	///
	/// Should it not desallocate the arrays ?
	virtual ~HerculeParticleIODumpStruct(){
	}
	int nbCell=0; ///< Number of cells
	Array<double> cells_xmin; ///< Lower bound in the x dimension for each cell
	Array<double> cells_xmax; ///< Upper bound in the x dimension for each cell
	Array<double> cells_ymin; ///< Lower bound in the y dimension for each cell
	Array<double> cells_ymax; ///< Upper bound in the y dimension for each cell
	Array<double> cells_zmin; ///< Lower bound in the z dimension for each cell
	Array<double> cells_zmax; ///< Upper bound in the z dimension for each cell
	Array<int> cells_nbParticles; ///< Number of particles in each cell
	
	int nbParticle=0; ///< Number of particles

  Array<double> particles_coordinates[3]; ///< Coordinates for each particle
  Array<double> particles_johnDoe[3]; ///< Some unused variable
  Array<double> particles_velocity[3]; ///< Velocity for each particle
  Array<double> particles_quaternion[4]; ///< Quaternion for each particle (not implemented in ExaStamp for now)
  Array<double> particles_momentumAngular[3]; ///< Angular momentum for each particle
  Array<double> particles_oriantation[3]; ///< Orientation for each particle (not implemented in ExaStamp for now)
  Array<int>    particles_type; ///< Type for each particle
  Array<int>    particles_iD; ///< Index for each particle
};


/// @brief [ECom] Particles storage structure for a Hercule dump file
///
/// Store the particles and cells
class HerculeDPDEParticleIODumpStruct {
public:

  /// @brief Default constructor
  HerculeDPDEParticleIODumpStruct(){}
  
  /// @brief Initialize the structure
  ///
  /// Allocate the arrays of the structure with the right size
  /// @param [in] _nbCell Number of cells
  /// @param [in] _nbParticle Number of particles
  void init(int _nbCell,int _nbParticle){
    nbCell=_nbCell;
    if (nbCell){
      cells_xmin.alloc(nbCell);
      cells_xmax.alloc(nbCell);
      cells_ymin.alloc(nbCell);
      cells_ymax.alloc(nbCell);
      cells_zmin.alloc(nbCell);
      cells_zmax.alloc(nbCell);
      cells_nbParticles.alloc(nbCell);
    }
    nbParticle=_nbParticle;
    if (nbParticle){
      particles_coordinates[0].alloc(nbParticle);
      particles_coordinates[1].alloc(nbParticle);
      particles_coordinates[2].alloc(nbParticle);
      particles_johnDoe[0].alloc(nbParticle);
      particles_johnDoe[1].alloc(nbParticle);
      particles_johnDoe[2].alloc(nbParticle);
      particles_velocity[0].alloc(nbParticle);
      particles_velocity[1].alloc(nbParticle);
      particles_velocity[2].alloc(nbParticle);
      particles_internalEnergy.alloc(nbParticle);
      particles_internalTemperature.alloc(nbParticle);
      particles_progress.alloc(nbParticle);
      particles_quaternion[0].alloc(nbParticle);
      particles_quaternion[1].alloc(nbParticle);
      particles_quaternion[2].alloc(nbParticle);
      particles_quaternion[3].alloc(nbParticle);
      particles_momentumAngular[0].alloc(nbParticle);
      particles_momentumAngular[1].alloc(nbParticle);
      particles_momentumAngular[2].alloc(nbParticle);
      particles_oriantation[0].alloc(nbParticle);
      particles_oriantation[1].alloc(nbParticle);
      particles_oriantation[2].alloc(nbParticle);
      particles_type.alloc(nbParticle);
      particles_iD.alloc(nbParticle);
    }
  }
  /// @brief [ECom] Destructor (nothing to do)
  ///
  /// Should it not desallocate the arrays ?
  virtual ~HerculeDPDEParticleIODumpStruct(){
  }
  int nbCell=0; ///< Number of cells
  Array<double> cells_xmin; ///< Lower bound in the x dimension for each cell
  Array<double> cells_xmax; ///< Upper bound in the x dimension for each cell
  Array<double> cells_ymin; ///< Lower bound in the y dimension for each cell
  Array<double> cells_ymax; ///< Upper bound in the y dimension for each cell
  Array<double> cells_zmin; ///< Lower bound in the z dimension for each cell
  Array<double> cells_zmax; ///< Upper bound in the z dimension for each cell
  Array<int> cells_nbParticles; ///< Number of particles in each cell
  
  int nbParticle=0; ///< Number of particles
  
  Array<double> particles_coordinates[3]; ///< Coordinates for each particle
  Array<double> particles_johnDoe[3]; ///< Some unused variable
  Array<double> particles_velocity[3]; ///< Velocity for each particle
  Array<double> particles_internalEnergy; ///< Internal energy for each particle
  Array<double> particles_internalTemperature; ///< Internal temperature for each particle (not used in ExaStamp)
  Array<double> particles_progress; ///< Progress variable for each particle (for chemical reactions)
  Array<double> particles_quaternion[4]; ///< Quaternion for each particle (not implemented in ExaStamp for now)
  Array<double> particles_momentumAngular[3]; ///< Angular momentum for each particle
  Array<double> particles_oriantation[3]; ///< Orientation for each particle (not implemented in ExaStamp for now)
  Array<int>    particles_type; ///< Type for each particle
  Array<int>    particles_iD; ///< Index for each particle
};


#endif // __PARTICLE_BUFFER_HPP_INCLUDED
