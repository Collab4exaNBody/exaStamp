/// @file 
/// @brief Definition of types to exchange particles

#ifndef __MPI_PARTICLE_HPP_INCLUDED
#define __MPI_PARTICLE_HPP_INCLUDED


#include "utils/vec3/vec3.hpp"


class MPI__InMol;


/// @brief Base type to exchange a particle
class MPI__ParticleBase {

public:

  /// @brief Default constructor
  MPI__ParticleBase() : id(), ti(), r() {}

  /// @brief Destructor (nothing to do)
  virtual ~MPI__ParticleBase() {}

  /// @brief Copy constructor
  /// @param [in] particle MPI__ParticleBase to copy
  MPI__ParticleBase(const MPI__ParticleBase& particle)
    : id(particle.id), ti(particle.ti), r(particle.r) {}

  /// @brief Assignment operator
  /// @param [in] particle MPI__ParticleBase to copy
  MPI__ParticleBase& operator = (const MPI__ParticleBase& particle) {
    id = particle.id;
    ti = particle.ti;
    r  = particle.r;
    return *this;
  }

	/// @brief Accessor to the index
	inline uint64_t& index() {return id;}
	/// @brief Constant accessor to the index
	inline uint64_t index() const {return id;}
	/// @brief Accessor to the type
	inline uint8_t& type() {return ti;}
	/// @brief Constant accessor to the type
	inline uint8_t type() const {return ti;}
	/// @brief Accessor to the position
	inline vec3<double>& position() {return r;}
	/// @brief Constant accessor to the position
	inline vec3<double> position() const {return r;}

	/// @brief Shift a particle from the base cell to a replicate
	/// @param [in] idShift Shift for the index
	/// @param [in] posShift Shift for the position
	/// @param [in] notUsed For compatibility with molecules, not used, default 0
	virtual void shift(uint64_t idShift, vec3<double> posShift, uint64_t notUsed=0) {
		id+=idShift;
		r+=posShift;
	}

  /// @brief Define of the base to this particle ghost
	typedef MPI__ParticleBase Base;

  uint64_t id; ///< Index
  uint8_t  ti; ///< Type

  vec3<double> r; ///< Position

};


/// @brief Type to exchange a particle with its velocity
class MPI__Particle : public MPI__ParticleBase {

public:

  /// @brief Default constructor
  MPI__Particle() : MPI__ParticleBase(), v() {}

  /// @brief Destructor (nothing to do)
  ~MPI__Particle() {}

  /// @brief Copy constructor
	/// @param [in] particle MPI__Particle to copy
  MPI__Particle(const MPI__Particle& particle)
    : MPI__ParticleBase(particle), v(particle.v) {}

  /// @brief Assignment operator
  /// @param [in] particle MPI__Particle to copy
  MPI__Particle& operator = (const MPI__Particle& particle) {
    id = particle.id;
    ti = particle.ti;
    r  = particle.r;
    v  = particle.v;
    return *this;
  }

  /// @brief Accessor to the velocity
  inline vec3<double>& velocity() {return v;}
  /// @brief Constant accessor to the velocity
  inline vec3<double> velocity() const {return v;}

  /// @brief Define of the base to this particle ghost
  typedef MPI__ParticleBase Base;

  vec3<double> v; ///< Velocity

};

/// @brief Short name for MPI__Particle
typedef MPI__Particle Particle;


#endif // __MPI_PARTICLE_HPP_INCLUDED
