/*
 * MPI_atomCharged.hpp
 *
 *  Created on: May 30, 2016
 *      Author: giarda
 */

/// @file
/// @brief Class MPI__AtomCharged

#ifndef MPI_ATOMCHARGED_HPP_
#define MPI_ATOMCHARGED_HPP_


#include "parallel/types/MPI_atom.hpp"


/// @brief Base class for atoms with a charge
///
/// Used for calculations with coulomb interactions
class MPI__AtomChargedBase : public MPI__ParticleBase {

public:

	/// @brief Default constructor
	MPI__AtomChargedBase() :
		MPI__ParticleBase(),
		m_charge(0.)
	{}
	/// @brief Copy constructor
	/// @param [in] other Other atom
	MPI__AtomChargedBase(const MPI__AtomChargedBase &other) :
		MPI__ParticleBase(other),
		m_charge(other.m_charge)
	{}

	/// @brief Assignement operator
	/// @param [in] other Atom to copy
	MPI__AtomChargedBase& operator = (const MPI__AtomChargedBase &other) {
		id = other.id;
		ti = other.ti;
		r  = other.r;
		m_charge=other.m_charge;
		return *this;
	}
	/// @brief Destructor (nothing to do)
	~MPI__AtomChargedBase() {}

	/// @brief Accessor to the charge
	double& charge() {return m_charge;}
	/// @brief Constant accessor to the charge
	double charge() const {return m_charge;}

  /// @brief Define of the base to this particle ghost
	typedef MPI__AtomChargedBase Base;

protected :

	double m_charge; ///< Charge

};


/// @brief Class for atoms with a charge
///
/// Used for calculations with coulomb interactions
class MPI__AtomCharged : public MPI__Atom {

public:

	/// @brief Default constructor
	MPI__AtomCharged() :
		MPI__Particle(),
		m_charge(0.)
	{}
	/// @brief Copy constructor
	/// @param [in] other Other atom
	MPI__AtomCharged(const MPI__AtomCharged &other) :
		MPI__Particle(other),
		m_charge(other.m_charge)
	{}

	/// @brief Assignement operator
	/// @param [in] other Atom to copy
	MPI__AtomCharged& operator = (const MPI__AtomCharged &other) {
		id = other.id;
		ti = other.ti;
		r  = other.r;
		v  = other.v;
		m_charge=other.m_charge;
		return *this;
	}
	/// @brief Destructor (nothing to do)
	~MPI__AtomCharged() {}

	/// @brief Accessor to the charge
	double& charge() {return m_charge;}
	/// @brief Constant accessor to the charge
	double charge() const {return m_charge;}

  /// @brief Define of the base to this particle ghost
	typedef MPI__AtomChargedBase Base;

protected :

	double m_charge; ///< Charge

};

#endif /* MPI_ATOMCHARGED_HPP_ */
