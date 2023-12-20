/* MPI_inMol.hpp
 *
 * Created on: May 30, 2016
 * from MPI_atomInMol.hpp, created on: Feb 22, 2016
 *
 *      Author: giarda
 */

/// @file
/// @brief Classes MPI__InMol and MPI__PInMol

#ifndef MPI_INMOL_HPP_
#define MPI_INMOL_HPP_


#include <set>
#include <iostream>

#include "utils/vec3/vec3.hpp"


/// @brief Maximum number of bonds for a molecule
#define MAX_NUM_BONDS 6


/// @brief Class to transfer the part to add to a point particle to use it in a molecule
/// (bonds, molecule and subtype
class MPI__InMol {

public :
	/// @brief Default constructor
	///
	///
	MPI__InMol() :
		m_molecule(0),
		m_subtype(0)
		{
			for(uint i=0;i<MAX_NUM_BONDS;i++) m_bonded[i]=0;
		}
/// @brief Copy constructor
/// @param [in] other Other atom
	MPI__InMol(const MPI__InMol &other) :
		m_molecule(other.m_molecule),
		m_subtype(other.m_subtype)
	{
		for(uint i=0;i<MAX_NUM_BONDS;i++) m_bonded[i]=other.m_bonded[i];
	}
	/// @brief Assignement operator
	/// @param [in] other Atom to copy
  MPI__InMol& operator = (const MPI__InMol &other) {
    m_molecule = other.m_molecule;
    m_subtype=other.m_subtype;
		for(uint i=0;i<MAX_NUM_BONDS;i++) m_bonded[i]=other.m_bonded[i];
    return *this;
  }
	/// @brief Destructor (nothing to do)
	virtual ~MPI__InMol() {}
	/// @brief Accessor to the owner molecule index
	inline uint64_t& molecule() {return m_molecule;}
	/// @brief Constant accessor to the owner molecule index
	inline uint64_t molecule() const {return m_molecule;}
	/// @brief Constant accessor to the indexes of the atoms bonded to this one
	/// @param [in] i Neighbor to access
	inline uint64_t bonded(uint i) const {return m_bonded[i];}
	/// @brief Accessor to the subtype
	inline uint8_t& subtype() {return m_subtype;}
	/// @brief Constant accessor to the subtype
	inline uint8_t subtype() const {return m_subtype;}
	/// @brief Get the number of bonds for this molecule
	/// @return Number of bonds
	inline uint nbBonds() const {
		std::set<uint64_t> bonds;
		for(uint i(0); i<MAX_NUM_BONDS;++i) if(m_bonded[i]>0&&m_bonded[i]<s_nbPointParticles+1) bonds.insert(m_bonded[i]);
		return bonds.size();
	}
	/// @brief Reset the bonds
	inline void resetBonds() {
		for(uint i(0); i<MAX_NUM_BONDS;++i) m_bonded[i]=0;
	}
	/// @brief Get a the next empty place in the bonded array
	/// @return Reference to the next empty place of the last place if the array is full
	inline uint64_t& nextBond() {
		for(uint i(0); i<MAX_NUM_BONDS;++i) if(m_bonded[i]==0) return m_bonded[i];
		std::cerr << "Error in MPI__InMol::nextBond : more than " << MAX_NUM_BONDS << " bonds ";
		std::cerr << "! Previous last bond will be deleted." << std::endl;
		return m_bonded[MAX_NUM_BONDS-1];
	}
	/// @brief Get the order of the bond between this atom and a neighbor
	/// @param [in] neighbor Index of the neighbor
	/// @return Bond order or -1 if partial bond order (aromatic or amide)
	inline int bondOrder(const uint64_t neighbor) const {
		int bondOrder(0);
		for(uint64_t atomId :  m_bonded) {
			if(atomId==neighbor+s_nbPointParticles) return -1;
			else if(atomId==neighbor) bondOrder++;
		}
		return bondOrder;
	}

	static uint64_t s_nbPointParticles; ///< Number of point particles in the system

protected :

	uint64_t m_molecule; ///< Owner molecule index
	uint8_t m_subtype; ///< Integer reference to the subtype of the atom
	// There shouldn't be a case where the is more than 256 atom types, if such case arise,
	// consider upgrading m_subtype to an uint16_t
	uint64_t m_bonded[MAX_NUM_BONDS]; ///< Indexes of the atoms bonded to this one

};


/// @brief Class to transfer a PointParticle that's inside a molecule
/// @tparam MPI__PointParticle The point particle the new particle inherit from
template<class MPI__PointParticle> class MPI__PInMol : public MPI__PointParticle, public MPI__InMol {

public :

	/// @brief Shortcut for the type of the point particle
	typedef MPI__PointParticle PP;

	/// @brief Default constructor
	MPI__PInMol() :
		MPI__PointParticle(),
		MPI__InMol()
	{}

	/// @brief Constructor
	/// @param [in] other Particle to copy
	MPI__PInMol(const MPI__PInMol<MPI__PointParticle> &other) :
		MPI__PointParticle(other),
		MPI__InMol(other)
	{}

  /// @brief Assignment operator
  /// @param [in] other MPI__ParticleBase to copy
  MPI__PInMol<MPI__PointParticle>& operator = (const MPI__PInMol<MPI__PointParticle>& other) {
  	this->MPI__PointParticle::operator=(other);
  	this->MPI__InMol::operator=(other);
    return *this;
  }

	/// @brief Destructor (nothing to do)
	virtual ~MPI__PInMol() {}

	/// @brief Shift a particle from the base cell to a replicate
	/// @param [in] idShift Shift for the index
	/// @param [in] posShift Shift for the position
	/// @param [in] molShift Shift for the molecule index
	virtual void shift(uint64_t idShift, vec3<double> posShift, uint64_t molShift) {
		MPI__PointParticle::shift(idShift,posShift);
		for(uint i(0); i<MAX_NUM_BONDS; ++i) if(m_bonded[i]!=0) m_bonded[i]+=idShift;
		if(m_molecule!=0) m_molecule+=molShift;
	}

  /// @brief Define of the base to this particle ghost
  typedef MPI__PInMol<typename MPI__PointParticle::Base> Base;

};

#endif /* MPI_INMOL_HPP_ */
