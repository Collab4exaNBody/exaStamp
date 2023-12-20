/*
 * MoleculeReader.hpp
 *
 *  Created on: 10 févr. 2016
 *      Author: giarda
 */

/// @file
/// @brief Class MoleculeReader

#ifndef MOLECULEREADER_HPP_
#define MOLECULEREADER_HPP_


#include <map>
#include <vector>
#include <string>

#ifdef __use_lib_openbabel
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#endif

#include "particle/molecule/configInMol.hpp"

#include "parallel/types/MPI_atom.hpp"
#include "parallel/types/MPI_atomCharged.hpp"
#include "parallel/types/MPI_inMol.hpp"

#include "utils/stampUnits.hpp"
#include "utils/stringUtils.hpp"


class MPI__InMol;
template <class P> class MPI__PInMol;


/// @brief Tool to read molecules from a file
class MoleculeReader {

public :

	/// @brief Default constructor
	MoleculeReader() {};
	/// @brief Destructor (nothing to do)
	virtual ~MoleculeReader() {}
	/// @brief Accessor to the system boundaries
	virtual void getSystemDimensions(vec3<double>& minBounds, vec3<double>& maxBounds)=0;
	/// @brief Accessor to the list of all possible subtypes in the system
	virtual std::vector<std::string> getSubtypes()=0;
	/// @brief Get all atoms
	virtual void getAtoms(std::vector<MPI__Particle*>& atoms, uint64_t& nMol)=0;

};


#ifdef __use_lib_openbabel
/// @brief Tool to read the charge
/// @tparam PointParticle Class of the atom whose charge is read
/// @tparam chM
template <class PointParticle, FF::ChargeMethod chM> class ChargeReader {

public :

	/// @brief Default constructor
	ChargeReader() {}
	/// @brief Destructor (nothing to do)
	~ChargeReader() {}
	/// @brief Get the partial charge for an atom
	/// @param [in] atom Iterator on the atom
	/// @param [in,out] toComplete Atom to complete with the charge
	void get(MPI__PInMol<PointParticle>* toComplete, OpenBabel::OBMolAtomIter atom);

};


/// @brief Tool to read the atoms subtypes
template <FF::Type ff> class TypeReader {

public :

	/// @brief Default constructor
	TypeReader() {}
	/// @brief Destructor (nothing to do)
	~TypeReader() {}
	/// @brief Get the subtype for an atom of a molecule
	/// @param [in] atom Iterator on the atom
	/// @return Subtype
	std::string get(OpenBabel::OBMolAtomIter atom);
};


/// @brief Shortcut for a template that depend on the force field and charge method
#define TMPLMR template<class PointParticle, FF::Type ff, FF::ChargeMethod chM>


/// @brief Tool to read molecules from a file with open babel
TMPLMR class MoleculeReaderOB : MoleculeReader {

public :

	MoleculeReaderOB(std::string fileName="", std::string fileDir=".", std::string format="pdb", int multiplier=1, bool add=false, int conformer=1, int rank=0);
	/// @brief Destructor (nothing to do)
	~MoleculeReaderOB() {}
	void getSystemDimensions(vec3<double>& minBounds, vec3<double>& maxBounds);
	std::vector<std::string> getSubtypes();
	void getAtoms(std::vector<MPI__Particle*>& atoms, uint64_t& nMol);

private :

	void readFile(int multiplier, int conformer, int rank);
	void build();
	/// @brief Get the global index for an atom
	/// @param [in] indexInMol Index of the atom in its molecule
	/// @param [in] iMol Index of the molecule
	inline uint64_t index(uint64_t indexInMol, uint64_t iMol) {return m_indexes[std::make_pair(iMol,indexInMol)];}
	void makeBondOrders(MPI__InMol* toComplete, OpenBabel::OBMolAtomIter atom, uint64_t numMol);
	// Parameters
	std::string m_fileName; ///< Name of the file to read
	std::string m_fileDir; ///< Directory of the molecular input file
	std::string m_format; ///< Format of the molecular input file
	bool m_addH; ///< Flag to indicate that hydrogens must be added to the molecule
	// Tools
	ChargeReader<PointParticle, chM> m_chargeReader; ///< Templated tool to read the charges
	TypeReader<ff> m_typeReader; ///< Templated tool to read the subtypes
	// Build particles and collected data
	std::vector<OpenBabel::OBMol> m_molecules; ///< Molecules read in the molecular input file
	std::map<std::pair<uint64_t,uint64_t>,uint64_t> m_indexes; ///< Global indexes for the atoms of each molecule
	std::vector<MPI__PInMol<PointParticle> > m_atoms; ///< Atoms in construction
	vec3<double> m_min; ///< Min bounds for atoms in the system
	vec3<double> m_max; ///< Max bounds for atoms in the system
	uint64_t m_nMol; ///< Number of true molecules
	// Initialize the subtypes
	std::vector<std::string> m_listOfSubtypes; ///< List of all possible subtypes in the system

};


#include "moleculeReader.hxx"


#undef TMPLMR

#endif
#endif /* MOLECULEREADER_HPP_ */
