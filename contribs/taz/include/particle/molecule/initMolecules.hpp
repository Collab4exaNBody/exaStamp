/*
 * initMolecules.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: giarda
 */

/// @file
/// @brief Declaration of functions to initialize some molecules and forcefield associated data

#ifndef INITMOLECULES_HPP_
#define INITMOLECULES_HPP_


#include "referenceMap.hpp"

#include "domain/domainInfo.hpp"

#include "forceField/forceField.hpp"

#ifdef __use_lib_openbabel
#include "io/moleculeReader.hpp"
#endif

#include "particle/molecule/configInMol.hpp"

#include "parallel/types/MPI_atomCharged.hpp"


void configureMol(Configuration<MPI__InMol>& moleConfig, Configuration<DomainInterface>& domainConfig, Configuration<IPotential>& potentialConfig, int rank);
void initDofs(const std::vector<MPI__Particle*>& particles);


/// @brief Configure the referenceMap and forceField for a molecule input case
/// @tparam ff Forcefield
/// @tparam chM Partial charges collection method
/// @param [in,out] moleConfig Configuration of the molecule input
/// @param [in,out] domainConfig Configuration of the domains
/// @param [in,out] potentialConfig Configuration of the potentials
/// @param [in] rank Node rank
template <class PointParticle, FF::Type ff, FF::ChargeMethod chM> void configureMol(Configuration<MPI__InMol>& moleConfig, Configuration<DomainInterface>& domainConfig, Configuration<IPotential>& potentialConfig, int rank) {
#ifdef __use_lib_openbabel
	// Open the molecule reader
	MoleculeReaderOB<PointParticle,ff,chM> reader(moleConfig.m_fileName, moleConfig.m_fileDir, moleConfig.m_format, product(moleConfig.m_nCells), moleConfig.m_add, moleConfig.m_conformer, rank);
	// Get the size of the system, add the margins and multiply with the number of lattice cells
	vec3<double> minBounds, maxBounds, extension;
	reader.getSystemDimensions(minBounds, maxBounds);
	maxBounds+=moleConfig.m_margins;
	extension=(maxBounds-minBounds);
	maxBounds=minBounds+extension*moleConfig.m_nCells;
	moleConfig.m_cellSize=extension;
	// Check size with what is required from the input file and correct to fit both
	vec3<double> minAsked(domainConfig.origin), maxAsked(minAsked+domainConfig.extension);
	if(minAsked.x>minBounds.x||minAsked.y>minBounds.y||minAsked.z>minBounds.z||maxAsked.x<maxBounds.x||maxAsked.y<maxBounds.y||maxAsked.z<maxBounds.z) {
		if(domainConfig.extension.x>0&&domainConfig.extension.y>0&&domainConfig.extension.z>0) {
			std::cerr << "\nError in configureMol : ";
			std::cerr << "The extension and origin asked for the system are unsuitable to encompass both molecules and asked margins. ";
			std::cerr << "System dimensions will be set from the molecular file and margins." << std::endl;
		}
		domainConfig.origin=auxMin(minBounds,minAsked);
		maxBounds=auxMax(maxAsked,maxBounds);
		domainConfig.extension=maxBounds-domainConfig.origin;
	}
	// Get the subtypes
	moleConfig.m_listOfSubtypes=reader.getSubtypes();
	// Get the atoms and number of molecules into the molecules configuration
	reader.getAtoms(moleConfig.m_particles, moleConfig.m_numMol);
	// Initialize the force field and reference map
	Global::ffield->beginInit();
	ForceField::s_rcut=3*moleConfig.m_maxBond;
	Global::reference.configure(moleConfig, potentialConfig);
	initDofs(moleConfig.m_particles);
	Global::ffield->endInit();
#endif
}


/// @brief Configure the referenceMap and forceField for a molecule input case
///
/// Switch trough the partial charges collection method and call full-templated configureMol/// @tparam ff Forcefield
/// @param [in,out] moleConfig Configuration of the molecule input
/// @param [in,out] domainConfig Configuration of the domains
/// @param [in,out] potentialConfig Configuration of the potentials
/// @param [in] rank Node rank
template <FF::Type ff> void configureMol(Configuration<MPI__InMol>& moleConfig, Configuration<DomainInterface>& domainConfig, Configuration<IPotential>& potentialConfig, int rank) {
	switch(moleConfig.m_charges) {
	case FF::NONE :
		domainConfig.submode=Configuration<DomainInterface>::ATOM;
		configureMol<MPI__Atom, ff, FF::NONE>(moleConfig, domainConfig, potentialConfig, rank);
		break;
	case FF::FROMCOMPASS :
		configureMol<MPI__AtomCharged, ff, FF::FROMCOMPASS>(moleConfig, domainConfig, potentialConfig, rank);
		break;
	case FF::OPENBABEL :
		configureMol<MPI__AtomCharged, ff, FF::OPENBABEL>(moleConfig, domainConfig, potentialConfig, rank);
		break;
	}
}

#endif /* INITMOLECULES_HPP_ */
