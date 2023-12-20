/*
 * initMolecules.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: giarda
 */

/// @file
/// @brief Functions to initialize some molecules and forcefield associated data


#include "globals.hpp"

#include "forceField/forceFieldUFF.hpp"
#include "forceField/forceFieldExpert.hpp"

#include "particle/molecule/initMolecules.hpp"
#include "particle/molecule/dof.hpp"

#include "parallel/types/MPI_inMol.hpp"


uint64_t MPI__InMol::s_nbPointParticles=0;


/// @brief Configure the referenceMap and forceField for a molecule input case
///
/// Switch trough the forcefields and call semi-templated configureMol
/// @param [in] domainConfig Configuration of the domains
/// @param [in] moleConfig Configuration of the molecule input
/// @param [in] rank Node rank
void configureMol(Configuration<MPI__InMol>& moleConfig, Configuration<DomainInterface>& domainConfig, Configuration<IPotential>& potentialConfig, int rank) {
  // =========================================================================
	if(rank==0) std::cout<< "  Reading molecular data";
  // =========================================================================
	switch(moleConfig.m_forceField) {
	case FF::UFF :
		Global::ffield=new ForceFieldUFF();
		configureMol<FF::UFF>(moleConfig, domainConfig, potentialConfig, rank);
		break;
	case FF::COMPASS :
//		Not implemented yet
//		Global::ffield=new ForceFieldCOMPASS;
//		configureMol<FF::COMPASS>(moleConfig, domainConfig, potentialConfig);
		break;
	case FF::AMBER :
//		Not implemented yet
//		Global::ffield=new ForceFieldAMBER;
//		configureMol<FF::AMBER>(moleConfig, domainConfig, potentialConfig);
		break;
	case FF::EXPERT :
// Not implemented yet
//		Global::ffield=new ForceFieldExpert;
		/* Debug print */ std::cerr << "Not implemented yet" << std::endl;
		break;
	}
  // =========================================================================
  if(rank==0) std::cout<< "\r  Reading molecular data ........... ok" << std::endl;
  // =========================================================================
}


/// @brief Build and register the dofs for all the atoms of the system
/// @param [in] particles Vector of the particles used to create the dofs
void initDofs(const std::vector<MPI__Particle*>& particles) {
	// Functions to get an atom from its index
	auto getSubtype = [&] (uint64_t id) -> std::string {
		MPI__InMol* ptrParticle=dynamic_cast<MPI__InMol*>(particles[id-1]);
		return Global::reference.findSubtype(ptrParticle->subtype());
	};
	auto getBondOrder = [&] (uint64_t id1, uint64_t id2) -> int {
		MPI__InMol* ptrParticle1=dynamic_cast<MPI__InMol*>(particles[id1-1]);
		return ptrParticle1->bondOrder(id2);
	};
	auto getBonds = [&] (uint64_t id) -> std::vector<uint64_t> {
		std::set<uint64_t> uniqBonds;
		MPI__InMol* ptrParticle=dynamic_cast<MPI__InMol*>(particles[id-1]);
		for(uint i(0); i<ptrParticle->nbBonds(); ++i) {
			if(ptrParticle->bonded(i)>0) uniqBonds.insert(ptrParticle->bonded(i));
		}
		std::vector<uint64_t> outBonds;
		for(uint64_t bond : uniqBonds) outBonds.push_back(bond);
		return outBonds;
	};
	// For each atom of the system
	for(uint i(0); i<particles.size(); ++i) {
		// Get the atom molecular data
		uint indexi(i+1);
		MPI__InMol* datai=dynamic_cast<MPI__InMol*>(particles[i]);
		// Loop on the bonded atoms
		for(uint j(0); j<datai->nbBonds(); ++j) {
			uint indexj=datai->bonded(j);
			MPI__InMol* dataj=dynamic_cast<MPI__InMol*>(particles[indexj-1]);
			// For each unique bond i-j (bond considered only for the lower index atom)
			if(indexj>indexi){
				// Make a bond dof
				Dof bond(indexi,indexj);
				bond.setType(Global::ffield->type(), getSubtype, getBondOrder, getBonds);
				if(!is_equal(bond.type(),"")) Global::ffield->setParam(bond.kind(),bond.type());
				// Look over neighbor of i and j to get the torsion k-i-j-l
				for(uint k(0); k<datai->nbBonds(); ++k) { for(uint l(0); l<dataj->nbBonds(); ++l) {
					if(datai->bonded(k)!=indexj&&dataj->bonded(l)!=indexi&&datai->bonded(k)!=dataj->bonded(l)) {
						// Make a tosrion dof
						Dof torsion(datai->bonded(k),indexi,indexj,dataj->bonded(l));
						torsion.setType(Global::ffield->type(), getSubtype, getBondOrder, getBonds);
						if(!is_equal(torsion.type(),"")) Global::ffield->setParam(torsion.kind(),torsion.type());
					}
				} }// End work torsion
			}// End work unique bond
			// Look each pair of neighbors of i to create the angle j-i-k
			for(uint k(j+1); k<datai->nbBonds(); ++k) {
				// Make an angle dof
				Dof angle(indexj,indexi,datai->bonded(k));
				angle.setType(Global::ffield->type(), getSubtype, getBondOrder, getBonds);
				if(!is_equal(angle.type(),"")) Global::ffield->setParam(angle.kind(),angle.type());
				// Look one more neighbor of i to create the improper torsion i-jkl
				for(uint l(k+1); l<datai->nbBonds(); ++l) {
					// Make an improper torsion dof
					Dof improper(indexi,indexj,datai->bonded(k),datai->bonded(l),true);
					improper.setType(Global::ffield->type(), getSubtype, getBondOrder, getBonds);
					if(!is_equal(improper.type(),"")) Global::ffield->setParam(improper.kind(),improper.type());
				}// End work on improper torsion
			}// End work on angle
		}// End loop on neighbors
	}// End loop on atoms
}
