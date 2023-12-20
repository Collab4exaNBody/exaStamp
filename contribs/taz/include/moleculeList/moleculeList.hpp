/*
 * moleculeList.hpp
 *
 *  Created on: Sep 15, 2016
 *      Author: giarda
 */

/// @file
/// @brief Define class MoleculeList, an interface to handle multiple molecules


#ifndef MOLECULELIST_HPP_
#define MOLECULELIST_HPP_


#include <unordered_map>

#include "grid/correcter.hpp"

#include "moleculeList/vectorizationMBuffer.hpp"

#include "parallel/thread/thread.hpp"

#include "particle/molecule/molecule.hpp"


/// @brief Base chunk for the vectorization on the molecules
#define MOL_CHUNK 128
/// @brief Chunk adjustement for the vectorization on the molcules
#define MOL_ADJUST 8


/// @brief Interface to access the molecules
class MoleculeInterface {

public :

	/// @brief Default constructor
	MoleculeInterface() {}
	/// @brief Destructor (nothing to do)
	virtual ~MoleculeInterface() {}
	/// @brief Interface to update the positions of the molecules
	virtual void updatePositions() {}
	/// @brief Interface to update the velocities of the molecules
	virtual void updateVelocities() {}
	/// @brief Interface to check that total molecular forces are null
	virtual void checkForcesConservation() {}
	/// @brief Interface to update the neighbors of the point particles according to their bonds
	virtual void updateNeighbors() {}
	/// @brief Interface to update the neighbors of the point particles according to their bonds
	/// for the particles in the inside cells
	virtual void updateNeighborsInside() {}
	/// @brief Interface to update the neighbors of the point particles according to their bonds
	/// for the particles in the edge cells
	virtual void updateNeighborsEdge() {}
	/// @brief Interface to compute the intra-molecular forces
	virtual void computeForces() {}
	/// @brief Interface to check the size of local arrays
	virtual void checkLocals()	{}
	/// @brief Interface to set the distance corrector
	virtual void setCorrector(Correcter* corrector) {}

protected :

};


/// @brief Tool to handle the list of molecules
/// @tparam cList Class of the associated cellLists
template <class CList, size_t align> class MoleculeList : public MoleculeInterface {

public :

	/// @brief Shortcut for the molecule type
	typedef Molecule<CList> mol_t;

	/// @brief Constructor
	/// @param [in] cells Pointer the associated cellLists
	MoleculeList(CList* cells) : m_cells(cells), m_mtx(), m_mapMtx() {
		m_mtx.check(1);
	}
	/// @brief Destructor (nothing to do)
	virtual ~MoleculeList() {
		m_cells=nullptr;
		m_molecules.clear();
	}
	virtual void updatePositions();
	virtual void updateVelocities();
	virtual void checkForcesConservation();
	virtual void updateNeighbors();
	virtual void updateNeighborsInside();
	virtual void updateNeighborsEdge();
	virtual void computeForces();
	virtual void checkLocals();
	virtual void setCorrector(Correcter* corrector);
	void add(uint64_t iMol, uint64_t iPP, uint iCell, uint lId);
	void remove(uint64_t iMol, uint64_t iPP);
	void shift(uint newPos, uint64_t iMol, uint64_t iPP);
	uint countNeighbors(uint64_t iPP, const uint level, uint64_t iMol);

protected :

	template <typename Lambda> void for_all_molecules(Lambda lambda);
	void makeChunksAndStacks(ExtArray<std::tuple<uint64_t, uint, uint> >& chunks, ExtArray<uint>& stacks);
	template <uint8_t kind> VectMBuffer<align>& getVectMBuffer(const uint& n);

	static Th_<VectMBuffer<align> > s_shared_vectMBuffer; ///< Local memory-aligned arrays to perform vectorized operations

	CList* m_cells; ///< Pointer to the associated cellLists
	std::unordered_map<uint64_t, mol_t > m_molecules; ///< Map of the molecules and their indexes
	MMutex m_mtx; ///< Array of mutexes
	MrwMutex m_mapMtx; ///< Read-write mutex used to lock the map during molecule add
	Correcter* m_corrector; ///< Tool to correct the positions in periodic boundaries case
};

/// @brief Initialization of m_shared_vectBuffer
template <class CList, size_t align> Th_<VectMBuffer<align> > MoleculeList<CList, align>::s_shared_vectMBuffer = Th_ <VectMBuffer<align> >();


/// @brief Update the position of each molecule
template <class CList, size_t align> void MoleculeList<CList, align>::updatePositions() {
	for_all_molecules([&](mol_t& molecule){
		molecule.updatePosition();
	});
}


/// @brief Update the velocity of each molecule
template <class CList, size_t align> void MoleculeList<CList, align>::updateVelocities() {
	for_all_molecules([&](mol_t& molecule){
		molecule.updateVelocity();
	});
}


/// @brief Calculate and print the total force on each molecule
template <class CList, size_t align> void MoleculeList<CList, align>::checkForcesConservation() {
	for_all_molecules([&](mol_t& molecule){
		molecule.checkForcesConservation();
	});
}


/// @brief For each molecule, update the neighbors of the point particles according to their bonds
///
///
template <class CList, size_t align> void MoleculeList<CList, align>::updateNeighbors() {
	for_all_molecules([&](mol_t& molecule){
		molecule.template updateNbr<0>();
	});
}


/// @brief For each molecule, update the neighbors of the point particles in the inside cells according to their bonds
///
///
template <class CList, size_t align> void MoleculeList<CList, align>::updateNeighborsInside() {
	for_all_molecules([&](mol_t& molecule){
		molecule.template updateNbr<1>();
	});
}


/// @brief For each molecule, update the neighbors of the point particles in the edge cells according to their bonds
///
///
template <class CList, size_t align> void MoleculeList<CList, align>::updateNeighborsEdge() {
	for_all_molecules([&](mol_t& molecule){
		molecule.template updateNbr<2>();
	});
}


/// @brief Compute the intra-molecular forces
template <class CList, size_t align> void MoleculeList<CList, align>::computeForces() {
	// I make chunks of molecules and give stacks of chunks to the threads
	// A chunk contain :
	//  - a molecule id
	//  - the number of atoms to take in that molecule
	//  - the position of the first atom to take
	// The stacks container contain the first chunk of each stack
	ExtArray<std::tuple<uint64_t,uint,uint> > atomChunks;
	ExtArray<uint> atomStacks;
	atomStacks.push_back(0);
	makeChunksAndStacks(atomChunks, atomStacks);
	// Big computation loop
  parallel_region(0, atomStacks.size()-1, [&](const uint begin, const uint end) {
  	// For each well sized stack
  	for(uint iStack(begin); iStack<end; ++iStack) {
  		// Local variables :
  		// 	Size and position of the chunks
  		uint size, posInMol, posInBuffer(0);
  		//	Molecule ID for each chunk
  		uint64_t iMol;
  		//	Dofs for all the chunks
  		std::vector<Dof> bonds, angles, dihedrals, impropers;
  		//  Number of each dof in each chunk
  		std::vector<uint> nBChunk(atomStacks[iStack+1]-atomStacks[iStack]+1,0), nAChunk(atomStacks[iStack+1]-atomStacks[iStack]+1,0),
  				nDChunk(atomStacks[iStack+1]-atomStacks[iStack]+1,0), nIChunk(atomStacks[iStack+1]-atomStacks[iStack]+1,0);
  		// Build dofs
  		//	Loop on the chunks
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			size=std::get<1>(atomChunks[iChunk+atomStacks[iStack]]);
  			posInMol=std::get<2>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].buildDofs(posInMol, size, bonds, angles, dihedrals, impropers);
  			nBChunk[iChunk+1]=bonds.size();
  			nAChunk[iChunk+1]=angles.size();
  			nDChunk[iChunk+1]=dihedrals.size();
  			nIChunk[iChunk+1]=impropers.size();
  		}
			// Compute forces on bonds
    	VectMBuffer<align>& buffer=getVectMBuffer<0>(bonds.size());
  		//	Loop on the chunks to fill the buffer
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			size=nBChunk[iChunk+1]-nBChunk[iChunk];
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].fillVectMBufferBond(buffer,size, posInBuffer, nBChunk[iChunk], bonds);
  			posInBuffer+=size;
  		}
  		m_corrector->correct(posInBuffer, buffer.rx1(), buffer.ry1(), buffer.rz1());
  		// 	Compute forces
  		Global::ffield->computeForceBonds(posInBuffer, buffer.ep(), buffer.fx1(), buffer.fy1(), buffer.fz1(),
  				buffer.rx1(), buffer.ry1(), buffer.rz1(), buffer.params(std::get<0>(ForceField::paramBuild.at(Global::ffield->type()))), FF::BHARM);
  		// 	Loop on the chunks to store the forces
  		posInBuffer=0;
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			size=nBChunk[iChunk+1]-nBChunk[iChunk];
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].writeBondForces(buffer, size, posInBuffer, nBChunk[iChunk], bonds);
  			m_molecules[iMol].writeBondViriel(buffer, size, posInBuffer, nBChunk[iChunk], bonds);
  			posInBuffer+=size;
  		}
  		bonds.clear();
			// Compute forces on angles
  		buffer.template resize<1>(angles.size());
  		posInBuffer=0;
  		//	Loop on the chunks to fill the buffer
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			size=nAChunk[iChunk+1]-nAChunk[iChunk];
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].fillVectMBufferAngle(buffer,size, posInBuffer, nAChunk[iChunk], angles);
  			posInBuffer+=size;
  		}
  		m_corrector->correct(posInBuffer, buffer.rx1(), buffer.ry1(), buffer.rz1());
  		m_corrector->correct(posInBuffer, buffer.rx2(), buffer.ry2(), buffer.rz2());
  		// 	Compute forces
  		Global::ffield->computeForceAngles(posInBuffer, buffer.ep(), buffer.fx1(), buffer.fy1(), buffer.fz1(), buffer.fx2(), buffer.fy2(), buffer.fz2(),
  				buffer.rx1(), buffer.ry1(), buffer.rz1(), buffer.rx2(), buffer.ry2(), buffer.rz2(), buffer.params(std::get<1>(ForceField::paramBuild.at(Global::ffield->type()))), FF::ACOSTO3);
  		// 	Loop on the chunks to store the forces
  		posInBuffer=0;
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			size=nAChunk[iChunk+1]-nAChunk[iChunk];
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].writeAngleForces(buffer, size, posInBuffer, nAChunk[iChunk], angles);
  			m_molecules[iMol].writeAngleViriel(buffer, size, posInBuffer, nAChunk[iChunk], angles);
  			posInBuffer+=size;
  		}
  		angles.clear();
			// Compute forces on dihedral angles
  		buffer.template resize<2>(dihedrals.size());
  		posInBuffer=0;
  		//	Loop on the chunks to fill the buffer
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			size=nDChunk[iChunk+1]-nDChunk[iChunk];
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].fillVectMBufferDihedral(buffer,size, posInBuffer, nDChunk[iChunk], dihedrals);
  			posInBuffer+=size;
  		}
  		m_corrector->correct(posInBuffer, buffer.rx1(), buffer.ry1(), buffer.rz1());
  		m_corrector->correct(posInBuffer, buffer.rx2(), buffer.ry2(), buffer.rz2());
  		m_corrector->correct(posInBuffer, buffer.rx3(), buffer.ry3(), buffer.rz3());
  		// 	Compute forces
  		Global::ffield->computeForceDihedrals(posInBuffer, buffer.ep(), buffer.fx1(), buffer.fy1(), buffer.fz1(), buffer.fx2(), buffer.fy2(), buffer.fz2(), buffer.fx3(), buffer.fy3(), buffer.fz3(),
  				buffer.rx1(), buffer.ry1(), buffer.rz1(), buffer.rx2(), buffer.ry2(), buffer.rz2(), buffer.rx3(), buffer.ry3(), buffer.rz3(),
  				buffer.params(std::get<2>(ForceField::paramBuild.at(Global::ffield->type()))), FF::DCOSN);
  		// 	Loop on the chunks to store the forces
  		posInBuffer=0;
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			size=nDChunk[iChunk+1]-nDChunk[iChunk];
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].writeDihedralForces(buffer, size, posInBuffer, nDChunk[iChunk], dihedrals);
  			m_molecules[iMol].writeDihedralViriel(buffer, size, posInBuffer, nDChunk[iChunk], dihedrals);
  			posInBuffer+=size;
  		}
  		dihedrals.clear();
			// Compute forces on improper torsions
  		buffer.template resize<3>(impropers.size());
  		posInBuffer=0;
  		//	Loop on the chunks to fill the buffer
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			size=nIChunk[iChunk+1]-nIChunk[iChunk];
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].fillVectMBufferImproper(buffer,size, posInBuffer, nIChunk[iChunk], impropers);
  			posInBuffer+=size;
  		}
  		m_corrector->correct(posInBuffer, buffer.rx1(), buffer.ry1(), buffer.rz1());
  		m_corrector->correct(posInBuffer, buffer.rx2(), buffer.ry2(), buffer.rz2());
  		m_corrector->correct(posInBuffer, buffer.rx3(), buffer.ry3(), buffer.rz3());
  		// 	Compute forces
  		Global::ffield->computeForceImpropers(posInBuffer, buffer.ep(), buffer.fx1(), buffer.fy1(), buffer.fz1(), buffer.fx2(), buffer.fy2(), buffer.fz2(), buffer.fx3(), buffer.fy3(), buffer.fz3(),
  				buffer.rx1(), buffer.ry1(), buffer.rz1(), buffer.rx2(), buffer.ry2(), buffer.rz2(), buffer.rx3(), buffer.ry3(), buffer.rz3(),
  				buffer.params(std::get<3>(ForceField::paramBuild.at(Global::ffield->type()))), FF::ICOSTO2);
  		// 	Loop on the chunks to store the forces
  		posInBuffer=0;
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			size=nIChunk[iChunk+1]-nIChunk[iChunk];
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_molecules[iMol].writeImproperForces(buffer, size, posInBuffer, nIChunk[iChunk], impropers);
  			m_molecules[iMol].writeImproperViriel(buffer, size, posInBuffer, nIChunk[iChunk], impropers);
  			posInBuffer+=size;
  		}
  		impropers.clear();
  		// Clear finished molecules
  		for(uint iChunk(0); iChunk<atomStacks[iStack+1]-atomStacks[iStack]; ++iChunk) {
  			iMol=std::get<0>(atomChunks[iChunk+atomStacks[iStack]]);
  			size=std::get<1>(atomChunks[iChunk+atomStacks[iStack]]);
  			m_mtx.lock(iMol%m_mtx.size());
  			m_molecules[iMol].finishedWithDofs(size);
  			m_mtx.unlock(iMol%m_mtx.size());
  		}
  	}
  });
}


/// @brief Add a point particle to a molecule
/// @tparam CList Class of the associated cellLists
/// @param [in] iMol Index of the molecule
/// @param [in] iPP Index of the point particle
/// @param [in] iCell Cell of the point particle
/// @param [in] lId Position of the point particle in the cell
template <class CList, size_t align> void MoleculeList<CList, align>::add(uint64_t iMol, uint64_t iPP, uint iCell, uint lId) {
	if(iMol!=0) {
		m_mtx.lock(iMol%m_mtx.size());
		m_mapMtx.lock_read();
			if(m_molecules.find(iMol)==m_molecules.end()) {
				Molecule<CList> mol(iMol, m_cells);
				m_mapMtx.switch_to_write();
				m_molecules.insert(std::make_pair(iMol,mol));
				m_mapMtx.switch_to_read();
			}
			m_molecules[iMol].addAtom(iPP, iCell, lId);
		m_mapMtx.unlock();
		m_mtx.unlock(iMol%m_mtx.size());
	}
}


/// @brief Check and resize the array of mutexes
template <class CList, size_t align> void MoleculeList<CList, align>::checkLocals() {
	uint64_t maxMol(0), minMol(0);
  for(auto it=m_molecules.begin(); it!=m_molecules.end(); ++it) {
		if(minMol==0) minMol=it->first;
		else minMol=auxMin(minMol, it->first);
		maxMol=auxMax(maxMol, it->first);
	}
	m_mtx.check(uint(maxMol-minMol)+1);
}


/// @brief Set the corrector
/// @param [in] corrector The grid corrector
template <class CList, size_t align> void MoleculeList<CList, align>::setCorrector(Correcter* corrector) {
	m_corrector=corrector;
}


/// @brief Remove a point particle from a molecule
/// @tparam CList Class of the associated cellLists
/// @param [in] iMol Index of the molecule
/// @param [in] iPP Index of the point particle
template <class CList, size_t align> void MoleculeList<CList, align>::remove(uint64_t iMol, uint64_t iPP) {
	if(iMol!=0) {
		m_mtx.lock(iMol%m_mtx.size());
		m_mapMtx.lock_read();
		if(m_molecules.find(iMol)!=m_molecules.end()) {
			if(m_molecules[iMol].getSize()!=0) m_molecules[iMol].removeAtom(iPP);
			else {
				m_mapMtx.switch_to_write();
				m_molecules.erase(iMol);
				m_mapMtx.switch_to_read();
			}
		}
		m_mapMtx.unlock();
		m_mtx.unlock(iMol%m_mtx.size());
	}
}


/// @brief Correct the position of an atom in its cell
/// @tparam CList Class of the associated cellLists
/// @param [in] newPos Shift value for the last atom of the cell
/// @param [in] iMol Index of the molecule
/// @param [in] iPP Index of the point particle
template <class CList, size_t align> void MoleculeList<CList, align>::shift(uint newPos, uint64_t iMol, uint64_t iPP) {
	if(iMol!=0) m_molecules[iMol].shift(iPP,newPos);
}


/// @brief Get the number of atoms linked to the specified atom by the specified number of bonds
/// @param [in] iPP Index of the specified atom
/// @param [in] level Number of linking bonds
/// @param [in] iMol Index of the particle of the specified bond
/// @return Number of linked atoms
template <class CList, size_t align> uint MoleculeList<CList, align>::countNeighbors(uint64_t iPP, const uint level, uint64_t iMol) {
	return m_molecules[iMol].getNumberOfNeighbors(iPP, level);
}


/// @brief Apply a function to all the molecules
/// @tparam Lambda Class of the function
/// @param [in] lambda Function to apply
template <class CList, size_t align> template <typename Lambda> void MoleculeList<CList, align>::for_all_molecules(Lambda lambda) {
	// This vector creation enable the use of thread of the molecule map
	// Feel free to change this part if you have a better solution
	std::vector<uint64_t> molIds;
  for(auto it=m_molecules.begin(); it!=m_molecules.end(); ++it) molIds.push_back(it->first);
	// Parallel work
  parallel_region(0, molIds.size(), [&](const uint begin, const uint end) {
  	for(uint i=begin; i<end; ++i) {
  		lambda(m_molecules[molIds[i]]);}
    });
}


/// @brief Build chunks of molecules and give stacks of chunks to the threads
///
/// A chunk contain :
///  - a molecule id
///  - the number of atoms to take in that molecule
///  - the position of the first atom to take
/// The stacks container contain the first chunk of each stack
/// @param [out] chunks Chunks
/// @param [out] stacks Stacks
template <class CList, size_t align> void MoleculeList<CList, align>::makeChunksAndStacks(ExtArray<std::tuple<uint64_t, uint, uint> >& chunks, ExtArray<uint>& stacks) {
	uint64_t idMol;
	uint chunkSize, stackSize, position, lasting;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	stackSize=0;
  for(auto it=m_molecules.begin(); it!=m_molecules.end(); ++it) {
  	// Initialize the building of dofs in this molecule
  	it->second.initBuilding();
		// Get the id
		idMol=it->first;
		// Get the number of atoms
		lasting=it->second.getRealSize();
		// Begin from the first a
		position=0;
		while(lasting!=0) { // While the are unadded atoms in this molecule
			if(stackSize+lasting<MOL_CHUNK+MOL_ADJUST) { // If there is enough room, add all the lasting atoms
				stackSize+=lasting; // Increase the number of atoms in the stack
				chunks.push_back(std::make_tuple(idMol, lasting, position)); // Add the chunk
				lasting=0; // Decrease the number of lasting atoms accordingly
				if(stackSize>=MOL_CHUNK) { // If the stack is full
					stacks.push_back(chunks.size()); // Add to the list
					stackSize=0; // Reset the stack size
				}
			}
			else { // If there isn't enough room
				if(stackSize<=MOL_CHUNK-MOL_ADJUST) {	// And what room is left is not ridiculously small, add a chunk of the molecule
					chunkSize=MOL_CHUNK-stackSize; // Add to fill the stack
					chunks.push_back(std::make_tuple(idMol, chunkSize, position)); // Add the chunk
					lasting-=chunkSize; // Decrease the number of lasting atoms accordingly
					position+=chunkSize; // Adjust the position accordingly
					stackSize+=chunkSize;
				}
				// Then make a new stack
				stacks.push_back(chunks.size()); // Add to the list
				stackSize=0; // Reset the stack size
			}
		} // No more of those dofs in the molecule
	} // End of the loop on the molecules
	if(stackSize!=0) stacks.push_back(chunks.size()); // Add to the list
}


/// @brief Get a local version of the vectorization buffer and resize it
/// @param [in] n Size
/// @return Local vectorization buffer
template <class CList, size_t align> template <uint8_t kind> VectMBuffer<align>& MoleculeList<CList, align>::getVectMBuffer(const uint& n)  {
  // Get the local version
	auto& vb = s_shared_vectMBuffer.local();
  // Resize to n
	vb.template resize<kind>(n);
  return vb;
}

#endif /* MOLECULELIST_HPP_ */
