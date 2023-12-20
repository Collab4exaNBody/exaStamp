/*
 * molecule.hpp
 *
 *  Created on: Feb 5, 2016
 *      Author: giarda
 */

/// @file
/// @brief Definition of the molecule

#ifndef MOLECULE_HPP_
#define MOLECULE_HPP_


#include <unordered_map>
#include <vector>
#include <set>

#include "moleculeList/vectorizationMBuffer.hpp"

#include "parallel/types/MPI_particle.hpp"

#include "particle/molecule/dof.hpp"


/// @brief Type for a molecule
#define TYPEMOL 0


/// @brief Class to store the position of an point particle in the grid
class Position {

public :

	/// @brief Default constructor
	Position() :
		m_cell(0),
		m_index(0)
	{}
	/// @brief Constructor
	/// @param [in] other Position to copy
	Position(const Position& other) :
		m_cell(other.m_cell),
		m_index(other.m_index)
	{}
	/// @brief Constructor
	/// @param [in] cell Cell
	/// @param [in] index Position in the cell
	Position(const uint cell, const uint index) :
		m_cell(cell),
		m_index(index)
	{}

	/// @brief Destructor (nothing to do)
	~Position() {}

	/// @brief Reference accessor to the cell
	uint& cell() {return m_cell;}
	/// @brief Constant accsessor to the cell
	uint cell() const {return m_cell;}
	/// @brief Reference accessor to the local index
	uint& index() {return m_index;}
	/// @brief Constant accsessor to the local index
	uint index() const {return m_index;}

private :

	uint m_cell; ///< Cell
	uint m_index; ///< Position of the point particule in the cell

};


/// @brief Shortcut for a template that depend on a point particle, a chunk and an alignment
#define TMPLM template <class CList>
/// @brief Shortcut for a cellList with molecules
#define TMPL_Molecule Molecule<CList>


/// @brief Class to handle a molecule
TMPLM class Molecule : public Particle {

public :

	/// @brief Default constructor
	Molecule() :
	Particle(),
	m_pointParticles(),
	m_allIDs(),
	m_idUsed(0),
	m_cells(nullptr)
	{}
	/// @brief Copy constructor
	/// @param [in] other Molecule to copy
	Molecule(const Molecule<CList>& other) :
	Particle(other),
	m_pointParticles(other.m_pointParticles),
	m_allIDs(),
	m_idUsed(0),
	m_cells(other.m_cells)
	{}
	Molecule(const uint64_t index, CList* cellList);
	virtual ~Molecule();
	void initBuilding();
	void finishedWithDofs(const uint size);
	void buildDofs(const uint begin, const uint size,
			std::vector<Dof>& bonds, std::vector<Dof>& angles, std::vector<Dof>& dihedrals, std::vector<Dof>& impropers);
	template <uint level> void updateNbr();
	void getNeighbors(const uint64_t start, std::set<uint64_t>& neighbors, const uint level);
	uint getNumberOfNeighbors(const uint64_t start, const uint level);
	/// @brief Get the number of point particles in the molecule
	uint getSize() {return m_pointParticles.size();}
	/// @brief Get the number of real point particles in the molecule
	///
	/// Use only after initBuilding()
	uint getRealSize() {return m_allIDs.size();}
	/// @brief Add a point particle to the molecule
	/// @param [in] gIndex Global index of the point particle
	/// @param [in] cell Cell for the point particle
	/// @param [in] lIndex Index of the point particle in the cell
	void addAtom(const uint64_t gIndex, const uint cell, const uint lIndex) {
		Position pos(cell,lIndex);
		m_pointParticles[gIndex]=pos;
	}
	/// @brief Add an point particle from the molecule
	/// @param [in] gIndex Global index of the point particle
	void removeAtom(const uint64_t gIndex) {
		m_pointParticles.erase(gIndex);
	}
	void updatePosition();
	void updateVelocity();
	void checkForcesConservation();

	/// @brief Shift the position of an atom in its cell
	/// @param [in] gIndex Global index of the atom
	/// @param [in] newPos New position
	inline void shift(const uint64_t gIndex, const uint newPos) {
		m_pointParticles[gIndex].index()=newPos;
	}
	template <size_t align> void fillVectMBufferBond(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& bonds);
	template <size_t align> void fillVectMBufferAngle(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& angles);
	template <size_t align> void fillVectMBufferDihedral(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& dihedrals);
	template <size_t align> void fillVectMBufferImproper(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& impropers);
	template <size_t align> void writeBondForces(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& bonds);
	template <size_t align> void writeAngleForces(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& angles);
	template <size_t align> void writeDihedralForces(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& dihedrals);
	template <size_t align> void writeImproperForces(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& impropers);
	template <size_t align> void writeBondViriel(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& bonds);
	template <size_t align> void writeAngleViriel(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& angles);
	template <size_t align> void writeDihedralViriel(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& dihedrals);
	template <size_t align> void writeImproperViriel(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& impropers);

protected :

	/// @brief Get the bonds for a point particle of the molecule
	/// @param [in] gIndex Global index of the particle
	/// @return Indexes vector of the particles bonded
	inline std::vector<uint64_t> getBonds(const uint64_t gIndex) {
		return m_cells[getCell(gIndex)].getBonds(getLocal(gIndex));
	}
	/// @brief Get the bond order between two particle of the molecule
	/// @param [in] gIndex Global index of the first particle
	/// @param [in] nbrIndex Global index of the second particle
	/// @return Bond order (integer, -1 for a partial bond order)
	inline int getBondOrder(const uint64_t gIndex, const uint64_t nbrIndex) {
		return m_cells[getCell(gIndex)].getBondOrder(getLocal(gIndex), nbrIndex);
	}
	/// @brief Get the subtype for a point particle of the molecule
	/// @param [in] gIndex Global index of the particle
	/// @return Subtype
	inline uint8_t getSubtype(const uint64_t gIndex) {
		return m_cells[getCell(gIndex)].getSubtype(getLocal(gIndex));
	}
	/// @brief Get the cell of a point particle
	/// @param [in] gIndex Global index of the point particle
	/// @return Cell
	inline uint getCell(const uint64_t gIndex) {
		return m_pointParticles[gIndex].cell();
	}
	/// @brief Get the local index of a point particle
	/// @param [in] gIndex Global index of the point particle
	/// @return Local index
	inline uint getLocal(const uint64_t gIndex) {
		return m_pointParticles[gIndex].index();
	}
	/// @brief Check if a point particle is in the ghosts
	/// @param [in] gIndex Index of the particle
	/// @return Ghost flag
	inline bool isGhost(const uint64_t gIndex) {
		return m_cells[getCell(gIndex)].isGhost();
	}

	/// @brief Check if a point particle not on domain edges
	/// @param [in] gIndex Index of the particle
	/// @return Inside flag
	inline bool isInside(const uint64_t gIndex) {
		return m_cells[getCell(gIndex)].isInside();
	}

	template <size_t align> void fillVectMBufferParam(const std::string typeDof, VectMBuffer<align>& buffer, const uint posInBuffer);

	std::unordered_map<uint64_t, Position> m_pointParticles; ///< Location for each point particle of the molecule
	std::vector<uint64_t> m_allIDs; ///< Id of each point particle of the molecule in a vector
	uint m_idUsed; ///< Count the atoms whose dofs have already been built
	CList*  m_cells; ///< Pointer to the cellLists from the grid

};


#include "particle/molecule/molecule.hxx"


#undef TMPLM
#undef TMPL_Mol

#endif /* MOLECULE_HPP_ */

