/*
 * molecule.hxx
 *
 *  Created on: Apr 28, 2016
 *      Author: giarda
 */

/// @file
/// @brief Implementations for the class Molecule


#include "globals.hpp"


/// @brief Constructor
/// @tparam CList Class of the point particle inside the cellList
/// @tparam elem_chunk (not used)
/// @tparam align Alignment size for vectorization
/// @param [in] index Index of the particle
/// @param [in] ptrCellList Where the atoms are
TMPLM TMPL_Molecule::Molecule(const uint64_t index, CList* ptrCellList) :
		m_pointParticles(),
		m_allIDs(),
		m_idUsed(0),
		m_cells(ptrCellList)
{
	id=index;
	ti=TYPEMOL;
}


/// @brief Destructor (nothing to do)
TMPLM TMPL_Molecule::~Molecule() {}


/// @brief Build a vector with the id of the molecule non-ghost atoms (if not done)
///
///
TMPLM void TMPL_Molecule::initBuilding() {
	if(m_idUsed==0) {
		for(auto it=m_pointParticles.begin(); it!=m_pointParticles.end(); ++it) {
			if(!isGhost(it->first)) m_allIDs.push_back(it->first);
		}
  }
}


/// @brief Clean the vector with the id of the molecule atoms (if dof use complete)
///
///
TMPLM void TMPL_Molecule::finishedWithDofs(const uint size) {
	m_idUsed+=size;
	if(m_idUsed==m_allIDs.size()) {
		m_allIDs.clear();
		m_idUsed=0;
	}
}


/// @brief Build the dofs for a subset of atoms of the molecule
///
/// Use only when ghost have been added
/// @param [in] begin First atom to work on
/// @param [in] size Number of atoms to work on
/// @param [in,out] bonds To store the bonds built
/// @param [in,out] angles To store the angles built
/// @param [in,out] dihedrals To store the dihedral angles built
/// @param [in,out] impropers To store the improper torsions built
TMPLM void TMPL_Molecule::buildDofs(const uint begin, const uint size,
		std::vector<Dof>& bonds, std::vector<Dof>& angles, std::vector<Dof>& dihedrals, std::vector<Dof>& impropers) {
	// Functions to get an atom from its index
	auto strSubtype = [&] (uint64_t id) -> std::string {return Global::reference.findSubtype(getSubtype(id));};
	auto bondOrder = [&] (uint64_t id1, uint64_t id2) -> int {return getBondOrder(id1, id2);};
	auto gBonds = [&] (uint64_t id) -> std::vector<uint64_t> {return getBonds(id);};
	// For each point particles of the subset
  for(auto nAt=begin; nAt<begin+size; ++nAt) {
		uint64_t indexi(m_allIDs[nAt]);
		// Build the dof only for the real particles
		if(!isGhost(indexi)) {
			std::vector<uint64_t> bondsi=getBonds(indexi);
			// Loop on the bonded atoms
			for(uint j(0); j<bondsi.size(); ++j) {
				uint64_t indexj=bondsi[j];
				std::vector<uint64_t> bondsj=getBonds(indexj);
				// For each unique bond i-j (bond considered only for the lower index atom)
				if(indexj>indexi||isGhost(indexj)){
					// Make a bond dof
					Dof bond(indexi,indexj);
					bond.setType(Global::ffield->type(), strSubtype, bondOrder, gBonds);
					if(!is_equal(bond.type(),"")) bonds.push_back(bond);
					// Look over neighbor of i and j to get the torsion k-i-j-l
					for(uint k(0); k<bondsi.size(); ++k) { for(uint l(0); l<bondsj.size(); ++l) {
						if(bondsi[k]!=indexj&&bondsj[l]!=indexi&&bondsi[k]!=bondsj[l]) {
							// Make a torsion dof
							Dof torsion(bondsi[k],indexi,indexj,bondsj[l]);
							torsion.setType(Global::ffield->type(), strSubtype, bondOrder, gBonds);
							if(!is_equal(torsion.type(),"")) {
								dihedrals.push_back(torsion);
							}
						}
					} }// End work on torsion
				}// End work on unique bond
				// Look each pair of neighbors of i to create the angle j-i-k
				for(uint k(j+1); k<bondsi.size(); ++k) {
					// Make an angle dof
					Dof angle(indexj,indexi,bondsi[k]);
					angle.setType(Global::ffield->type(), strSubtype, bondOrder, gBonds);
					if(!is_equal(angle.type(),"")) angles.push_back(angle);
					// Look one more neighbor of i to create the improper torsion i-jkl
					for(uint l(k+1); l<bondsi.size(); ++l) {
						// Make an improper torsion dof
						Dof improper(indexi,indexj,bondsi[k],bondsi[l],true);
						improper.setType(Global::ffield->type(), strSubtype, bondOrder, gBonds);
						if(!is_equal(improper.type(),"")) impropers.push_back(improper);
					}// End work on improper torsion
				}// End work on angle
				// Look at angles and improper torsions where central atom is in the ghost and dihedral angles
				// where central bond is in the ghost
				if(isGhost(indexj)) {
					// For each other neighbor of the ghost central atom
					for(uint k(0); k<bondsj.size();++k) {
						uint64_t indexk=bondsj[k];
						std::vector<uint64_t> bondsk=getBonds(indexk);
						// For angles and impropers, consider the other neighbor only if
						// it's a ghost or strictly larger than indexi (for uniqueness)
						if(indexk>indexi||isGhost(indexk)){
							// Make an angle dof
							Dof angle(indexi,indexj,indexk);
							angle.setType(Global::ffield->type(), strSubtype, bondOrder, gBonds);
							if(!is_equal(angle.type(),"")) angles.push_back(angle);
							// Look one more neighbor of j to create the improper torsion j-ikl
							// (Only considered if ghost or strictly larger than indexi for uniqueness)
							for(uint l(k+1); l<bondsj.size(); ++l) { if(bondsj[l]>indexi||isGhost(bondsj[l])) {
								// Make an improper torsion dof
								Dof improper(indexj,indexi,indexk,bondsj[l],true);
								improper.setType(Global::ffield->type(), strSubtype, bondOrder, gBonds);
								if(!is_equal(improper.type(),"")) impropers.push_back(improper);
							} }// End work on improper torsion
						}// End work on angles
						// For dihedrals, consider the other neighbor only is it's a ghost
						if(isGhost(indexk)) {
							// Look one more neighbor of k to create the dihedral angle ijkl
							// (Only considered if ghost or strictly larger than indexi for uniqueness)
							for(uint l(0); l<bondsk.size(); ++l) { if((bondsk[l]>indexi||isGhost(bondsk[l]))&&bondsk[l]!=indexj) {
								// Make a tosrion dof
								Dof torsion(indexi,indexj,indexk,bondsk[l]);
								torsion.setType(Global::ffield->type(), strSubtype, bondOrder, gBonds);
								if(!is_equal(torsion.type(),"")) {
									dihedrals.push_back(torsion);
								}
							} }// End work k neighbors
						}// End work on dihedrals
					}// End loop on second neighbor of ghost central
				}// End work on ghost central
			}// End loop on neighbors
		}// End not ghost condition
	}// End loop on atoms
}


/// @brief Update the celllist by deleting bonded neighbors
///
/// Those should not interact through Van der Waals or Coulomb interactions
/// This function will also check that all connected neighbors are in the cell or ghosts
/// @tparam level Indicate on which cells the correction is to be done :
/// 0 for all the cells, 1 for the inside cells only, 2 for the edge cells only
TMPLM template <uint level>  void TMPL_Molecule::updateNbr() {

	// Neighbor of a point particle will be uptated if it's not a ghost and either level is 0 or level is 1 and
	// the particle is inside the domain or level is 2 (not 1) and the particle is on the edge (not inside)
	auto consider = [&] (uint64_t id) -> bool {return !isGhost(id)&&(level==0||((level==1)==(isInside(id))));};

	// For each non-ghost atom
  for(auto it=m_pointParticles.begin(); it!=m_pointParticles.end(); ++it)	{
  	if(consider(it->first)) {
		// Useful variables
		std::tuple<uint, uint16_t, uint8_t>* nbrList;
		uint nbrSize;
		std::set<uint64_t> allLayers, currentLayer, addLayer;
		// Search for bonded neighbors
		addLayer.insert(it->first);
		// For each layer until the third
		for(uint layer(0); layer<3; ++layer) {
			// Update storages
			allLayers.insert(currentLayer.begin(),currentLayer.end());
			currentLayer=addLayer;
			addLayer.clear();
			// Get the next neighbors
			for(uint64_t particle : currentLayer) for(uint64_t nbr : getBonds(particle)) {
				if(m_pointParticles.count(nbr)!=0) addLayer.insert(nbr);
				// Stop if the neighbor is not in the grid
				else {
					std::cerr << "Error in Molecule::updateNbr : atom " << nbr << " bound to atom " << it->first;
					std::cerr << " should be on the domain. Please check that the distance between this two atoms";
					std::cerr << " has not exceeded the physical limit set in the input file." << std::endl;
					exit(1);
				}
			}
			// Delete those atoms from the neighbors list
			if(layer<Global::ffield->deletedNeighbors()) for(uint64_t nbr : addLayer) {
				// Get the neighbor list
				m_cells[getCell(it->first)].neighborList.getNeighbors(getLocal(it->first), nbrList, nbrSize);
				// Search for all corresponding indexes in the neighbors list and delete those neighbors
				for(uint i(0); i<nbrSize; ++i) {
					uint64_t nbrIndex=m_cells[std::get<0>(nbrList[i])].id[std::get<1>(nbrList[i])];
					m_cells[getCell(it->first)].lock(getLocal(it->first));
						if(nbrIndex==nbr) {
							m_cells[getCell(it->first)].neighborList.delNeighbor(getLocal(it->first),i);
						}
					m_cells[getCell(it->first)].unlock(getLocal(it->first));
				}
			}// Delete done
		}// End loop on layers
  	}
  }// End loop on atoms
}


/// @brief Get the nth layer of bonded neighbors for a point particle
///
/// The nth layer of bonded neighbors for a particle are the particles that
/// are linked to that particle by a path of n bonds and no path of fewer bonds
/// @param [in] start Index of the target point particle
/// @param [out] neighbors Index of the point particles in the layer
/// @param [in] level Level of the layer
TMPLM void TMPL_Molecule::getNeighbors(const uint64_t start, std::set<uint64_t>& neighbors, const uint level) {
	for(uint64_t nbr : getBonds(start)) neighbors.insert(nbr);
	std::set<uint64_t> oldLayers;
	oldLayers.insert(start);
	for(uint i(1); i<level; ++i) {
		std::set<uint64_t> oldLayer=neighbors;
		oldLayers.insert(neighbors.begin(),neighbors.end());
		neighbors.clear();
		for(uint64_t pointParticle : oldLayer) for(uint64_t nbr : getBonds(pointParticle)) {
			if(oldLayers.count(nbr)==0) neighbors.insert(nbr);
		}
	}
}


/// @brief Get then number of particles in the nth layer of
/// bonded neighbors for a point particle
///
/// The nth layer of bonded neighbors for a particle are the particles that
/// are linked to that particle by a path of n bonds and no path of fewer bonds
/// @param [in] start Index of the target point particle
/// @param [in] level Level of the layer
TMPLM uint TMPL_Molecule::getNumberOfNeighbors(const uint64_t start, const uint level) {
	std::set<uint64_t> neighbors;
	getNeighbors(start, neighbors, level);
	return neighbors.size();
}


/// @brief Update the position center of mass of the molecule according
/// to the position of it's atoms
///
///
TMPLM void TMPL_Molecule::updatePosition(){
	double totalMass(0), mass(0);
	this->position()=0;
  for(auto it=m_pointParticles.begin(); it!=m_pointParticles.end(); ++it) {
		mass=Global::reference.getMass(m_cells[it->second.cell()].ti[it->second.index()]);
		totalMass+=mass;
		this->position().x+=mass*m_cells[it->second.cell()].rx[it->second.index()];
		this->position().y+=mass*m_cells[it->second.cell()].ry[it->second.index()];
		this->position().z+=mass*m_cells[it->second.cell()].rz[it->second.index()];
	}
	this->position()/=totalMass;
}


/// @brief Update the velocity of the molecule according to the velocity of it's atoms
///
///
TMPLM void TMPL_Molecule::updateVelocity(){
	double totalMass(0), mass(0);
	this->velocity()=0;
  for(auto it=m_pointParticles.begin(); it!=m_pointParticles.end(); ++it) {
		mass=Global::reference.getMass(m_cells[it->second.cell()].ti[it->second.index()]);
		totalMass+=mass;
		this->velocity().x+=mass*m_cells[it->second.cell()].vx[it->second.index()];
		this->velocity().y+=mass*m_cells[it->second.cell()].vy[it->second.index()];
		this->velocity().z+=mass*m_cells[it->second.cell()].vz[it->second.index()];
	}
	this->velocity()/=totalMass;
}


/// @brief Calculate and print the total force on the molecule
///
///
TMPLM void TMPL_Molecule::checkForcesConservation(){
	double mass(0);
	vec3<double> force(0);
	this->velocity()=0;
  for(auto it=m_pointParticles.begin(); it!=m_pointParticles.end(); ++it) {
		mass=Global::reference.getMass(m_cells[it->second.cell()].ti[it->second.index()]);
		force.x+=mass*m_cells[it->second.cell()].fx[it->second.index()];
		force.y+=mass*m_cells[it->second.cell()].fy[it->second.index()];
		force.z+=mass*m_cells[it->second.cell()].fz[it->second.index()];
	}
	std::cerr << "forces on molecule " << this->id << " : " << force << std::endl;
}


/// @brief Fill the force buffer with the positions of the atoms for the specified bonds
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of bonds to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first bond to consider
/// @param [in] bonds Vector with the bonds to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::fillVectMBufferBond(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& bonds) {
	for(uint i(0); i<size; ++i) {
	  uint cell=getCell(bonds[i+posInDofs].index(0));
	  uint position=getLocal(bonds[i+posInDofs].index(0));
		buffer.rx1(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry1(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz1(posInBuffer+i)=m_cells[cell].rz[position];
	  cell=getCell(bonds[i+posInDofs].index(1));
	  position=getLocal(bonds[i+posInDofs].index(1));
		buffer.rx1(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry1(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz1(posInBuffer+i)-=m_cells[cell].rz[position];
		fillVectMBufferParam(bonds[i+posInDofs].type(), buffer, posInBuffer+i);
	}
}


/// @brief Fill the force buffer with the positions of the atoms for the specified angles
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of angles to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first angle to consider
/// @param [in] angles Vector with the angles to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::fillVectMBufferAngle(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& angles) {
	for(uint i(0); i<size; ++i) {
	  uint cell=getCell(angles[i+posInDofs].index(0));
	  uint position=getLocal(angles[i+posInDofs].index(0));
		buffer.rx1(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry1(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz1(posInBuffer+i)=m_cells[cell].rz[position];
	  cell=getCell(angles[i+posInDofs].index(2));
	  position=getLocal(angles[i+posInDofs].index(2));
		buffer.rx2(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry2(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz2(posInBuffer+i)=m_cells[cell].rz[position];
	  cell=getCell(angles[i+posInDofs].index(1));
	  position=getLocal(angles[i+posInDofs].index(1));
		buffer.rx1(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry1(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz1(posInBuffer+i)-=m_cells[cell].rz[position];
		buffer.rx2(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry2(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz2(posInBuffer+i)-=m_cells[cell].rz[position];
		fillVectMBufferParam(angles[i+posInDofs].type(), buffer, posInBuffer+i);
	}
}


/// @brief Fill the force buffer with the positions of the atoms for the specified dihedral angles
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of dihedral angles to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first dihedral angle to consider
/// @param [in] dihedrals Vector with the dihedral angles to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::fillVectMBufferDihedral(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& dihedrals) {
	for(uint i(0); i<size; ++i) {
	  uint cell=getCell(dihedrals[i+posInDofs].index(0));
	  uint position=getLocal(dihedrals[i+posInDofs].index(0));
		buffer.rx1(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry1(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz1(posInBuffer+i)=m_cells[cell].rz[position];
	  cell=getCell(dihedrals[i+posInDofs].index(3));
	  position=getLocal(dihedrals[i+posInDofs].index(3));
		buffer.rx3(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry3(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz3(posInBuffer+i)=m_cells[cell].rz[position];
	  cell=getCell(dihedrals[i+posInDofs].index(2));
	  position=getLocal(dihedrals[i+posInDofs].index(2));
		buffer.rx2(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry2(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz2(posInBuffer+i)=m_cells[cell].rz[position];
		buffer.rx3(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry3(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz3(posInBuffer+i)-=m_cells[cell].rz[position];
	  cell=getCell(dihedrals[i+posInDofs].index(1));
	  position=getLocal(dihedrals[i+posInDofs].index(1));
		buffer.rx1(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry1(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz1(posInBuffer+i)-=m_cells[cell].rz[position];
		buffer.rx2(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry2(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz2(posInBuffer+i)-=m_cells[cell].rz[position];
		fillVectMBufferParam(dihedrals[i+posInDofs].type(), buffer, posInBuffer+i);
	}
}


/// @brief Fill the force buffer with the positions of the atoms for the specified improper torsions
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of improper torsions to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first improper torsion to consider
/// @param [in] impropers Vector with the improper torsions to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::fillVectMBufferImproper(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& impropers) {
	for(uint i(0); i<size; ++i) {
		uint cell=getCell(impropers[i+posInDofs].index(0));
		uint position=getLocal(impropers[i+posInDofs].index(0));
		buffer.rx1(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry1(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz1(posInBuffer+i)=m_cells[cell].rz[position];
		buffer.rx2(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry2(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz2(posInBuffer+i)=m_cells[cell].rz[position];
		buffer.rx3(posInBuffer+i)=m_cells[cell].rx[position];
		buffer.ry3(posInBuffer+i)=m_cells[cell].ry[position];
		buffer.rz3(posInBuffer+i)=m_cells[cell].rz[position];
	  cell=getCell(impropers[i+posInDofs].index(1));
	  position=getLocal(impropers[i+posInDofs].index(1));
		buffer.rx1(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry1(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz1(posInBuffer+i)-=m_cells[cell].rz[position];
	  cell=getCell(impropers[i+posInDofs].index(2));
	  position=getLocal(impropers[i+posInDofs].index(2));
		buffer.rx2(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry2(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz2(posInBuffer+i)-=m_cells[cell].rz[position];
	  cell=getCell(impropers[i+posInDofs].index(3));
	  position=getLocal(impropers[i+posInDofs].index(3));
		buffer.rx3(posInBuffer+i)-=m_cells[cell].rx[position];
		buffer.ry3(posInBuffer+i)-=m_cells[cell].ry[position];
		buffer.rz3(posInBuffer+i)-=m_cells[cell].rz[position];
		fillVectMBufferParam(impropers[i+posInDofs].type(), buffer, posInBuffer+i);
	}
}


/// @brief Fiff the force buffer with the parameters for a dof
/// @param [in] typeDof Type of the dof
/// @param [in,out] buffer Buffer to fill
/// @param [in] posInBuffer Position to fill
TMPLM template <size_t align> void TMPL_Molecule::fillVectMBufferParam(const std::string typeDof, VectMBuffer<align>& buffer, const uint posInBuffer){
	std::vector<double> params;
	Global::ffield->getParam(typeDof, params);
	switch(params.size()) {
	case 4 :
		buffer.p4(posInBuffer)=params[3];
	case 3 :
		buffer.p3(posInBuffer)=params[2];
	case 2 :
		buffer.p2(posInBuffer)=params[1];
	case 1 :
		buffer.p1(posInBuffer)=params[0];
	default :
		break;
	}
}


/// @brief Write the forces back into the cellList (for bonds)
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of bonds to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first bond to consider
/// @param [in] bonds Vector with the bonds to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::writeBondForces(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& bonds) {
	for(uint i(0); i<size; ++i) {
	  uint cell=getCell(bonds[i+posInDofs].index(0));
	  uint position=getLocal(bonds[i+posInDofs].index(0));
	  double invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the first particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx1(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy1(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz1(posInBuffer+i);
		m_cells[cell].ep[position]+=.5*buffer.ep(posInBuffer+i);
		m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the second particle
	  cell=getCell(bonds[i+posInDofs].index(1));
	  position=getLocal(bonds[i+posInDofs].index(1));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the second particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]-=invMass*buffer.fx1(posInBuffer+i);
		m_cells[cell].fy[position]-=invMass*buffer.fy1(posInBuffer+i);
		m_cells[cell].fz[position]-=invMass*buffer.fz1(posInBuffer+i);
		m_cells[cell].ep[position]+=.5*buffer.ep(posInBuffer+i);
	  m_cells[cell].unlock(position);
	}
}


/// @brief Write the forces back into the cellList (for angles)
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of angles to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first angle to consider
/// @param [in] angles Vector with the angles to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::writeAngleForces(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& angles){
	for(uint i(0); i<size; ++i) {
	  // Get the cell, position and inverse mass for the right particle
	  uint cell=getCell(angles[i+posInDofs].index(0));
	  uint position=getLocal(angles[i+posInDofs].index(0));
	  double invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the first particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx1(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy1(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz1(posInBuffer+i);
		m_cells[cell].ep[position]+=buffer.ep(posInBuffer+i)/3.;
		m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the left particle
	  cell=getCell(angles[i+posInDofs].index(2));
	  position=getLocal(angles[i+posInDofs].index(2));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the third particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx2(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy2(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz2(posInBuffer+i);
		m_cells[cell].ep[position]+=buffer.ep(posInBuffer+i)/3.;
	  m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the middle particle
	  cell=getCell(angles[i+posInDofs].index(1));
	  position=getLocal(angles[i+posInDofs].index(1));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the second particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]-=invMass*(buffer.fx1(posInBuffer+i)+buffer.fx2(posInBuffer+i));
		m_cells[cell].fy[position]-=invMass*(buffer.fy1(posInBuffer+i)+buffer.fy2(posInBuffer+i));
		m_cells[cell].fz[position]-=invMass*(buffer.fz1(posInBuffer+i)+buffer.fz2(posInBuffer+i));
		m_cells[cell].ep[position]+=buffer.ep(posInBuffer+i)/3.;
	  m_cells[cell].unlock(position);
	}
}


/// @brief Write the forces back into the cellList (for dihedrals)
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of dihedrals to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first dihedral angle to consider
/// @param [in] dihedrals Vector with the dihedral angles to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::writeDihedralForces(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& dihedrals){
	for(uint i(0); i<size; ++i) {
	  // Get the cell, position and inverse mass for the first particle
		uint cell=getCell(dihedrals[i+posInDofs].index(0));
	  uint position=getLocal(dihedrals[i+posInDofs].index(0));
	  double invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the first particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx1(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy1(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz1(posInBuffer+i);
		m_cells[cell].ep[position]+=0.25*buffer.ep(posInBuffer+i);
		m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the second particle
	  cell=getCell(dihedrals[i+posInDofs].index(1));
	  position=getLocal(dihedrals[i+posInDofs].index(1));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the second particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx2(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy2(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz2(posInBuffer+i);
		m_cells[cell].ep[position]+=0.25*buffer.ep(posInBuffer+i);
	  m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the third particle
	  cell=getCell(dihedrals[i+posInDofs].index(3));
	  position=getLocal(dihedrals[i+posInDofs].index(3));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the second particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx3(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy3(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz3(posInBuffer+i);
		m_cells[cell].ep[position]+=0.25*buffer.ep(posInBuffer+i);
	  m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the forth particle
	  cell=getCell(dihedrals[i+posInDofs].index(2));
	  position=getLocal(dihedrals[i+posInDofs].index(2));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the second particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]-=invMass*(buffer.fx1(posInBuffer+i)+buffer.fx2(posInBuffer+i)+buffer.fx3(posInBuffer+i));
		m_cells[cell].fy[position]-=invMass*(buffer.fy1(posInBuffer+i)+buffer.fy2(posInBuffer+i)+buffer.fy3(posInBuffer+i));
		m_cells[cell].fz[position]-=invMass*(buffer.fz1(posInBuffer+i)+buffer.fz2(posInBuffer+i)+buffer.fz3(posInBuffer+i));
		m_cells[cell].ep[position]+=0.25*buffer.ep(posInBuffer+i);
	  m_cells[cell].unlock(position);
	}
}


/// @brief Write the forces back into the cellList (for dihedrals)
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of dihedrals to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first dihedral angle to consider
/// @param [in] impropers Vector with the improper torsions to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::writeImproperForces(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& impropers){
	for(uint i(0); i<size; ++i) {
	  // Get the cell, position and inverse mass for the central particle
		uint cell=getCell(impropers[i+posInDofs].index(0));
		uint position=getLocal(impropers[i+posInDofs].index(0));
		double invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the second particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]-=invMass*(buffer.fx1(posInBuffer+i)+buffer.fx2(posInBuffer+i)+buffer.fx3(posInBuffer+i));
		m_cells[cell].fy[position]-=invMass*(buffer.fy1(posInBuffer+i)+buffer.fy2(posInBuffer+i)+buffer.fy3(posInBuffer+i));
		m_cells[cell].fz[position]-=invMass*(buffer.fz1(posInBuffer+i)+buffer.fz2(posInBuffer+i)+buffer.fz3(posInBuffer+i));
		m_cells[cell].ep[position]+=0.25*buffer.ep(posInBuffer+i);
	  m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the second particle
	  cell=getCell(impropers[i+posInDofs].index(1));
	  position=getLocal(impropers[i+posInDofs].index(1));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the second particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx1(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy1(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz1(posInBuffer+i);
		m_cells[cell].ep[position]+=0.25*buffer.ep(posInBuffer+i);
	  m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the third particle
	  cell=getCell(impropers[i+posInDofs].index(2));
	  position=getLocal(impropers[i+posInDofs].index(2));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the second particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx2(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy2(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz2(posInBuffer+i);
		m_cells[cell].ep[position]+=0.25*buffer.ep(posInBuffer+i);
	  m_cells[cell].unlock(position);
	  // Get the cell, position and inverse mass for the forth particle
	  cell=getCell(impropers[i+posInDofs].index(3));
	  position=getLocal(impropers[i+posInDofs].index(3));
	  invMass = Global::reference.getInvMass(m_cells[cell].ti[position]);
	  // Write the forces for the first particle
	  m_cells[cell].lock(position);
		m_cells[cell].fx[position]+=invMass*buffer.fx3(posInBuffer+i);
		m_cells[cell].fy[position]+=invMass*buffer.fy3(posInBuffer+i);
		m_cells[cell].fz[position]+=invMass*buffer.fz3(posInBuffer+i);
		m_cells[cell].ep[position]+=0.25*buffer.ep(posInBuffer+i);
		m_cells[cell].unlock(position);
	}
}


/// @brief Write the viriel back into the cellList
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of bonds to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first bond to consider
/// @param [in] bonds Vector with the bonds to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::writeBondViriel(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& bonds) {
	for(uint i(0); i<size; ++i) {
	  // TODO
	}
}


/// @brief Write the viriel back into the cellList
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of bonds to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first bond to consider
/// @param [in] angles Vector with the angles to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::writeAngleViriel(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& angles){
	for(uint i(0); i<size; ++i) {
	  // TODO
	}
}


/// @brief Write the viriel back into the cellList
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of bonds to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first bond to consider
/// @param [in] dihedrals Vector with the dihedral angles to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::writeDihedralViriel(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& dihedrals){
	for(uint i(0); i<size; ++i) {
	  // TODO
	}
}


/// @brief Write the viriel back into the cellList
/// @param [in,out] buffer Buffer to fill
/// @param [in] size Number of bonds to complete
/// @param [in] posInBuffer Position of the data in the buffer
/// @param [in] posInDofs Position of the first bond to consider
/// @param [in] impropers Vector with the improper torsions to consider (and others)
TMPLM template <size_t align> void TMPL_Molecule::writeImproperViriel(VectMBuffer<align>& buffer, const uint size, const uint posInBuffer, const uint posInDofs, const std::vector<Dof>& impropers){
	for(uint i(0); i<size; ++i) {
	  // TODO
	}
}
