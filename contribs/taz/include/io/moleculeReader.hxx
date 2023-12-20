/*
 * moleculeReader.hxx
 *
 *  Created on: 10 févr. 2016
 *      Author: giarda
 */

/// @file
/// @brief Implementations for the class MoleculeReaderOB


/// @brief Constructor
/// @param [in] fileName Initialization file name
/// @param [in] fileDir Initialization file directory
/// @param [in] format Format for the initialization file
/// @param [in] multiplier Number of cells in the system (to get the total number of atoms)
/// @param [in] add Flag indicating if hydrogens must be added to the system, default=false
/// @param [in] conformer Number of the conformer to read in the initialization file, default=1
/// @param [in] rank MPI rank
TMPLMR MoleculeReaderOB<PointParticle,ff,chM>::MoleculeReaderOB(std::string fileName, std::string fileDir, std::string format, int multiplier, bool add, int conformer, int rank) :
		m_fileName(fileName),
		m_fileDir(fileDir),
		m_format(format),
		m_addH(add),
		m_atoms()
	{
	readFile(multiplier, conformer, rank);
	build();
}


/// @brief On big loop on all atoms of all molecules to build the system
///
/// This will build and store :
///  - all atoms
///  - system dimensions
///  - list of all subtypes
TMPLMR void MoleculeReaderOB<PointParticle,ff,chM>::build() {
	// Initialize the subtypes
	std::string currentSubtype;
	// Get the position of the first atom
	m_min.x=convert(m_molecules[0].GetFirstAtom()->GetX(),SI_Units_base::angstrom,Stamp_Units::length);
	m_min.y=convert(m_molecules[0].GetFirstAtom()->GetY(),SI_Units_base::angstrom,Stamp_Units::length);
	m_min.z=convert(m_molecules[0].GetFirstAtom()->GetZ(),SI_Units_base::angstrom,Stamp_Units::length);
	m_max=m_min;
	// Initialize the number of true molecules (not momoatomic) and
	// the molecules index
	uint64_t nbMol(0);
	uint indexMol(0);
	// For each molecule
	for (uint64_t iMol=0;iMol<m_molecules.size();iMol++) {
		// Distinguish between monoatomic and true molecules
		if(m_molecules[iMol].NumAtoms()==1) {
			indexMol=0;
		}
		else {
			++nbMol;
			indexMol=nbMol+MPI__InMol::s_nbPointParticles;
		}
		// For each atom
		FOR_ATOMS_OF_MOL(atom,m_molecules[iMol]) {
			// Get atom data
			// Pointer
			MPI__PInMol<PointParticle>* atomToComplete(&m_atoms[index(atom->GetIdx(),iMol)-1]);
			// Index
			atomToComplete->index()=index(atom->GetIdx(),iMol);
			// Position
			atomToComplete->position().x=convert(atom->GetX(),SI_Units_base::angstrom,Stamp_Units::length);
			atomToComplete->position().y=convert(atom->GetY(),SI_Units_base::angstrom,Stamp_Units::length);
			atomToComplete->position().z=convert(atom->GetZ(),SI_Units_base::angstrom,Stamp_Units::length);
			// Charge
			m_chargeReader.get(atomToComplete, atom);
			// Molecule index
			atomToComplete->molecule()=indexMol;
			// Bonds
			for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
				atomToComplete->nextBond()=index((*itBond)->GetNbrAtom(&(*atom))->GetIdx(),iMol);
			}
			if(ff==FF::UFF) makeBondOrders(atomToComplete,atom,iMol);

			// Update min and max coordinate
			m_min=auxMin(m_min,vec3<double>(atomToComplete->position().x,atomToComplete->position().y, atomToComplete->position().z));
			m_max=auxMax(m_max,vec3<double>(atomToComplete->position().x,atomToComplete->position().y, atomToComplete->position().z));

			// Get subtype
			currentSubtype=m_typeReader.get(atom)+"$"+std::to_string(atom->GetAtomicNum())+"$"+std::to_string(atom->GetExactMass());
			bool there(false);
			// Check if this subtype is in the list. If already there, get index
			for(uint8_t iType(0); iType<m_listOfSubtypes.size(); ++iType) {
				if(is_equal(m_listOfSubtypes[iType],currentSubtype)) {
					there=true;
					m_atoms[index(atom->GetIdx(),iMol)-1].subtype()=iType;
					break;
				}
			}
			// Else add to the list of subtypes
			if(!there) {
				atomToComplete->subtype()=m_listOfSubtypes.size();
				m_listOfSubtypes.push_back(currentSubtype);
			}

		} // End of atoms
	} // End of molecules
	m_nMol=nbMol;
}


/// @brief Get the size of the system
/// @param [out] minBounds Lower boundaries
/// @param [out] maxBounds Upper boundaries
TMPLMR void MoleculeReaderOB<PointParticle,ff,chM>::getSystemDimensions(vec3<double>& minBounds, vec3<double>& maxBounds) {
	minBounds=m_min;
	maxBounds=m_max;
}


/// @brief Return the list all possible subtypes in the system
/// @return List
TMPLMR std::vector<std::string> MoleculeReaderOB<PointParticle, ff,chM>::getSubtypes() {
	return m_listOfSubtypes;
}


/// @brief Build the atoms
/// @tparam PointParticle Class of the atoms in the molecules
/// @return Vector of atoms in transferable form
TMPLMR void MoleculeReaderOB<PointParticle,ff,chM>::getAtoms(std::vector<MPI__Particle*>& atoms, uint64_t& numMol)
{
	for(MPI__PInMol<PointParticle> atomInMol : m_atoms) {
		atoms.push_back(new MPI__PInMol<PointParticle>(atomInMol));
	}
	numMol=m_nMol;
}


/// @brief Read the input file to get the molecules
///
/// Does also count the point particles in the system and split them
/// into disconnected molecules
/// @param [in] multiplier Number of cells in the system (to get the total number of atoms)
/// @param [in] conformer Number of the conformer to read in the file
/// @param [in] rank MPI rank
TMPLMR void MoleculeReaderOB<PointParticle,ff,chM>::readFile(int multiplier, int conformer, int rank) {
	//
	// File reading
	//
	// Get an OpenBabel molecular formats conversion tool and set the in and out formats
	OpenBabel::OBConversion converter;
	converter.SetInAndOutFormats(m_format.c_str(),m_format.c_str());
	// Set the input stream
	std::fstream file;
	std::string fileName=m_fileDir+m_fileName;
	file.open(fileName.c_str(),std::fstream::in);
	converter.SetInStream(&file);
	// New OpenBabel molecule and its number of atoms
	OpenBabel::OBMol newMol;
	uint64_t nbAtoms;
	// Read molecules until the specified conformer, keep only the last
	bool read=false;
	for(int i=0;i<conformer;i++) {
		read=converter.Read(&newMol);	}
	if(!read) {
		std::cerr << "Error in MoleculeReaderOB::readFile : molecule not read !" << std::endl;
		exit(-1);
	}
	// Close the file
	file.close();
	// Get the number of atoms
	nbAtoms=newMol.NumAtoms();
	//
	// Hydrogens addition
	//
	// If adding hydrogens is required
	if(m_addH) {
		// Store the number of atoms
		uint64_t oldNbAtoms=nbAtoms;
		// Add hydrogens
		newMol.AddHydrogens();
		// Get the new number of atoms
		nbAtoms=newMol.NumAtoms();
		// While the number of atoms keep changing
		while (oldNbAtoms!=nbAtoms) {
			// Print the molecule and read it again
			converter.WriteFile(&newMol,m_fileDir+"AddH"+std::to_string(rank)+"."+m_format);
			newMol.Clear();
			converter.ReadFile(&newMol,m_fileDir+"AddH"+std::to_string(rank)+"."+m_format);
			// Then store the number of atoms and add new hydrogens
			oldNbAtoms=nbAtoms;
			newMol.AddHydrogens();
			nbAtoms=newMol.NumAtoms();
		}
	}
	//
	// Point particle count
	//
	m_atoms.resize(newMol.NumAtoms());
	MPI__InMol::s_nbPointParticles=newMol.NumAtoms()*multiplier;
	//
	// Molecules separation
	//
	// Divide the molecule into its non-bonded fragments
	m_molecules=newMol.Separate();
	// Get a global index for each atom of each molecule (stored in m_indexes)
	uint64_t globalIndexAtom(0);
	for(uint64_t iMol(0);iMol<m_molecules.size();++iMol) { FOR_ATOMS_OF_MOL(atom,m_molecules[iMol]) {
			++globalIndexAtom;
			m_indexes[std::make_pair(iMol,atom->GetIdx())]=globalIndexAtom;
	} }
}


/// @brief Add data to an atom bonded list so that it's possible to get the bond order from it
/// @param [in,out] toComplete Pointer to the atom to complete
/// @param [in] atom Iterator to the Open Babel atom containing the data
/// @param [in] numMol Index of the current molecule
TMPLMR void MoleculeReaderOB<PointParticle,ff,chM>::makeBondOrders(MPI__InMol* toComplete, OpenBabel::OBMolAtomIter atom, uint64_t numMol){
	// For each bonded neighbor of the atom
	for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
		// If amide or aromatic bond add -indexOfTheNeighbor to the bonded list (non integral bond order)
		if((*itBond)->IsAmide()||(*itBond)->IsAromatic())
			toComplete->nextBond()=index((*itBond)->GetNbrAtom(&(*atom))->GetIdx(),numMol)+MPI__InMol::s_nbPointParticles;
		// Else add the index of the neighbor bondOrder-1 times so the number of occurrences of this index is now the bond order
		else for(uint i(1); i<(*itBond)->GetBondOrder(); ++i)
			toComplete->nextBond()=index((*itBond)->GetNbrAtom(&(*atom))->GetIdx(),numMol);
	} // End of the loop on neighbors
}
