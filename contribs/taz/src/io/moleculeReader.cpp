/*
 * typesReader.hxx
 *
 *  Created on: 25 févr. 2016
 *      Author: giarda
 */

/// @file
/// @brief Implementations for the class TypeReader


#include <iostream>

#include "io/moleculeReader.hpp"


	/// @brief Default constructor
Configuration<MPI__InMol>::Configuration() :
	m_forceField(FF::UFF),
	m_charges(FF::NONE),
	m_fileName(""),
	m_fileDir("./"),
	m_format("pdb"),
	m_add(false),
	m_conformer(1),
	m_margins(0.),
	m_initialTemperature(300),
	m_maxBond(1.05),
	m_nCells(1)
	{}


	/// @brief Constructor from the input
	/// @param [in] input Input data
Configuration<MPI__InMol>::Configuration(const Input& input) :
	m_fileName(input.initFile),
	m_fileDir(input.initDir+"/"),
	m_format(input.m_format),
	m_add(input.m_addHydrogens),
	m_conformer(input.m_numConformer),
	m_margins(input.m_margins),
	m_initialTemperature(input.initialTemperature),
	m_maxBond(input.m_maxBonds),
	m_nCells(input.numberOfCells)
	{
		switch(input.m_forceField) {
		case Input::UFF :
			m_forceField=FF::UFF;
			break;
		case Input::COMPASS :
			m_forceField=FF::COMPASS;
			break;
		case Input::AMBER :
			m_forceField=FF::AMBER;
			break;
		}
		switch(input.m_chargeMethod) {
		case Input::NOT_SPECIFIED :
			m_charges=ForceField::defaultCharge.at(m_forceField);
			break;
		case Input::NO_CHARGES :
			m_charges=FF::NONE;
			break;
		case Input::FROM_COMPASS :
			m_charges=FF::FROMCOMPASS;
			break;
		case Input::OPENBABEL :
			m_charges=FF::OPENBABEL;
			break;
		}
	}


#ifdef __use_lib_openbabel
/// @brief Get the subtype for an atom of a molecule (case UFF)
/// @param [in] atom Iterator on the atom
/// @return Subtype
template <> std::string TypeReader<FF::UFF>::get(OpenBabel::OBMolAtomIter atom) {
	std::string type="";
	// Subtype construction depend on the element
	switch(atom->GetAtomicNum()) {
	case 1 :
	{
		int nbBoron(0);
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
			uint8_t Zlinked=(*itBond)->GetNbrAtom(&(*atom))->GetAtomicNum();
			if(Zlinked==5) nbBoron++;
		}
		if(nbBoron>1) type="H_b";
		else type="H_";
		break;
	}
	case 2 :
		type="He4+4";
		break;
	case 3 :
		type="Li";
		break;
	case 4 :
		type="Be3+2";
		break;
	case 5 :
	{
		int nbBonds(0);
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) nbBonds++;
		if(nbBonds==4) type="B_3";
		else type="B_2";
		break;
	}
	case 6 :
		if(atom->IsAromatic()) type="C_R";
		else {
			int nbBonds(0);
			for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) nbBonds++;
			if(nbBonds>3) type="C_3";
			else if(nbBonds==3) type="C_2";
			else type="C_1";
		}
		break;
	case 7 :
		if(atom->IsAromatic()) type="C_R";
		else {
			int nbBonds(0);
			for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) nbBonds++;
			if(nbBonds>2) type="N_3";
			else if(nbBonds==2) type="N_2";
			else type="N_1";
		}
		break;
	case 8 :
		if(atom->IsAromatic()) type="O_R";
		else {
			int nbBonds(0), nbSilicon(0), nbTriple(0);
			for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
				nbBonds++;
				uint8_t Zlinked=(*itBond)->GetNbrAtom(&(*atom))->GetAtomicNum();
				if(Zlinked==14) nbSilicon++;
				if((*itBond)->GetBondOrder()==3) nbTriple++;
			}
			if(nbSilicon>1) type="O_3_z";
			else if(nbTriple>0) type="O_1";
			else if(nbBonds>1) type="O_3";
			else type="O_2";
		}
		break;
	case 9 :
		type="F_";
		break;
	case 10 :
		type="Ne4+4";
		break;
case 11 :
	type="Na";
	break;
	case 12 :
		type="Mg3+2";
		break;
	case 13 :
		type="Al3";
		break;
	case 14 :
		type="Si3";
		break;
	case 15 :
	{
		int nbOxidation(atom->GetFormalCharge());
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
			uint8_t Zlinked=(*itBond)->GetNbrAtom(&(*atom))->GetAtomicNum();
			if(Zlinked!=15) {
				bool PGetElec(true);
				for(uint8_t Z : {1,6,7,8,9,16,17,34,35,36,44,45,46,53,54,74,76,77,78,79,85,86}) {
					if(Zlinked==Z) PGetElec=false;
				}
				if(PGetElec) nbOxidation-=(*itBond)->GetBondOrder();
				else nbOxidation+=(*itBond)->GetBondOrder();
			}
		}
		if(nbOxidation==3) type="P_3+3";
		else if(nbOxidation==5) type="P_3+5";
		else type="P_3+q";
		break;
	}
	case 16 :
		if(atom->IsAromatic()) type="S_R";
		else {
			int nbDouble(0);
			int nbOxidation(atom->GetFormalCharge());
			for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
				if((*itBond)->GetBondOrder()==2) nbDouble++;
				uint8_t Zlinked=(*itBond)->GetNbrAtom(&(*atom))->GetAtomicNum();
				if(Zlinked!=16) {
					bool PGetElec(true);
					for(uint8_t Z : {7,8,9,17,35,36,53,54}) {
						if(Zlinked==Z) PGetElec=false;
					}
					if(PGetElec) nbOxidation-=(*itBond)->GetBondOrder();
					else nbOxidation+=(*itBond)->GetBondOrder();
				}
			}
			if(nbDouble==3) type="S_2";
			else if(nbOxidation<3) type="S_3+2";
			else if(nbOxidation<5) type="S_3+4";
			else type="S_3+6";
		}
		break;
	case 17 :
		type="Cl";
		break;
	case 18 :
		type="Ar4+4";
		break;
	case 19 :
		type="K_";
		break;
	case 20 :
		type="Ca6+2";
		break;
	case 21 :
		type="Sc3+3";
		break;
	case 22 :
	{
		int nbBonds=0;
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) nbBonds++;
		if(nbBonds<5) type="Ti3+4";
		else type="Ti6+4";
		break;
	}
	case 23 :
		type="V_3+5";
		break;
	case 24 :
		type="Cr6+3";
		break;
	case 25 :
		type="Mn6+2";
		break;
	case 26 :
	{
		int nbBonds=0;
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) nbBonds++;
		if(nbBonds<5) type="Fe3+2";
		else type="Fe6+2";
		break;
	}
	case 27 :
		type="Co6+3";
		break;
	case 28 :
		type="Ni4+2";
		break;
	case 29 :
		type="Cu3+1";
		break;
	case 30 :
		type="Zn3+2";
		break;
	case 31 :
		type="Ga3+3";
		break;
	case 32 :
		type="Ge3";
		break;
	case 33 :
		type="As3+3";
		break;
	case 34 :
		type="Se3+2";
		break;
	case 35 :
		type="Br";
		break;
	case 36 :
		type="Kr4+4";
		break;
	case 37 :
		type="Rb";
		break;
	case 38 :
		type="Sr6+2";
		break;
	case 39 :
		type="Y_3+3";
		break;
	case 40 :
		type="Zr3+4";
		break;
	case 41 :
		type="Nb3+5";
		break;
	case 42 :
	{
		int nbBonds=0;
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) nbBonds++;
		if(nbBonds<5) type="Mo3+6";
		else type="Mo6+6";
		break;
	}
	case 43 :
		type="Tc6+5";
		break;
	case 44 :
		type="Ru6+2";
		break;
	case 45 :
		type="Rh6+3";
		break;
	case 46 :
		type="Pd4+2";
		break;
	case 47 :
		type="Ag1+1";
		break;
	case 48 :
		type="Cd3+2";
		break;
	case 49 :
		type="In3+3";
		break;
	case 50 :
		type="Sn3";
		break;
	case 51 :
		type="Sb3+3";
		break;
	case 52 :
		type="Te3+2";
		break;
	case 53 :
		type="I_";
		break;
	case 54 :
		type="Xe4+4";
		break;
	case 55 :
		type="Cs";
		break;
	case 56 :
		type="Ba6+2";
		break;
	case 57 :
		type="La3+3";
		break;
	case 58 :
		type="Ce6+3";
		break;
	case 59 :
		type="Pr6+3";
		break;
	case 60 :
		type="Nd6+3";
		break;
	case 61 :
		type="Pm6+3";
		break;
	case 62 :
		type="Sm6+3";
		break;
	case 63 :
		type="Eu6+3";
		break;
	case 64 :
		type="Gd6+3";
		break;
	case 65 :
		type="Tb6+3";
		break;
	case 66 :
		type="Dy6+3";
		break;
	case 67 :
		type="Ho6+3";
		break;
	case 68 :
		type="Er6+3";
		break;
	case 69 :
		type="Tm6+3";
		break;
	case 70 :
		type="Yb6+3";
		break;
	case 71 :
		type="Lu6+3";
		break;
	case 72 :
		type="Hf3+4";
		break;
	case 73 :
		type="Ta3+5";
		break;
	case 74 :
	{
		int nbBonds=0;
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) nbBonds++;
		if(nbBonds>4) type="W_6+6";
		else {
			int nbOxidation(atom->GetFormalCharge());
			for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
				uint8_t Zlinked=(*itBond)->GetNbrAtom(&(*atom))->GetAtomicNum();
				if(Zlinked!=74) {
					bool WGetElec(true);
					for(uint8_t Z : {6,7,8,9,16,17,34,35,36,53,54,79}) {
						if(Zlinked==Z) WGetElec=false;
					}
					if(WGetElec) nbOxidation-=(*itBond)->GetBondOrder();
					else nbOxidation+=(*itBond)->GetBondOrder();
				}
			}
			if(nbOxidation<5) type="W_3+4";
			else type="W_3+6";
		}
		break;
	}
	case 75 :
	{
		int nbBonds=0;
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) nbBonds++;
		if(nbBonds<5) type="Re3+7";
		else type="Re6+5";
		break;
	}
	case 76 :
		type="Os6+6";
		break;
	case 77 :
		type="Ir6+3";
		break;
	case 78 :
		type="Pt4+2";
		break;
	case 79 :
		type="Au4+3";
		break;
	case 80 :
		type="Hg1+2";
		break;
	case 81 :
		type="Tl3+3";
		break;
	case 82 :
		type="Pb3";
		break;
	case 83 :
		type="Bi3+3";
		break;
	case 84 :
		type="Po3+2";
		break;
	case 85 :
		type="At";
		break;
	case 86 :
		type="Rn4+4";
		break;
	case 87 :
		type="Fr";
		break;
	case 88 :
		type="Ra6+2";
		break;
	case 89 :
		type="Ac6+3";
		break;
	case 90 :
		type="Th6+4";
		break;
	case 91 :
		type="Pa6+4";
		break;
	case 92 :
		type="U_6+4";
		break;
	case 93 :
		type="Np6+4";
		break;
	case 94 :
		type="Pu6+4";
		break;
	case 95 :
		type="Am6+4";
		break;
	case 96 :
		type="Cm6+3";
		break;
	case 97 :
		type="Bk6+3";
		break;
	case 98 :
		type="Cf6+3";
		break;
	case 99 :
		type="Es6+3";
		break;
	case 100 :
		type="Fm6+3";
		break;
	case 101 :
		type="Md6+3";
		break;
	case 102 :
		type="No6+3";
		break;
	case 103 :
		type="Lw6+3";
		break;
	default :
		std::cerr << "Warning : atom type not available for atoms beyond Lawrencium" << std::endl;
		break;
	}
	return type;
}


/// @brief Get the subtype for an atom of a molecule (case COMPASS)
/// @param [in] atom Iterator on the atom
/// @return Subtype
template <> std::string TypeReader<FF::COMPASS>::get(OpenBabel::OBMolAtomIter atom) {
	std::string type;
	// Subtype construction depend on the element
	switch(atom->GetAtomicNum()) {
	case 18 :
		type="ar";
		break;
	case 35 :
		type="br1";
		break;
	case 6 :
	{
		int nbBonds(0), nbTriple(0), nbPolar(0), nbH(0), nbO(0), nbCl(0), nbOm(0);
		bool doubleO(false), doubleN(false), doubleC(false);
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
			OpenBabel::OBAtom* linkedAtom((*itBond)->GetNbrAtom(&(*atom)));
			uint8_t Zlinked=linkedAtom->GetAtomicNum();
			nbBonds++;
			if(Zlinked==1) nbH++;
			if(Zlinked==8) {
				nbO++;
				if(linkedAtom->GetFormalCharge()==-1) nbOm++;
			}
			if(Zlinked==17) nbCl++;
			if((*itBond)->GetBondOrder()==3) nbTriple++;
			if((*itBond)->GetBondOrder()==2) {
				if(Zlinked==6) doubleC=true;
				if(Zlinked==7) doubleN=true;
				if(Zlinked==8) doubleO=true;
			}
			for(uint8_t Z : {7,8,9,16,17,35,53}) {
				if(Zlinked==Z) nbPolar++;
			}
		}
		if(nbBonds==1) type="c1o";
		else if(nbBonds==2) {
			if(nbTriple>0) type="c2t";
			else type="c2=";
		}
		else if(nbBonds==3) {
			if (atom->IsAromatic()) type="c3a";
			else if(doubleO) {
				if(nbO==3) type="c3#";
				else if(nbPolar==3) type="c3''";
				else if(nbOm>0) type="c3-";
				else if(nbPolar==2) type="c3'";
				else type="c3o";
			}
			else if(doubleC) type="c3=";
			else if(doubleN) type="c3n";
			else type="c3";
		}
		else {
			if(nbH==0) type="c44";
			if(nbH==1) type="c43";
			if(nbO>0) type="c4o";
			if(nbCl>0) type="c4x";
			else type="c4";
		}
		break;
	}
	case 9 :
	case 17 :
	{
		uint8_t Z(atom->GetAtomicNum());
		if(Z==9) type="f";
		else type="cl";
		OpenBabel::OBAtom* linkedAtom=(*atom->BeginBonds())->GetNbrAtom(&(*atom));
		if(linkedAtom->GetAtomicNum()==15) {
			int nbBond(0);
			for(OpenBabel::OBBondIterator itBond=linkedAtom->BeginBonds();itBond!=linkedAtom->EndBonds();itBond++) {
				nbBond+=(*itBond)->GetBondOrder();
			}
			if (nbBond==5) type+="1p";
		}
		else if(linkedAtom->GetAtomicNum()==6) {
			int nbHalogen(0);
			for(OpenBabel::OBBondIterator itBond=linkedAtom->BeginBonds();itBond!=linkedAtom->EndBonds();itBond++) {
				uint8_t ZlinkedOfLinked=(*itBond)->GetNbrAtom(&(*linkedAtom))->GetAtomicNum();
				if(ZlinkedOfLinked==9||ZlinkedOfLinked==17||ZlinkedOfLinked==35||ZlinkedOfLinked==53) nbHalogen++;
			}
			if(nbHalogen==4) type+="14";
			else if(nbHalogen==3) type+="13";
			else if(nbHalogen==2) type+="12";
		}
		else type+="1";
		break;
	}
	case 1 :
		if(atom->GetFormalCharge()==1) type="h1+";
		else {
			uint8_t Zlinked=(*atom->BeginBonds())->GetNbrAtom(&(*atom))->GetAtomicNum();
			switch(Zlinked) {
			case 1 :
				type="h1h";
				break;
			case 7 :
			case 17 :
				type="h1n";
				break;
			case 8 :
			case 9 :
				type="h1o";
				break;
			default :
				type="h1";
				break;
			}
		}
		break;
	case 2 :
		type="he";
		break;
	case 53 :
		type="i1";
		break;
	case 36 :
		type="kr";
		break;
	case 7 :
	{
		int nbBond(0), nbH(0), nbN(0), nbO(0), nbTriple(0);
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
			nbBond++;
			if((*itBond)->GetBondOrder()==3) nbTriple++;
			uint8_t Zlinked=(*atom->BeginBonds())->GetNbrAtom(&(*atom))->GetAtomicNum();
			if(Zlinked==1) nbH++;
			else if(Zlinked==7) nbN++;
			else if(Zlinked==8) nbO++;
		}
		switch(nbBond) {
		case 1 :
			if(nbN>0) type="n1n";
			else if(nbO>0) type="n1o";
			else type="n1t";
			break;
		case 2 :
			if(atom->IsAromatic()) type="n2a";
			else if(nbTriple>0) type="n2t";
			else type="n2=";
			break;
		case 3 :
			if(atom->IsAromatic()) type="n3a";
			else if(atom->IsAmideNitrogen()) {
				if(nbH>0) type="n3mh";
				else type="n3m";
			}
			else if(nbO>1) type="n3o";
			else switch(nbH){
			case 3 :
				type="n3*";
				break;
			case 2 :
				type="n3h2";
				break;
			case 1 :
				type="n3h1";
				break;
			default :
				type="n3";
			}
			default :
				if(nbO>0) type="n4o";
				else type="n4+";
		}
		break;
	}
	case 10 :
		type="ne";
		break;
	case 8 :
	{
		int nbBond(0), nbH(0), nbC(0), nbN(0), nbO(0), nbSi(0), nbDouble(0), nbTriple(0);
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
			nbBond++;
			if((*itBond)->GetBondOrder()==2) nbDouble++;
			else if((*itBond)->GetBondOrder()==3) nbTriple++;
			uint8_t Zlinked=(*atom->BeginBonds())->GetNbrAtom(&(*atom))->GetAtomicNum();
			if(Zlinked==1) nbH++;
			else if(Zlinked==6) nbO++;
			else if(Zlinked==7) nbN++;
			else if(Zlinked==8) nbO++;
			else if(Zlinked==14) nbSi++;
		}
		switch(nbBond) {
		case 1 :
			if(nbN>0) {
				if(atom->IsNitroOxygen()) type="o12";
				else type="o1n";
			}
			else if(nbO>0) type="o1o";
			else if(nbTriple>0) type="o1c";
			else if(nbDouble>0) {
				OpenBabel::OBAtom* linkedAtom=(*atom->BeginBonds())->GetNbrAtom(&(*atom));
				int nbBondLinked(0);
				for(OpenBabel::OBBondIterator itBond=linkedAtom->BeginBonds();itBond!=linkedAtom->EndBonds();itBond++) {
					nbBondLinked++;
				}
				if(nbBondLinked>1) type="o1=*";
				else type="o1=";
			}
			else type="o1-";
			break;
		case 2 :
			if(atom->IsAromatic()) type="o2a";
			else if(nbSi>0) type="o2z";
			else if(nbO>0) type="o2b";
			else if(nbH+nbC==2) {
				bool caboxyl(false);
				for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
					OpenBabel::OBAtom* linkedAtom((*itBond)->GetNbrAtom(&(*atom)));
					for(OpenBabel::OBBondIterator itBondOfLinked=linkedAtom->BeginBonds();itBondOfLinked!=linkedAtom->EndBonds();itBondOfLinked++) {
						if((*itBondOfLinked)->GetBondOrder()==2&&(*itBondOfLinked)->GetNbrAtom(linkedAtom)->IsOxygen()) caboxyl=true;
					}
				}
				if(nbH>1) type="o2*";
				else if(caboxyl) {
					if(nbH>0) type="o2c";
					else type="o2s";
				}
				else {
					if(nbH>0) type="o2h";
					else type="o2e";
				}
			}
			else type="o2";
			break;
		default :
			type="o3";
			break;
		}
		break;
	}
	case 15 :
		type="p4=";
		break;
	case 16 :
		type="s2";
		break;
	case 14 :
	{
		int nbH(0),nbO(0);
		for(OpenBabel::OBBondIterator itBond=atom->BeginBonds();itBond!=atom->EndBonds();itBond++) {
			uint8_t Zlinked=(*itBond)->GetNbrAtom(&(*atom))->GetAtomicNum();
			if(Zlinked==1) nbH++;
			if(Zlinked==8) nbO++;
		}
		if(nbO>0&&nbH==0) type="si4c";
		else type="si4";
		break;
	}
	case 54 :
		type="xe";
		break;
	default :
		std::cerr << "Warning : no atom type for Z=" << atom->GetAtomicNum() <<"in Compass" << std::endl;
		break;
	}
	return type;
}


/// @brief Get the partial charge for an atom (case NONE)
/// @param [in,out] toComplete MPI atom to complete
/// @param [in] atom Iterator on the OpenBabel atom
template <> void ChargeReader<MPI__Atom, FF::NONE>::get(MPI__PInMol<MPI__Atom>* toComplete, OpenBabel::OBMolAtomIter atom){}


/// @brief Get the partial charge for an atom (case OPENBABEL)
/// @param [in,out] toComplete MPI atom to complete
/// @param [in] atom Iterator on the OpenBabel atom
template <> void ChargeReader<MPI__AtomCharged, FF::OPENBABEL>::get(MPI__PInMol<MPI__AtomCharged>* toComplete, OpenBabel::OBMolAtomIter atom){
	toComplete->charge()=atom->GetPartialCharge();
}


/// @brief Get the partial charge for an atom (case FROPMCOMPASS)
/// @param [in,out] toComplete MPI atom to complete
/// @param [in] atom Iterator on the OpenBabel atom
template <>void ChargeReader<MPI__AtomCharged, FF::FROMCOMPASS>::get(MPI__PInMol<MPI__AtomCharged>* toComplete, OpenBabel::OBMolAtomIter atom){
	std::cerr << "Charge method FROMCOMPASS not implemented yet" << std::endl;
}

#endif
