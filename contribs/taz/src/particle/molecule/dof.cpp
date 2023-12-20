/*
 * dof.cpp
 *
 *  Created on: May 31, 2016
 *      Author: giarda
 */

/// @file
/// @brief Implementations for the class Dof


#include "particle/molecule/dof.hpp"


/// @brief Constructor
/// @param [in] atom1 Index of the first atom
/// @param [in] atom2 Index of the second atom
/// @param [in] atom3 Index of the third atom, default=-1 (none)
/// @param [in] atom4 Index of the forth atom, default=-1 (none)
/// @param [in] improper Indicate if the torsion is an improper torsion, default=false
Dof::Dof(uint64_t atom1, uint64_t atom2, uint64_t atom3, uint64_t atom4, bool improper) {
	// Set the atoms
	m_atoms.push_back(atom1);
	m_atoms.push_back(atom2);
	for(int i : {atom3,atom4}) {
		if(i>0) m_atoms.push_back(i);
	}
	// Set the type
	m_kind=m_atoms.size()-2;
	if(improper) m_kind=3;
	m_type="";
}


/// @brief Stream insertion operator for a degree of freedom
/// @param [in,out] out Stream
/// @param [in] dof Degree of freedom
std::ostream& operator << (std::ostream& out, const Dof dof) {
	switch(dof.m_kind){
	case 0 :
		out << "Bond of type " << dof.m_type << " : ";
		break;
	case 1 :
		out << "Angle of type " << dof.m_type << " : ";
		break;
	case 2 :
		out << "Torsion of type " << dof.m_type << " : ";
		break;
	case 3 :
		out << "Improper torsion of type " << dof.m_type << " : ";
		break;
	}
	for (uint64_t atom : dof.m_atoms){
		out << atom << " ";
	}
  return out;
}
