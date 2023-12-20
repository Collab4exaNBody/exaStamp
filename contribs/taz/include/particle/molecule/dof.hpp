/*
 * dof.hpp
 *
 *  Created on: May 31, 2016
 *      Author: giarda
 */

/// @file
/// @brief Definition of the class Dof

#ifndef DOF_HPP_
#define DOF_HPP_


#include <string>
#include <vector>
#include <iostream>

#include "utils/stringUtils.hpp"

#include "forceField/forceField.hpp"


/// @brief Class to handle a degree of freedom
class Dof {

public :

	/// @brief Default constructor
	Dof() {}
	Dof(uint64_t atom1, uint64_t atom2, uint64_t atom3=0, uint64_t atom4=0, bool improper=false);
	/// @brief Destructor (nothing to do)
	~Dof() {}
	/// @brief Constant accessor to the atoms
	const std::vector<uint64_t> atoms() {return m_atoms;}
	/// @brief Constant accessor to one atom of the dof
	/// @param [in] i Atom to get
	/// @return Index of the atom
	uint64_t index(int i) const {return m_atoms[i];}
	/// @brief Accessor to the type
	std::string& type() {return m_type;}
	/// @brief Constant accessor to the type
	std::string type() const {return m_type;}
	/// @brief Constant accessor to the kind
	uint8_t kind() const {return m_kind;}
	/// @brief Comparison operator
	/// @param [in] otherDof Dof to compare to
	bool operator<(Dof const& otherDof) const {
		return m_kind<otherDof.m_kind;
	}
	template<typename Lambda1, typename Lambda2, typename Lambda3>
	void setType(const FF::Type ff, Lambda1 getSubtype, Lambda2 getBondOrder, Lambda3 getBonds);

private :

	uint8_t m_kind;						///< Kind of degree of freedom (0 for bond, 1 for angle, 2 for torsion and 3 for improper)
	std::vector<uint64_t> m_atoms;	///< Atoms of the degree of freedom
	std::string m_type;				///< Type of the degree of freedom (depend on the subtypes of its atoms)


	template<typename Lambda1, typename Lambda2, typename Lambda3>
	void setTypeUFF(Lambda1 getSubtype, Lambda2 getBondOrder, Lambda3 getBonds);

  friend std::ostream& operator << (std::ostream& out, const Dof dof);

};


/// @brief Set the type of the dof
///
/// Switch through the possible force field and call the appropriate specialized function
/// @tparam Lambda1 Type of the function used to get the subtype of an atom
/// @tparam Lambda2 Type of the function used to get the bond order between two atoms
/// @tparam Lambda3 Type of the function used to get the bonds of an atom
/// @param [in] ff Force field
/// @param [in] getSubtype Function used to get the subtype of an atom
/// @param [in] getBondOrder Function used to get the bond order between two atoms
/// @param [in] getBonds Function used to get the bonds of an atom
template<typename Lambda1, typename Lambda2, typename Lambda3>
void Dof::setType(const FF::Type ff, Lambda1 getSubtype, Lambda2 getBondOrder, Lambda3 getBonds) {
	switch(ff) {
	case FF::UFF :
		setTypeUFF<Lambda1, Lambda2, Lambda3>(getSubtype, getBondOrder, getBonds);
		break;
	default :
		break;
	}
}


/// @brief Set the type of the dof (case of the UFF force field)
/// @tparam Lambda1 Type of the function used to get the subtype of an atom
/// @tparam Lambda2 Type of the function used to get the bond order between two atoms
/// @tparam Lambda3 Type of the function used to get the bonds of an atom
/// @param [in] getSubtype Function used to get the subtype of an atom
/// @param [in] getBondOrder Function used to get the bond order between two atoms
/// @param [in] getBonds Function used to get the bonds of an atom
template<typename Lambda1, typename Lambda2, typename Lambda3>
void Dof::setTypeUFF(Lambda1 getSubtype, Lambda2 getBondOrder, Lambda3 getBonds) {
	// Initialize the dof type
	m_type="";
	// Lambda function to calculate bond order between to atoms
	// Parameters : atom1 -> index of the first atom of the bond
	// 							atom2 -> index of the second atom of the bond
	auto order = [&] (uint64_t atom1, uint64_t atom2) -> double {
		// The nonintegral cases
		if(getBondOrder(atom1,atom2)==-1) {
			bool isAmide(false);
			// Get the elements of the two atoms
			std::string element1=getSubtype(atom1).substr(0,2);
			std::string element2=getSubtype(atom2).substr(0,2);
			// If 1 is carbon and 2 is nitrogen
			if(is_equal(element1,"C_")&&is_equal(element2,"N_")) {
				std::vector<uint64_t> bondsOf1=getBonds(atom1);
				// And one of the atoms bonded to 1
				for(uint64_t neighbor : bondsOf1) {
					std::string elementN=getSubtype(neighbor).substr(0,2);
					// Is an oxygen and has a double bond to 1, then its an amide bond
					if(is_equal(elementN,"O_")&&getBondOrder(atom1,neighbor)==2) isAmide=true;
				} // End loop bonds
			} // End case carbon/nitrogen
			// If 2 is carbon and 1 is nitrogen
			else if(is_equal(element2,"C_")&&is_equal(element1,"N_")) {
				std::vector<uint64_t> bondsOf2=getBonds(atom2);
				// And one of the atoms bonded to 2
				for(uint64_t neighbor : bondsOf2) {
					std::string elementN=getSubtype(neighbor).substr(0,2);
					// Is an oxygen and has a double bond to 2, then its an amide bond
					if(is_equal(elementN,"O_")&&getBondOrder(atom2,neighbor)==2) isAmide=true;
				} // End loop bonds
			} // End case nitrogen/carbon
			// Amide bonds have bond order 1,41
			if(isAmide) return 1.41;
			// Else it's an aromatic bond and bond order is 1,5
			else return 1.5;
		} // End of nonintegral cases
		// Simple case : convert integer reminder to double
		else return double(getBondOrder(atom1,atom2));
	};
	// Type construction depend on the kind of the dof
	switch(m_kind) {
	// Case bond
	case 0 :
	{
		// Get the two atoms subtypes and bond order
		std::string subtype1=getSubtype(m_atoms[0]);
		std::string subtype2=getSubtype(m_atoms[1]);
		double bondOrder=order(m_atoms[0],m_atoms[1]);
		// Type is typeOfSmallerTypeAtom*floatingPointBondOrder*typeOfBiggerTypeAtom
		if(subtype1.compare(subtype2)<0) m_type=subtype1+"*"+std::to_string(bondOrder).substr(0,4)+"*"+subtype2;
		else m_type=subtype2+"*"+std::to_string(bondOrder).substr(0,4)+"*"+subtype1;
		break;
	} // End of bond case
	// Case angle
	case 1 :
	{
		// Get the atoms subtypes and bond orders
		std::string subtype1=getSubtype(m_atoms[0]);
		std::string subtype2=getSubtype(m_atoms[1]);
		std::string subtype3=getSubtype(m_atoms[2]);
		double bondOrder1=order(m_atoms[0],m_atoms[1]);
		double bondOrder2=order(m_atoms[1],m_atoms[2]);
		// Compare the subtypes and exchange atom 1 and 3 if necessary
		if(subtype1.compare(subtype3)<0) {
			std::string subtypeTmp=subtype1;
			subtype1=subtype3;
			subtype3=subtypeTmp;
			double bondOrderTmp(bondOrder1);
			bondOrder1=bondOrder2;
			bondOrder2=bondOrderTmp;
		}
		else if(subtype1.compare(subtype3)==0) {
			if(bondOrder2<bondOrder1) {
				std::string subtypeTmp=subtype1;
				subtype1=subtype3;
				subtype3=subtypeTmp;
				double bondOrderTmp(bondOrder1);
				bondOrder1=bondOrder2;
				bondOrder2=bondOrderTmp;
			}
		}
		// Type is typeOfSmallerTypeBorderAtom*floatingPointBondOrder*typeOfCentralAtom*floatingPointBondOrder*typeOfBiggerTypeBorderAtom
		m_type=subtype1+"*"+std::to_string(bondOrder1).substr(0,4)+"*"+subtype2+"*"+std::to_string(bondOrder2).substr(0,4)+"*"+subtype3;
		break;
	} // End of angle case
	// Case torsion
	case 2 :
	{
		// Get the two central atoms subtypes
		std::string subtype1=getSubtype(m_atoms[1]);
		std::string subtype2=getSubtype(m_atoms[2]);
		// The torsion is considered only if both atoms are sp2 or sp3
		if((subtype1[2]=='2'||subtype1[2]=='R'||subtype1[2]=='3')&&(subtype2[2]=='2'||subtype2[2]=='R'||subtype2[2]=='3')) {
			// Get the bond order
			double bondOrder=order(m_atoms[1], m_atoms[2]);
			// Type is typeOfSmallerTypeAtom*floatingPointBondOrder*typeOfBiggerTypeAtom
			if(subtype1.compare(subtype2)<0) m_type="t"+subtype1+"*"+std::to_string(bondOrder).substr(0,4)+"*"+subtype2;
			else m_type="t"+subtype2+"*"+std::to_string(bondOrder).substr(0,4)+"*"+subtype1;
			// If one atom is sp3 and the other is sp2 bonded to another sp2 add *pc to the type
			bool propCase(false);
			if((subtype1[2]=='2'||subtype1[2]=='R')&&subtype2[2]=='3') {
				std::vector<uint64_t> bondsOf1=getBonds(m_atoms[1]);
				for(uint64_t neighbor : bondsOf1) {
					if(getSubtype(neighbor)[2]=='2'||getSubtype(neighbor)[2]=='R') propCase=true;
				}
			}
			if((subtype2[2]=='2'||subtype2[2]=='R')&&subtype1[2]=='3') {
				std::vector<uint64_t> bondsOf2=getBonds(m_atoms[2]);
				for(uint64_t neighbor : bondsOf2) {
					if(getSubtype(neighbor)[2]=='2'||getSubtype(neighbor)[2]=='R') propCase=true;
				}
			}
			if(propCase) m_type+="*pc";
		} // End of the considered cases
		break;
	} // End of torsion case
	// Case improper torsion
	case 3 :
	{
		// Get the central atom subtype and bonds
		std::string subtype=getSubtype(m_atoms[0]);
		std::vector<uint64_t> bonds=getBonds(m_atoms[0]);
		// The improper torsion is considered only if the central atom is bonded to exactly 3 atoms
		if(bonds.size()==3) {
			// Get the element
			std::string element=subtype.substr(0,2);
			// Type is the element for the elements considered (C, N, O, P, S, As, Se, Sb, Te, Bi, Po)
			if(is_equal(element,"C_")) {
				// If its a carbon
				bool linkedToO2(false);
				// Bonded to
				for(uint64_t neighbor : bonds) {
					std::string elementN=getSubtype(neighbor).substr(0,3);
					// An sp2 oxygen
					if(is_equal(elementN,"O_2")) linkedToO2=true;
				}
				// Add *O to the type
				if(linkedToO2) element+="*O";
				m_type=element;
			} // End of carbon
			else {
				bool group5(is_equal(element,"N_")||is_equal(element,"P_")||is_equal(element,"As")||is_equal(element,"Sb")||is_equal(element,"Bi"));
				bool group6(is_equal(element,"O_")||is_equal(element,"S_")||is_equal(element,"Se")||is_equal(element,"Te")||is_equal(element,"Po"));
				if(group5||group6) {
					m_type=element;
				}
			} // End of group 5 and 6
		} // End of the considered cases
		break;
	} // End of improper torsion case
	} // End of construction
}

#endif /* DOF_HPP_ */
