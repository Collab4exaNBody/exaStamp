/*
 * periodicTable.h
 *
 *  Created on: Feb 19, 2016
 *      Author: giarda
 */

/// @file
/// @brief Class to store and use the periodic table

#ifndef PERIODICTABLE_H_
#define PERIODICTABLE_H_


#include <vector>
#include <string>


/// @brief A class for the periodic table
class PeriodicTable {

public :

	/// @brief Default constructor
	PeriodicTable() :
		m_symbols({"H","He",
			"Li","Be","B","C","N","O","F","Ne",
			"Na","Mg","Al","Si","P","S","Cl","Ar",
			"K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
			"Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
			"Cs","Ba",
			"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
			"Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
			"Fr","Ra",
			"Ac","Th","Pa","U","Np","Pu","Am","Cu","Bk","Cf","Es","Fm","Md","No","Lr",
			"Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"})
	{}
	/// @brief Destructor
	///
	///
	~PeriodicTable() {
		m_symbols.clear();
		std::vector<std::string>().swap(m_symbols);
	}
	/// @brief Get an element of the table
	/// @param [in] atNum Atomic number
	/// @return Symbol of the element
	std::string operator[](int atNum) {return m_symbols[atNum-1];}

private :

	std::vector<std::string> m_symbols; ///< Symbol for each element of the periodic table (indexed by atomic number)

};

#endif /* PERIODICTABLE_H_ */
