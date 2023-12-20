/*
 * forceFieldExpert.hpp
 *
 *  Created on: Sep 19, 2017
 *      Author: giarda
 */

/// @file
/// @brief Class for the expert forcefields

#ifndef FORCEFIELDEXPERT_HPP_
#define FORCEFIELDEXPERT_HPP_


#include <map>

#include "forceField/forceField.hpp"
#include "forceField/formComputer.hpp"


/// @brief Class used to calculate the intramolecular forces with expert mode forcefield
class ForceFieldExpert: public ForceField {

public:

	ForceFieldExpert();
	virtual ~ForceFieldExpert();
	/// @brief Compute the forces and energy on bonds
	/// @param [in] size Number of bonds
	/// @param [out] e Energies
	/// @param [out] fx X components of the forces
	/// @param [out] fy Y components of the forces
	/// @param [out] fz Z components of the forces
	/// @param [in] rx X component of the distance
	/// @param [in] ry Y component of the distance
	/// @param [in] rz Z component of the distance
	/// @param [in] params Pointer to the parameters
	/// @param [in] form Functional form for the force calculation (unused in non expert force fields)
	virtual void computeForceBonds(const uint size, double* e, double* fx, double* fy, double* fz,
			const double* rx, const double* ry, const double* rz, double** params, FF::FunctForm form) const {
		m_bVecto.at(form)->compute(size,e,fx,fy,fz,rx,ry,rz,params);
	}
	/// @brief Compute the forces and energy on angles (123)
	///
	/// This calculate only the forces for the side atoms as the forces on the
	/// middle atom can be easily deducted by considering that total force is
	/// zero
	/// @param [in] size Number of angles
	/// @param [out] e Energies
	/// @param [out] fx1 X components of the forces on 1
	/// @param [out] fy1 Y components of the forces on 1
	/// @param [out] fz1 Z components of the forces on 1
	/// @param [out] fx2 X components of the forces on 3
	/// @param [out] fy2 Y components of the forces on 3
	/// @param [out] fz2 Z components of the forces on 3
	/// @param [in] rx21 X component of the vector between 2 and 1
	/// @param [in] ry21 Y component of the vector between 2 and 1
	/// @param [in] rz21 Z component of the vector between 2 and 1
	/// @param [in] rx23 X component of the vector between 2 and 3
	/// @param [in] ry23 Y component of the vector between 2 and 3
	/// @param [in] rz23 Z component of the vector between 2 and 3
	/// @param [in] params Pointer to the parameters
	/// @param [in] form Functional form for the force calculation
	virtual void computeForceAngles(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, double** params, FF::FunctForm form) const {
		m_aVecto.at(form)->compute(size,e,fx1,fy1,fz1,fx2,fy2,fz2,rx21,ry21,rz21,rx23,ry23,rz23,params);
	}
	/// @brief Compute the forces and energy on dihedral angles (1234)
	///
	/// This calculate only the forces for the atoms 1,2 and 4 as the forces on the
	/// atom 4 can easily be deducted by considering that total force is zero
	/// @param [in] size Number of dihedral angles
	/// @param [out] e Energies
	/// @param [out] fx1 X components of the forces on 1
	/// @param [out] fy1 Y components of the forces on 1
	/// @param [out] fz1 Z components of the forces on 1
	/// @param [out] fx2 X components of the forces on 2
	/// @param [out] fy2 Y components of the forces on 2
	/// @param [out] fz2 Z components of the forces on 2
	/// @param [out] fx4 X components of the forces on 4
	/// @param [out] fy4 Y components of the forces on 4
	/// @param [out] fz4 Z components of the forces on 4
	/// @param [in] rx21 X component of the vector between 2 and 1
	/// @param [in] ry21 Y component of the vector between 2 and 1
	/// @param [in] rz21 Z component of the vector between 2 and 1
	/// @param [in] rx23 X component of the vector between 2 and 3
	/// @param [in] ry23 Y component of the vector between 2 and 3
	/// @param [in] rz23 Z component of the vector between 2 and 3
	/// @param [in] rx34 X component of the vector between 3 and 4
	/// @param [in] ry34 Y component of the vector between 3 and 4
	/// @param [in] rz34 Z component of the vector between 3 and 4
	/// @param [in] params Pointer to the parameters
	/// @param [in] form Functional form for the force calculation
	virtual void computeForceDihedrals(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, const double* rx34, const double* ry34, const double* rz34,
			double** params, FF::FunctForm form) const {
		m_dVecto.at(form)->compute(size,e,fx1,fy1,fz1,fx2,fy2,fz2,fx4,fy4,fz4,rx21,ry21,rz21,rx23,ry23,rz23,rx34,ry34,rz34,params);
	}
	/// @brief Compute the forces and energy on improper torsion angles (1-234)
	///
	/// This calculate only the forces for the atoms 2,3 and 4 as the forces on the
	/// atom 3 can easily be deducted by considering that total force is zero
	/// @param [in] size Number of dihedral angles
	/// @param [out] e Energies
	/// @param [out] fx2 X components of the forces on 2
	/// @param [out] fy2 Y components of the forces on 2
	/// @param [out] fz2 Z components of the forces on 2
	/// @param [out] fx3 X components of the forces on 3
	/// @param [out] fy3 Y components of the forces on 3
	/// @param [out] fz3 Z components of the forces on 3
	/// @param [out] fx4 X components of the forces on 4
	/// @param [out] fy4 Y components of the forces on 4
	/// @param [out] fz4 Z components of the forces on 4
	/// @param [in] rx21 X component of the vector between 2 and 1
	/// @param [in] ry21 Y component of the vector between 2 and 1
	/// @param [in] rz21 Z component of the vector between 2 and 1
	/// @param [in] rx31 X component of the vector between 3 and 1
	/// @param [in] ry31 Y component of the vector between 3 and 1
	/// @param [in] rz31 Z component of the vector between 3 and 1
	/// @param [in] rx41 X component of the vector between 4 and 1
	/// @param [in] ry41 Y component of the vector between 4 and 1
	/// @param [in] rz41 Z component of the vector between 4 and 1
	/// @param [in] params Pointer to the parameters
	/// @param [in] form Functional form for the force calculation
	virtual void computeForceImpropers(const uint size, double* e, double* fx2, double* fy2, double* fz2, double* fx3, double* fy3, double* fz3, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx31, const double* ry31, const double* rz31, const double* rx41, const double* ry41, const double* rz41,
			double** params, FF::FunctForm form) const {
		m_iVecto.at(form)->compute(size,e,fx2,fy2,fz2,fx3,fy3,fz3,fx4,fy4,fz4,rx21,ry21,rz21,rx31,ry31,rz31,rx41,ry41,rz41,params);
	}
	/// @brief Create the parameters for a type of dof
	/// @param [in] kindDof Kind of dof
	/// @param [in] typeDof Type of the dof
	virtual void setParam(const uint8_t kindDof, const std::string typeDof);
	/// @brief Get the parameters for a dof in a vectorization buffer
	/// @param [in] typeDof Type of the dof
	/// @param [out] params Parameters
	virtual void getParam(const std::string typeDof, std::vector<double>& params) const;
	/// @brief Begin the initialization of the force field
	virtual void beginInit();
	/// @brief End the initialization of the force field
	virtual void endInit();
	/// @brief Get the type of an atom
	/// @param [in] subtype Subtype of said atom
	inline virtual std::string getAtomType(std::string subtype) const {return (*m_initOfTypes)[subtype];}
	/// @brief Create the potential between to given atom types
	/// @param [in] typei First atom type
	/// @param [in] typej Second atom type
	/// @return Created potential
	virtual IPotential* makePotential(std::string typei, std::string typej) const;

private :

	std::unordered_map<std::string,ForceFieldParam> m_parameters; ///< Parameters for this forcefield
//	std::unordered_map<std::string,std::vector<double>>* m_initOfBondedParams; ///< Pointer to a temporary array to store the data used to build the parameters for bonded interactions
//	std::unordered_map<std::string,std::vector<double>>* m_initOfNonBondedParams; ///< Pointer to a temporary array to store the data used to build the parameters for non bonded interactions
	std::unordered_map<std::string,std::string>* m_initOfTypes; ///< Pointer to a temporary array to store the type/subtype correspondence
//	std::ifstream m_initFile; ///< Parameters data file

	std::map<FF::FunctForm,BondComputer*> m_bVecto; ///< Vectorized bond forces calculations
	std::map<FF::FunctForm,AngleComputer*> m_aVecto; ///< Vectorized angle forces calculations
	std::map<FF::FunctForm,DihedralComputer*> m_dVecto; ///< Vectorized dihedral forces calculations
	std::map<FF::FunctForm,ImproperComputer*> m_iVecto; ///< Vectorized improper forces calculations

};

#endif /* FORCEFIELDEXPERT_HPP_ */
