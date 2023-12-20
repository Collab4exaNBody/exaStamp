/*
 * forceField.hpp
 *
 *  Created on: Feb 26, 2016
 *      Author: giarda
 */

/// @file
/// @brief Class for the forcefields and associated enumerations

#ifndef FORCEFIELD_HPP_
#define FORCEFIELD_HPP_


#include <map>
#include <vector>
#include <cstdint>
#include <string>

using uint = unsigned int;

class IPotential;


/// @brief Namespace for the force field related properties
namespace FF {


	/// @brief Enumeration of the possible force fields
	enum Type {
		UFF, ///< General force field with parameters for the full periodic table, open source, no Coulomb interaction
		COMPASS, ///< Condensed phase optimized force field, not fully available, simple charges
		AMBER, ///< Force field widely used for proteins and DNA, open source
		EXPERT ///< Custom force field
	};


	/// @brief Enumeration of the method to assign the partial charges
	enum ChargeMethod {
//		NOT_SPECIFIED, ///< Use the default charge method for the force field
		NONE, ///< Do not use any charges (all charges to 0.)
		OPENBABEL, ///< Use OpenBabel default charges
		FROMCOMPASS ///< Build Compass charges
	};

	/// @brief Enumeration of the available functional forms to compute intramolecular forces
	enum FunctForm {
		BHARM, ///< Harmonic for bonds
		ACOSTO3, ///< For angles, E=C0+C1*costheta+C2*cos2theta+C3*cos3theta
		DCOSN, ///< For dihedrals, E=C0+C1*cos(n*phi), n=1,2,3 or 6
		ICOSTO2 ///< For impropers, E=C0+C1*cospsi+C2*cos2psi
	};

}  // namespace ForceField


/// @brief Base class the store the parameters of a force field
class ForceFieldParam {

public :

	/// @brief Constructor
	///
	///
	ForceFieldParam() :
		m_params()
	{}
	/// @brief Destructor (nothing to do)
	virtual ~ForceFieldParam() {}

	/// @brief Add a parameter
	/// @param [in] param Parameter to add
	void addParam(double param) {
		m_params.push_back(param);
	}

	/// @brief Get a parameter
	/// @param [in] i Index of the parameter
	/// @return Parameter
	double get(int i) const {
		return m_params[i];
	}

	/// @brief Get all parameters
	/// @param [out] outD Vector of the parameters
	void getAll(std::vector<double>& outD) const {
			outD=m_params;
	}


private :
	std::vector<double> m_params; ///< Parameters

};


/// @brief Class used to calculate the intramolecular forces
class ForceField {

public :

	static double s_rcut; ///< Cutoff radius for the intra-molecular interactions
	static const std::map<FF::Type,std::tuple<uint,uint,uint,uint>> paramBuild; ///< Number of parameters for each possible force field (by dofs)
	static const std::map<FF::FunctForm,int> numParam; ///< Number of parameters for each possible fonctionnal form
	static const std::map<FF::Type,FF::ChargeMethod> defaultCharge; ///< Default charge method for each possible force field

	/// @brief Constructor
	/// @param [in] type Force field type
	/// @param [in] deleted Minimum number of bonds between two atoms to consider a non-bonded interaction
	ForceField(FF::Type type, uint deleted) :
	m_type(type),
	m_deletedNeighbors(deleted)
	{}
	/// @brief Destructor (nothing to do)
	virtual ~ForceField() {}
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
			const double* rx, const double* ry, const double* rz, double** params, FF::FunctForm form) const =0;
	/// @brief Compute the forces and energy on angles (123)
	/// @param [in] size Number of bonds
	/// @param [out] e Energies
	/// @param [out] fx1 X components of the forces on 1
	/// @param [out] fy1 Y components of the forces on 1
	/// @param [out] fz1 Z components of the forces on 1
	/// @param [out] fx2 X components of the forces on 2
	/// @param [out] fy2 Y components of the forces on 2
	/// @param [out] fz2 Z components of the forces on 2
	/// @param [in] rx21 X component of the vector between 2 and 1
	/// @param [in] ry21 Y component of the vector between 2 and 1
	/// @param [in] rz21 Z component of the vector between 2 and 1
	/// @param [in] rx23 X component of the vector between 2 and 3
	/// @param [in] ry23 Y component of the vector between 2 and 3
	/// @param [in] rz23 Z component of the vector between 2 and 3
	/// @param [in] params Pointer to the parameters
	/// @param [in] form Functional form for the force calculation (unused in non expert force fields)
	virtual void computeForceAngles(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, double** params, FF::FunctForm form) const =0;
	/// @brief Compute the forces and energy on dihedral angles (1234)
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
	/// @param [in] form Functional form for the force calculation (unused in non expert force fields)
	virtual void computeForceDihedrals(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, const double* rx34, const double* ry34, const double* rz34,
			double** params, FF::FunctForm form) const =0;
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
	/// @param [in] form Functional form for the force calculation (unused in non expert force fields)
	virtual void computeForceImpropers(const uint size, double* e, double* fx2, double* fy2, double* fz2, double* fx3, double* fy3, double* fz3, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx31, const double* ry31, const double* rz31, const double* rx41, const double* ry41, const double* rz41,
			double** params, FF::FunctForm form) const =0;
	/// @brief Create the parameters for a type of dof
	/// @param [in] kindDof Kind of dof
	/// @param [in] typeDof Type of the dof
	virtual void setParam(const uint8_t kindDof, const std::string typeDof)=0;
	/// @brief Get the parameters for a dof in a vectorization buffer
	/// @param [in] typeDof Type of the dof
	/// @param [out] params Parameters
	virtual void getParam(const std::string typeDof, std::vector<double>& params) const =0;
	/// @brief Create the potential between to given atom types
	/// @param [in] type1 First atom type
	/// @param [in] type2 Second atom type
	/// @return Created potential
	virtual IPotential* makePotential(std::string type1, std::string type2) const = 0;
	/// @brief Begin the initialization of the force field
	virtual void beginInit()=0;
	/// @brief End the initialization of the force field
	virtual void endInit()=0;
	/// @brief Accessor to the type of the force field
	inline FF::Type type() const {return m_type;}
	/// @brief Accessor to the minimum number of bonds between non-bonded interacting atoms
	inline uint deletedNeighbors() const {return m_deletedNeighbors;}
	/// @brief Accessor to the forcefield costs
	/// @param [in] type Type of the dof whose cost is accessed,
	/// 0 for the bonds, 1 for the angles and 2 for the dihedral angles
	inline double cost(const uint type) const {return m_costs[type];}
	/// @brief Get the type of an atom
	/// @param [in] subtype Subtype of said atom
	virtual std::string getAtomType(std::string subtype) const = 0;

protected :

	FF::Type m_type; ///< Type of the force field
	uint m_deletedNeighbors; ///< Minimum number of bonds between two atoms to consider a non-bonded interaction
	std::vector<double> m_costs; ///< Cost of the forcefield for the bonds, angles and dihedral angles

};


#endif /* FORCEFIELD_HPP_ */
