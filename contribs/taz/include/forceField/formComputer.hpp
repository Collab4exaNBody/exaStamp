/*
 * formComputer.hpp
 *
 *  Created on: Oct 17, 2017
 *      Author: giarda
 */

/// @file
/// @brief Tools for the expert force field

#ifndef FORMCOMPUTER_HPP_
#define FORMCOMPUTER_HPP_

#include <unordered_map>

#include "simd/forceField/bharm.hpp"
#include "simd/forceField/acosto3.hpp"
#include "simd/forceField/dcosn.hpp"
#include "simd/forceField/icosto2.hpp"


/// @brief Generic class used to calculate the bond forces with expert mode forcefield
class BondComputer {

public :

	/// @brief Default constructor
	BondComputer() {}
	/// @brief Destructor (nothing to do)
	~BondComputer() {}

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
	virtual void compute(const uint size, double* e, double* fx, double* fy, double* fz,
			const double* rx, const double* ry, const double* rz, double** params) const =0;

private :

};


/// @brief Class used to calculate the bond forces with the bharm functional form
class BHarmComputer : BondComputer {

public :

	BHarmComputer() {}
	/// @brief Destructor (nothing to do)
	~BHarmComputer() {}

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
	virtual void compute(const uint size, double* e, double* fx, double* fy, double* fz,
			const double* rx, const double* ry, const double* rz, double** params) const {
		m_compute(fx, fy, fz, e, rx, ry, rz, params[0], params[1], size);
	}

private :

	simd::kernels::bharm<double> m_compute; ///< Vectorized forces calculation

};


/// @brief Generic class used to calculate the angle forces with expert mode forcefield
class AngleComputer {

public :

	/// @brief Default constructor
	AngleComputer() {}
	/// @brief Destructor (nothing to do)
	~AngleComputer() {}

	/// @brief Compute the forces and energy on angles
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
	virtual void compute(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, double** params) const =0;

private :

};


/// @brief Class used to calculate the angle forces with the acosto3 functional form
class ACosTo3Computer : AngleComputer {

public :

	ACosTo3Computer() {}
	/// @brief Destructor (nothing to do)
	~ACosTo3Computer() {}

	/// @brief Compute the forces and energy on angles
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
	virtual void compute(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2,
				const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, double** params) const {
		m_compute(fx1, fy1, fz1, fx2, fy2, fz2, e, rx21, ry21, rz21, rx23, ry23, rz23, params[0], params[1], params[2], params[3], size);
	}

private :

	simd::kernels::acosto3<double> m_compute; ///< Vectorized forces calculation

};


/// @brief Generic class used to calculate the dihedral forces with expert mode forcefield
class DihedralComputer {

public :

	DihedralComputer() {}
	/// @brief Destructor (nothing to do)
	~DihedralComputer() {}

	/// @brief Compute the forces and energy on dihedrals
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
	virtual void compute(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, const double* rx34, const double* ry34, const double* rz34,
			double** params) const =0;

private :

};


/// @brief Class used to calculate the dihedral forces with the dcosn functional form
class DCosNComputer : DihedralComputer {

public :

	DCosNComputer() {}
	/// @brief Destructor (nothing to do)
	~DCosNComputer() {}

	/// @brief Compute the forces and energy on dihedrals
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
	virtual void compute(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, const double* rx34, const double* ry34, const double* rz34,
			double** params) const {
		m_compute(fx1, fy1, fz1, fx2, fy2, fz2, fx4, fy4, fz4, e, rx21, ry21, rz21, rx23, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
	}

private :

	simd::kernels::dcosn<double> m_compute; ///< Vectorized forces calculation

};


/// @brief Generic class used to calculate the improper forces with expert mode forcefield
class ImproperComputer {

public :

	ImproperComputer() {}
	/// @brief Destructor (nothing to do)
	~ImproperComputer() {}

	/// @brief Compute the forces and energy on impropers
	/// @param [in] size Number of improper torsions
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
	virtual void compute(const uint size, double* e, double* fx2, double* fy2, double* fz2, double* fx3, double* fy3, double* fz3, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx31, const double* ry31, const double* rz31, const double* rx41, const double* ry41, const double* rz41,
			double** params) const =0;

private :

};


/// @brief Class used to calculate the improper forces with the icosto2 functional form
class ICosTo2Computer : ImproperComputer {

public :

	ICosTo2Computer() {}
	/// @brief Destructor (nothing to do)
	~ICosTo2Computer() {}

	/// @brief Compute the forces and energy on impropers
	/// @param [in] size Number of improper torsions
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
	virtual void compute(const uint size, double* e, double* fx2, double* fy2, double* fz2, double* fx3, double* fy3, double* fz3, double* fx4, double* fy4, double* fz4,
				const double* rx21, const double* ry21, const double* rz21, const double* rx31, const double* ry31, const double* rz31, const double* rx41, const double* ry41, const double* rz41,
				double** params) const {
		m_compute(fx2, fy2, fz2, fx3, fy3, fz3, fx4, fy4, fz4, e, rx21, ry21, rz21, rx31, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
	}

private :

	simd::kernels::icosto2<double> m_compute; ///< Vectorized forces calculation

};

#endif /* FORMCOMPUTER_HPP_ */
