/*
 * forceFieldUFF.hpp
 *
 *  Created on: Mar 21, 2016
 *      Author: giarda
 */

/// @file
/// @brief Class for the UFF forcefield

#ifndef FORCEFIELDUFF_HPP_
#define FORCEFIELDUFF_HPP_


#include <fstream>
#include <unordered_map>

#include "forceField/forceField.hpp"

#include "simd/forceField/bharm.hpp"
#include "simd/forceField/acosto3.hpp"
#include "simd/forceField/dcosn.hpp"
#include "simd/forceField/icosto2.hpp"


/// @brief Class used to calculate the intramolecular forces when UFF is used
class ForceFieldUFF: public ForceField {

public:

	ForceFieldUFF();
	virtual ~ForceFieldUFF();
	virtual void computeForceBonds(const uint size, double* e, double* fx, double* fy, double* fz,
			const double* rx, const double* ry, const double* rz, double** params, FF::FunctForm form) const;
	virtual void computeForceAngles(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, double** params, FF::FunctForm form) const;
	virtual void computeForceDihedrals(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, const double* rx34, const double* ry34, const double* rz34,
			double** params, FF::FunctForm form) const;
	virtual void computeForceImpropers(const uint size, double* e, double* fx2, double* fy2, double* fz2, double* fx3, double* fy3, double* fz3, double* fx4, double* fy4, double* fz4,
			const double* rx21, const double* ry21, const double* rz21, const double* rx31, const double* ry31, const double* rz31, const double* rx41, const double* ry41, const double* rz41,
			double** params, FF::FunctForm form) const;
	virtual void setParam(const uint8_t kindDof, const std::string typeDof);
	virtual void getParam(const std::string typeDof, std::vector<double>& params) const;
	virtual void beginInit();
	virtual void endInit();
	/// @brief Get the type of an atom
	/// @param [in] subtype Subtype of said atom
	inline virtual std::string getAtomType(std::string subtype) const {return (*m_initOfTypes)[subtype];}
	virtual IPotential* makePotential(std::string typei, std::string typej) const;

private :

	std::unordered_map<std::string,ForceFieldParam> m_parameters; ///< Parameters for this forcefield
	std::unordered_map<std::string,std::vector<double>>* m_initOfBondedParams; ///< Pointer to a temporary array to store the data used to build the parameters for bonded interactions
	std::unordered_map<std::string,std::vector<double>>* m_initOfNonBondedParams; ///< Pointer to a temporary array to store the data used to build the parameters for non bonded interactions
	std::unordered_map<std::string,std::string>* m_initOfTypes; ///< Pointer to a temporary array to store the type/subtype correspondence
	std::ifstream m_initFile; ///< Parameters data file
	simd::kernels::bharm<double> m_computeBond; ///< Vectorized bond forces calculation
	simd::kernels::acosto3<double> m_computeAngle; ///< Vectorized angle forces calculation
	simd::kernels::dcosn<double> m_computeDihedral; ///< Vectorized dihedral forces calculation
	simd::kernels::icosto2<double> m_computeImproper; ///< Vectorized improper forces calculation

};

#endif /* FORCEFIELDUFF_HPP_ */
