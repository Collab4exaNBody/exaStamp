/*
 * forceFieldUFF.cpp
 *
 *  Created on: Mar 21, 2016
 *      Author: giarda
 */


#include <math.h>

#include "forceField/forceFieldUFF.hpp"

#include "potential/lennardJonesPotential.hpp"

#include "utils/stampUnits.hpp"
#include "utils/stringUtils.hpp"


#ifndef __xsp_compute_rep
/// @brief Defines default slice value
#define __xsp_compute_rep ""
#endif


/// @brief Default constructor
ForceFieldUFF::ForceFieldUFF() :
	ForceField(FF::UFF,3)
	{
	m_costs.push_back(1); // Cost for bonds
	m_costs.push_back(1); // Cost for angles
	m_costs.push_back(1); // Cost for dihedral angles
	}


/// @brief Destructor (nothing to do)
ForceFieldUFF::~ForceFieldUFF() {}


/// @brief Begin the initialization of the force field
void ForceFieldUFF::beginInit() {
	// Allocate the initialization maps
	m_initOfBondedParams=new std::unordered_map<std::string,std::vector<double>>();
	m_initOfNonBondedParams=new std::unordered_map<std::string,std::vector<double>>();
	m_initOfTypes=new std::unordered_map<std::string,std::string>();
	// Open and check the parameters file
	std::string repname=__xsp_compute_rep;
	m_initFile.open(repname+"/data/UFF.dat");
	if(!m_initFile.is_open()) {
		std::cerr << "Error in ForceFieldUFF::beginInit : file UFF.dat not found !" << std::endl;
		exit(-1);
	}
	std::string line;
	// Pass legend line
	getline(m_initFile,line);
	// For each non-empty line of bonded parameters part
	while(m_initFile.good()) {
		// Get the line
		getTrimedLine(m_initFile,line);
		// If the end of the part is found, get out of this loop
		if(is_equal(line,"endbonded")) break;
		// Else do the work (if the line is not empty)
		else if(line.size()!=0) {
			std::string subtype, type;
			std::vector<double> values;
			// Get the string until the next space or tabulation as subtype
			subtype=line.substr(0,line.find_first_of(" \t"));
			line=line.substr(line.find_first_of(" \t"));
			// Get 6 times the string until next space or tabulation as value
			for (int iValue(0);iValue<6;++iValue) {
				std::string strValue;
				line=line.substr(line.find_first_not_of(" \t")); // Delete extra spaces and tabulations before the next value
				strValue=line.substr(0,line.find_first_of(" \t"));
				line=line.substr(line.find_first_of(" \t"));
				values.push_back(stod(strValue));
			} // End of values
			// Set these values as parameters for this subtype
			m_initOfBondedParams->insert(make_pair(subtype,values));
			// Get the last string of the line as type for this subtype
			trim_spaces(line, "\t");
			type=line;
			m_initOfTypes->insert(make_pair(subtype,type));
		} // End of the work on the line
	} // End of the bonded parameters part
	// Pass legend line
	getline(m_initFile,line);
	// For each non-empty line of non bonded parameters part
	while(m_initFile.good()) {
		// Get the line
		getTrimedLine(m_initFile,line);
		// If the end of the part is found, get out of this loop
		if(is_equal(line,"endnonbonded")) break;
		// Else do the work (if the line is not empty)
		else if(line.size()!=0) {
			std::string type;
			std::vector<double> values(3,0.);
			// Get the string until the next space or tabulation as type
			type=line.substr(0,line.find_first_of(" \t"));
			line=line.substr(line.find_first_of(" \t"));
			// Get 3 times the string until next space or tabulation as value
			for (int iValue(0);iValue<3;++iValue) {
				std::string strValue;
				strValue=line.substr(line.find_last_of(" \t")+1);	// Get last value
				line=line.substr(0,line.find_last_of(" \t"));
				line=line.substr(0,line.find_last_not_of(" \t")+1); // Delete extra spaces at the end
				values[2-iValue]=stod(strValue);
			} // End of values
			// Set these values as parameters for this type
			m_initOfNonBondedParams->insert(make_pair(type,values));
		} // End of the work on the line
	} // End of the non bonded parameters part
	// Close the file
	m_initFile.close();
}


/// @brief End the initialization of the force field
void ForceFieldUFF::endInit() {
	// Delete the initialization maps
	delete m_initOfBondedParams;
	delete m_initOfNonBondedParams;
	delete m_initOfTypes;
}


/// @brief Create the parameters for a type of dof
/// @param [in] kindDof Kind of dof
/// @param [in] typeDof Type of the dof
void ForceFieldUFF::setParam(const uint8_t kindDof, const std::string typeDof) {
	// Create parameters only if there is no parameters for this dof
	if(m_parameters.find(typeDof)==m_parameters.end()) {
		// Initialize the parameters
		m_parameters[typeDof]=ForceFieldParam();
		// Lambda function that calculate the equilibrium bond length for two atoms
		// Parameters : r1 -> radius for atom 1
		//							r2 -> radius for atom 2
		//							X1 -> electronegativity for atom 1
		//							X2 -> electronegativity for atom 2
		// 							bOrder -> bond order
		auto dr = [&] (double r1, double r2, double X1, double X2, double bOrder) -> double {
			double rBO=-0.1332*(r1+r2)*auxLog(bOrder); // Bond order component
			double rEN=r1*r2*auxSq((auxSqrt(X1)-auxSqrt(X2)))/(X1*r1+X2*r2); // Electronegativity component
			return r1+r2+rBO+rEN;
		};
		// Lambda function to determine if an atom is from group 6
		// Parameter : typeAtom -> type of the atom
		auto group6 = [] (std::string typeAtom) -> bool {
			std::string element=typeAtom.substr(0,2);
			return (is_equal(element,"O_")||is_equal(element,"S_")||is_equal(element,"Se")||is_equal(element,"Te")||is_equal(element,"Po"));
		};
		// Parameters construction depends on the kind of dof
		switch(kindDof) {
		// Case bond
		case 0 :
		{
			// Data defining the type of dof
			std::string typeAtomi, typeAtomj, tmp(typeDof);
			double bondOrder;
			// Get them
			typeAtomi=tmp.substr(0,tmp.find('*'));
			tmp=tmp.substr(tmp.find('*')+1);
			bondOrder=stod(tmp.substr(0,tmp.find('*')));
			typeAtomj=tmp.substr(tmp.find('*')+1);
			// Initialization parameters from the file
			double ri((*m_initOfBondedParams)[typeAtomi][0]), rj((*m_initOfBondedParams)[typeAtomj][0]);
			double Zi((*m_initOfBondedParams)[typeAtomi][2]), Zj((*m_initOfBondedParams)[typeAtomj][2]);
			double Xi((*m_initOfBondedParams)[typeAtomi][5]), Xj((*m_initOfBondedParams)[typeAtomj][5]);
			// Parameters of the force field
			double rij, Kij;
			// Calculate rij
			rij=dr(ri,rj,Xi,Xj,bondOrder);
			// Calculate Kij
			Kij=664.12*Zi*Zj/(rij*rij*rij);
			// Add to the stored parameters
			m_parameters[typeDof].addParam(convert(Kij/2,SI_Units_base::kcalPerMol/(SI_Units_base::angstrom*SI_Units_base::angstrom),Stamp_Units::energy/(Stamp_Units::length*Stamp_Units::length)));
			m_parameters[typeDof].addParam(convert(rij,SI_Units_base::angstrom,Stamp_Units::length));
			break;
		} // End of bond case
		// Case angle
		case 1 :
		{
			// Data defining the type of dof
			std::string typeAtomi, typeAtomj, typeAtomk, tmp(typeDof);
			double bondOrderij, bondOrderjk;
			// Get them
			typeAtomi=tmp.substr(0,tmp.find('*'));
			tmp=tmp.substr(tmp.find('*')+1);
			bondOrderij=stod(tmp.substr(0,tmp.find('*')));
			tmp=tmp.substr(tmp.find('*')+1);
			typeAtomj=tmp.substr(0,tmp.find('*'));
			tmp=tmp.substr(tmp.find('*')+1);
			bondOrderjk=stod(tmp.substr(0,tmp.find('*')));
			typeAtomk=tmp.substr(tmp.find('*')+1);
			// Initialization parameters from the file
			double ri((*m_initOfBondedParams)[typeAtomi][0]), rj((*m_initOfBondedParams)[typeAtomj][0]), rk((*m_initOfBondedParams)[typeAtomk][0]);
			double Zi((*m_initOfBondedParams)[typeAtomi][2]), Zk((*m_initOfBondedParams)[typeAtomk][2]);
			double Xi((*m_initOfBondedParams)[typeAtomi][5]), Xj((*m_initOfBondedParams)[typeAtomj][5]), Xk((*m_initOfBondedParams)[typeAtomk][5]);
			double theta0((*m_initOfBondedParams)[typeAtomj][1]);
			theta0=Constant::pi*theta0/180;
			// Parameters of the force field
			double Kijk, C0, C1, C2, C3;
			// Calculate Kijk
			double rij=dr(ri,rj,Xi,Xj,bondOrderij);
			double rjk=dr(rj,rk,Xj,Xk,bondOrderjk);
			double rik=auxSqrt(auxSq(rij)+auxSq(rjk)-2*rij*rjk*auxCos(theta0));
			Kijk=664.12*Zi*Zk*(3*rij*rjk*(1-auxSq(auxCos(theta0)))-rik*auxCos(theta0))/(auxPow(rik,5));
			// Calculate de C0, C1, C2 and n
			// Depend on the hybridization or geometry of the central atom
			char type=typeAtomj[2];
			switch(type) {
			// Case sp1
			case '1' :
				C0=1;
				C1=1; // Fixed with open babel form
				C2=0;
				C3=0;
				break;
			// Case sp2 or aromatic
			case '2' :
			case 'R' :
				C0=2./3;
				C1=8./9;
				C2=4./9;
				C3=0;
				break;
			// Case square planar or octahedral
			case '4' :
			case '6' :
				C0=1./2;
				C1=3./4;
				C2=1./2;
				C3=1./4;
				break;
			// Other cases (including sp3 and trigonal bipyramidal)
			default :
				C2=1./(4*auxSq(auxSin(theta0)));
				C1=-4*C2*auxCos(theta0);
				C0=C2*(2*auxSq(auxCos(theta0))+1);
				C3=0;
				break;
			}
			// Add to the stored parameters
			m_parameters[typeDof].addParam(convert(Kijk*C0,SI_Units_base::kcalPerMol,Stamp_Units::energy));
			m_parameters[typeDof].addParam(convert(Kijk*C1,SI_Units_base::kcalPerMol,Stamp_Units::energy));
			m_parameters[typeDof].addParam(convert(Kijk*C2,SI_Units_base::kcalPerMol,Stamp_Units::energy));
			m_parameters[typeDof].addParam(convert(Kijk*C3,SI_Units_base::kcalPerMol,Stamp_Units::energy));
			break;
		} // End of angle case
		// Case torsion
		case 2 :
		{
			// Data defining the type of dof
			std::string typeAtomj, typeAtomk, tmp(typeDof);
			double bondOrder;
			bool pc;
			// Get them
			tmp=tmp.substr(1,tmp.size());
			typeAtomj=tmp.substr(0,tmp.find('*'));
			tmp=tmp.substr(tmp.find('*')+1);
			bondOrder=stod(tmp.substr(0,tmp.find('*')));
			typeAtomk=tmp.substr(tmp.find('*')+1);
			tmp=tmp.substr(tmp.find('*')+1);
			pc=(tmp.size()!=0);
			// Parameters of the force field
			double Vphi, cosnphi0, n;
			// Get the coordination for the two atoms involved
			char typej(typeAtomj[2]), typek(typeAtomk[2]);
			// Case sp2-sp2
			if((typej=='2'||typej=='R')&&(typek=='2'||typek=='R')) {
				// Initialization parameters from the file
				double Uj((*m_initOfBondedParams)[typeAtomj][4]), Uk((*m_initOfBondedParams)[typeAtomk][4]);
				// Calculate Vphi, cosnphi0 and n
				Vphi=5*auxSqrt(Uj*Uk)*(1+4.18*auxLog(bondOrder));
				n=2;
				cosnphi0=1; // phi0=180°
			}
			// Case sp3-sp3
			else if((typej=='3')&&(typek=='3')) {
				// Initialization parameters from the file
				double Vj((*m_initOfBondedParams)[typeAtomj][3]), Vk((*m_initOfBondedParams)[typeAtomk][3]);
				// Calculate cosnphi0 and n
				n=3;
				cosnphi0=-1; // phi0=180°
				if(group6(typeAtomj)&&group6(typeAtomk)) {
					Vj=6.8;
					Vk=6.8;
					if(is_equal(typeAtomj.substr(0,2),"O_")) Vj=2;
					if(is_equal(typeAtomk.substr(0,2),"O_")) Vk=2;
					// Reset n (phi0=90° so cosnphi0=-1)
					n=2;
				}
				// Calculate Vphi
				Vphi=auxSqrt(Vj*Vk);
			}
			// Case sp2-sp3
			else {
				if(pc) {
					// Calculate Vphi, cosnphi0 and n
					Vphi=2;
					n=3;
					cosnphi0=-1; // phi0=180°
				}
				else {
					// Calculate Vphi
					Vphi=1;
					if(group6(typeAtomj)&&group6(typeAtomk)) {
						// Calculate n and cosnphi0
						n=2;
						cosnphi0=-1; // phi0=90°
					}
					else {
						n=6;
						cosnphi0=1; // phi0=0°
					}
				}
			}
			// Add to the stored parameters
			m_parameters[typeDof].addParam(convert(Vphi/2,SI_Units_base::kcalPerMol,Stamp_Units::energy));
			m_parameters[typeDof].addParam(convert(-Vphi*cosnphi0/2,SI_Units_base::kcalPerMol,Stamp_Units::energy));
			m_parameters[typeDof].addParam(n);
			break;
		} // End of torsion case
		// Case improper torsion
		case 3 :
		{
			// Parameters of the force field
			double K, C0, C1, C2;
			// Case carbon
			if(is_equal(typeDof.substr(0,2),"C_")) {
				C0=1;
				C1=-1;
				C2=0;
				K=6;
				if(is_equal(typeDof,"C_*O")) K=50;
			}
			// Case group 5 and 6
			else {
				// Natural inversion angle and its cosine
				double w0, cosw0;
				// Oxygen and nitrogen
				if(is_equal(typeDof,"N_")||is_equal(typeDof,"O_")) {
					K=6;
					w0=59.9;
				}
				// Phosphorus and sulfur
				else if(is_equal(typeDof,"P_")||is_equal(typeDof,"S_")) {
					K=22;
					w0=84.4339;
				}
				// Arsenic and selenium
				else if(is_equal(typeDof,"As")|is_equal(typeDof,"Se")) {
					K=22;
					w0=86.9735;
				}
				// Antimony and tellurium
				else if(is_equal(typeDof,"Sb")||is_equal(typeDof,"Te")) {
					K=22;
					w0=87.7047;
				}
				// Bismuth and polonium
				else if(is_equal(typeDof,"Bi")||is_equal(typeDof,"Po")) {
					K=22;
					w0=90;
				}
				else {
					K=0;
					w0=0;
				}
				cosw0=auxCos(w0);
				C2=1./(2*cosw0*cosw0-4*cosw0+2);
				C1=-4*cosw0*C2;
				C0=1-C1-C2;
			}
			// Add to the stored parameters
			m_parameters[typeDof].addParam(convert(K*C0,SI_Units_base::kcalPerMol,Stamp_Units::energy));
			m_parameters[typeDof].addParam(convert(K*C1,SI_Units_base::kcalPerMol,Stamp_Units::energy));
			m_parameters[typeDof].addParam(convert(K*C2,SI_Units_base::kcalPerMol,Stamp_Units::energy));

			break;
		} // End of improper torsion case
		} // All kinds of dof treated
	} // End of construction
}


/// @brief Create the potential between to given atom types
/// @param [in] typei First atom type
/// @param [in] typej Second atom type
/// @return Created potential
IPotential* ForceFieldUFF::makePotential(std::string typei, std::string typej) const {
	// Parameters of the potential
	double epsilon, sigma, rcut;
	// Get the initialization parameters
	double xi((*m_initOfNonBondedParams)[typei][0]), xj((*m_initOfNonBondedParams)[typej][0]);
	double Di((*m_initOfNonBondedParams)[typei][1]), Dj((*m_initOfNonBondedParams)[typej][1]);
	// Calculate the parameters and convert to stamp units
	epsilon=auxSqrt(Di*Dj);
	epsilon=convert(epsilon,SI_Units_base::kcalPerMol,Stamp_Units::energy);
	sigma=auxSqrt(xi*xj)/auxPow(2,1./6);
	sigma=convert(sigma,SI_Units_base::angstrom,Stamp_Units::length);
	rcut=2.5*sigma;
	// Build and return the potential
	return new LennardJonesPotential(rcut,LennardJonesParameters(epsilon,sigma));
}


/// @brief Get the parameters for a dof in a vectorization buffer
/// @param [in] typeDof Type of the dof
/// @param [out] params Parameters
void ForceFieldUFF::getParam(const std::string typeDof, std::vector<double>& params) const {
	m_parameters.at(typeDof).getAll(params);
}


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
void ForceFieldUFF::computeForceBonds(const uint size, double* e, double* fx, double* fy, double* fz,
		const double* rx, const double* ry, const double* rz, double** params, FF::FunctForm form) const {
	// Sequential version
//	double r, drx, dry, drz;
//	for(uint i(0); i<size; ++i) {
//		r=sqrt(rx*rx+ry*ry+rz*rz);
//		drx=rx/r;
//		dry=ry/r;
//		drz=rz/r;
//		e[i]=params[0][i]*(r-params[1][i])*(r-params[1][i]);
//		fx[i]=-2*params[0][i]*(r-params[1][i])*drx;
//		fy[i]=-2*params[0][i]*(r-params[1][i])*dry;
//		fz[i]=-2*params[0][i]*(r-params[1][i])*drz;
//	}
	// SIMD version
	m_computeBond(fx, fy, fz, e, rx, ry, rz, params[0], params[1], size);
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
/// @param [in] form Functional form for the force calculation (unused in non expert force fields)
void ForceFieldUFF::computeForceAngles(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2,
		const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, double** params, FF::FunctForm form) const {
	// Sequential version
		// Energy
//	auto Eseq = [&] (double rx1, double ry1, double rz1, double rx3, double ry3, double rz3, int i) -> double {
//		double scalar=rx1*rx3+ry1*ry3+rz1*rz3;
//		double r21=sqrt(rx1*rx1+ry1*ry1+rz1*rz1);
//		double r23=sqrt(rx3*rx3+ry3*ry3+rz3*rz3);
//		double theta=acos(scalar/(r21*r23));
//		return params[0][i]+params[1][i]*cos(theta)+params[2][i]*cos(params[3][i]*theta);
//	};
//	// Theta
//	auto theta = [&] (double rx1, double ry1, double rz1, double rx3, double ry3, double rz3) -> double {
//		double scalar=rx1*rx3+ry1*ry3+rz1*rz3;
//		double r21=sqrt(rx1*rx1+ry1*ry1+rz1*rz1);
//		double r23=sqrt(rx3*rx3+ry3*ry3+rz3*rz3);
//		return acos(scalar/(r21*r23));
//	};
	// SIMD version
	m_computeAngle(fx1, fy1, fz1, fx2, fy2, fz2, e, rx21, ry21, rz21, rx23, ry23, rz23, params[0], params[1], params[2], params[3], size);
	// Numerical forces
//	double epsilon(0.0001), r[size], ep[size], em[size], ftrash1[size], ftrash2[size], ftrash3[size], ftrash4[size], ftrash5[size], ftrash6[size], fx1n[size], fy1n[size], fz1n[size], fx2n[size], fy2n[size], fz2n[size];
//	for(uint i(0); i<size; ++i){
//		r[i]=rx21[i]+epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ep, r, ry21, rz21, rx23, ry23, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rx21[i]-epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, em, r, ry21, rz21, rx23, ry23, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		fx1n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=ry21[i]+epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ep, rx21, r, rz21, rx23, ry23, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=ry21[i]-epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, em, rx21, r, rz21, rx23, ry23, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		fy1n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rz21[i]+epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ep, rx21, ry21, r, rx23, ry23, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rz21[i]-epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, em, rx21, ry21, r, rx23, ry23, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		fz1n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rx23[i]+epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ep, rx21, ry21, rz21, r, ry23, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rx23[i]-epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, em, rx21, ry21, rz21, r, ry23, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		fx2n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=ry23[i]+epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ep, rx21, ry21, rz21, rx23, r, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=ry23[i]-epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, em, rx21, ry21, rz21, rx23, r, rz23, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		fy2n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rz23[i]+epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ep, rx21, ry21, rz21, rx23, ry23, r, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rz23[i]-epsilon;
//	}
//	m_computeAngle(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, em, rx21, ry21, rz21, rx23, ry23, r, params[0], params[1], params[2], params[3], size);
//	for(uint i(0); i<size; ++i){
//		fz2n[i]=(em[i]-ep[i])/(2*epsilon);
//		/* Debug print */ std::cerr << "Analytical forces : f1 = " << fx1[i] << ", " << fy1[i] << ", " << fz1[i] << ", f2 = " << fx2[i] << ", " << fy2[i] << ", " << fz2[i] << std::endl;
//		/* Debug print */ std::cerr << "Numerical forces : f1 = " << fx1n[i] << ", " << fy1n[i] << ", " << fz1n[i] << ", f2 = " << fx2n[i] << ", " << fy2n[i] << ", " << fz2n[i] << std::endl;
//	}
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
/// @param [in] form Functional form for the force calculation (unused in non expert force fields)
void ForceFieldUFF::computeForceDihedrals(const uint size, double* e, double* fx1, double* fy1, double* fz1, double* fx2, double* fy2, double* fz2, double* fx4, double* fy4, double* fz4,
		const double* rx21, const double* ry21, const double* rz21, const double* rx23, const double* ry23, const double* rz23, const double* rx34, const double* ry34, const double* rz34,
		double** params, FF::FunctForm form) const {
	// Sequential version
	// Phi
//	auto phi = [&] (double rx1, double ry1, double rz1, double rx2, double ry2, double rz2, double rx3, double ry3, double rz3) -> double {
//		double xA, yA, zA, xB, yB, zB;
//		xA=ry1*rz2-rz1*ry2;
//		yA=rz1*rx2-rx1*rz2;
//		zA=rx1*ry2-ry1*rx2;
//		xB=ry3*rz2-rz3*ry2;
//		yB=rz3*rx2-rx3*rz2;
//		zB=rx3*ry2-ry3*rx2;
//		double scalar=xA*xB+yA*yB+zA*zB;
//		double rA=sqrt(xA*xA+yA*yA+zA*zA);
//		double rB=sqrt(xB*xB+yB*yB+zB*zB);
//		return acos(scalar/(rA*rB));
//	};
	// SIMD version
	m_computeDihedral(fx1, fy1, fz1, fx2, fy2, fz2, fx4, fy4, fz4, e, rx21, ry21, rz21, rx23, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);

	// Numerical forces
//	double epsilon(0.00001), r[size], rr[size];
//	double ep[size], em[size], ftrash1[size], ftrash2[size], ftrash3[size], ftrash4[size], ftrash5[size], ftrash6[size], ftrash7[size], ftrash8[size], ftrash9[size];
//	double fx1n[size], fy1n[size], fz1n[size], fx2n[size], fy2n[size], fz2n[size], fx4n[size], fy4n[size], fz4n[size];
//	for(uint i(0); i<size; ++i){
//		r[i]=rx21[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, r, ry21, rz21, rx23, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rx21[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, r, ry21, rz21, rx23, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fx1n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=ry21[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, r, rz21, rx23, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=ry21[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, r, rz21, rx23, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fy1n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rz21[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, r, rx23, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rz21[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, r, rx23, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fz1n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rx23[i]-epsilon;
//		rr[i]=rx21[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rr, ry21, rz21, r, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rx23[i]+epsilon;
//		rr[i]=rx21[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rr, ry21, rz21, r, ry23, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fx2n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=ry23[i]-epsilon;
//		rr[i]=ry21[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, rr, rz21, rx23, r, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=ry23[i]+epsilon;
//		rr[i]=ry21[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, rr, rz21, rx23, r, rz23, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fy2n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rz23[i]-epsilon;
//		rr[i]=rz21[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rr, rx23, ry23, r, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rz23[i]+epsilon;
//		rr[i]=rz21[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rr, rx23, ry23, r, rx34, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fz2n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rx34[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, rx23, ry23, rz23, r, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rx34[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, rx23, ry23, rz23, r, ry34, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fx4n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=ry34[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, rx23, ry23, rz23, rx34, r, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=ry34[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, rx23, ry23, rz23, rx34, r, rz34, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fy4n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rz34[i]+epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, rx23, ry23, rz23, rx34, ry34, r, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rz34[i]-epsilon;
//	}
//	m_computeDihedral(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, rx23, ry23, rz23, rx34, ry34, r, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fz4n[i]=(em[i]-ep[i])/(2*epsilon);
//		/* Debug print */ std::cerr << i << " Analytical forces : f1 = " << fx1[i] << ", " << fy1[i] << ", " << fz1[i] << ", f2 = " << fx2[i] << ", " << fy2[i] << ", " << fz2[i] << ", f4 = " << fx4[i] << ", " << fy4[i] << ", " << fz4[i] << std::endl;
//		/* Debug print */ std::cerr << i << " Numerical forces : f1 = " << fx1n[i] << ", " << fy1n[i] << ", " << fz1n[i] << ", f2 = " << fx2n[i] << ", " << fy2n[i] << ", " << fz2n[i] << ", f4 = " << fx4n[i] << ", " << fy4n[i] << ", " << fz4n[i] << std::endl;
//	}
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
/// @param [in] form Functional form for the force calculation (unused in non expert force fields)
void ForceFieldUFF::computeForceImpropers(const uint size, double* e, double* fx2, double* fy2, double* fz2, double* fx3, double* fy3, double* fz3, double* fx4, double* fy4, double* fz4,
		const double* rx21, const double* ry21, const double* rz21, const double* rx31, const double* ry31, const double* rz31, const double* rx41, const double* ry41, const double* rz41,
		double** params, FF::FunctForm form) const {
	// Sequential version
	// SIMD version
	m_computeImproper(fx2, fy2, fz2, fx3, fy3, fz3, fx4, fy4, fz4, e, rx21, ry21, rz21, rx31, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);

//	// Numerical forces
//	double epsilon(0.00001), r[size];
//	double ep[size], em[size], ftrash1[size], ftrash2[size], ftrash3[size], ftrash4[size], ftrash5[size], ftrash6[size], ftrash7[size], ftrash8[size], ftrash9[size];
//	double fx2n[size], fy2n[size], fz2n[size], fx3n[size], fy3n[size], fz3n[size], fx4n[size], fy4n[size], fz4n[size];
//	for(uint i(0); i<size; ++i){
//		r[i]=rx21[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, r, ry21, rz21, rx31, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rx21[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, r, ry21, rz21, rx31, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fx2n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=ry21[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, r, rz21, rx31, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=ry21[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, r, rz21, rx31, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fy2n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rz21[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, r, rx31, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rz21[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, r, rx31, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fz2n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rx31[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, r, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rx31[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, r, ry31, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fx3n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=ry31[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, rx31, r, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=ry31[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, rx31, r, rz31, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fy3n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rz31[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, rx31, ry31, r, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rz31[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, rx31, ry31, r, rx41, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fz3n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rx41[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, rx31, ry31, rz31, r, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rx41[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, rx31, ry31, rz31, r, ry41, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fx4n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=ry41[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, rx31, ry31, rz31, rx41, r, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=ry41[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, rx31, ry31, rz31, rx41, r, rz41, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fy4n[i]=(em[i]-ep[i])/(2*epsilon);
//		r[i]=rz41[i]+epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, ep, rx21, ry21, rz21, rx31, ry31, rz31, rx41, ry41, r, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		r[i]=rz41[i]-epsilon;
//	}
//	m_computeImproper(ftrash1, ftrash2, ftrash3, ftrash4, ftrash5, ftrash6, ftrash7, ftrash8, ftrash9, em, rx21, ry21, rz21, rx31, ry31, rz31, rx41, ry41, r, params[0], params[1], params[2], size);
//	for(uint i(0); i<size; ++i){
//		fz4n[i]=(em[i]-ep[i])/(2*epsilon);
//		/* Debug print */ std::cerr << i << " Analytical forces : f2 = " << fx2[i] << ", " << fy2[i] << ", " << fz2[i] << ", f3 = " << fx3[i] << ", " << fy3[i] << ", " << fz3[i] << ", f4 = " << fx4[i] << ", " << fy4[i] << ", " << fz4[i] << std::endl;
//		/* Debug print */ std::cerr << i << " Numerical forces : f2 = " << fx2n[i] << ", " << fy2n[i] << ", " << fz2n[i] << ", f3 = " << fx3n[i] << ", " << fy3n[i] << ", " << fz3n[i] << ", f4 = " << fx4n[i] << ", " << fy4n[i] << ", " << fz4n[i] << std::endl;
//	}
}
