/*
 * vectorizationMBuffer.hpp
 *
 *  Created on: Jan 31, 2017
 *      Author: giarda
 */

/// @file
/// @brief Definition of the class VectMBuffer

#ifndef VECTORIZATIONMBUFFER_HPP_
#define VECTORIZATIONMBUFFER_HPP_


#include "globals.hpp"

#include "forceField/forceField.hpp"

#include "utils/array/extArray.hpp"


/// @brief Provide memory-aligned arrays to perform vectorized operations in
/// intramolecular force computation
///
/// Number of parameters for the intramolecular force calculation will change from one force
/// field to the other, therefore, feel free to add enough storage for your most expensive
/// force and to add cases in the resize functions accordingly
/// @tparam align Alignment
template <size_t align> class VectMBuffer {

public :

  /// @brief Default constructor
	VectMBuffer() {}
	/// @brief Destructor (nothing to do)
	~VectMBuffer() {}

  /// @brief Get a pointer to the x coordinate of the first atoms
  inline double* rx1() { return m_rx1.data();}
  /// @brief Get a constant pointer to the x coordinate of the first atoms
  inline const double* rx1() const { return m_rx1.data();}
  /// @brief Get a reference to the x coordinate of a specified first atom
  /// @param [in] i Particle index
  inline double& rx1(const uint i) { return m_rx1[i];}
  /// @brief Get a constant reference to the x coordinate of a specified first atom
  /// @param [in] i Particle index
  inline const double& rx1(const uint i) const { return m_rx1[i];}

  /// @brief Get a pointer to the y coordinate of the first atoms
  inline double* ry1() { return m_ry1.data();}
  /// @brief Get a constant pointer to the y coordinate of the first atoms
  inline const double* ry1() const { return m_ry1.data();}
  /// @brief Get a reference to the y coordinate of a specified first atom
  /// @param [in] i Particle index
  inline double& ry1(const uint i) { return m_ry1[i];}
  /// @brief Get a constant reference to the y coordinate of a specified first atom
  /// @param [in] i Particle index
  inline const double& ry1(const uint i) const { return m_ry1[i];}

  /// @brief Get a pointer to the z coordinate of the first atoms
  inline double* rz1() { return m_rz1.data();}
  /// @brief Get a constant pointer to the z coordinate of the first atoms
  inline const double* rz1() const { return m_rz1.data();}
  /// @brief Get a reference to the z coordinate of a specified first atom
  /// @param [in] i Particle index
  inline double& rz1(const uint i) { return m_rz1[i];}
  /// @brief Get a constant reference to the z coordinate of a specified first atom
  /// @param [in] i Particle index
  inline const double& rz1(const uint i) const { return m_rz1[i];}

  /// @brief Get a pointer to the x coordinate of the second atoms
  inline double* rx2() { return m_rx2.data();}
  /// @brief Get a constant pointer to the x coordinate of the second atoms
  inline const double* rx2() const { return m_rx2.data();}
  /// @brief Get a reference to the x coordinate of a specified second atom
  /// @param [in] i Particle index
  inline double& rx2(const uint i) { return m_rx2[i];}
  /// @brief Get a constant reference to the x coordinate of a specified second atom
  /// @param [in] i Particle index
  inline const double& rx2(const uint i) const { return m_rx2[i];}

  /// @brief Get a pointer to the y coordinate of the second atoms
  inline double* ry2() { return m_ry2.data();}
  /// @brief Get a constant pointer to the y coordinate of the second atoms
  inline const double* ry2() const { return m_ry2.data();}
  /// @brief Get a reference to the y coordinate of a specified second atom
  /// @param [in] i Particle index
  inline double& ry2(const uint i) { return m_ry2[i];}
  /// @brief Get a constant reference to the y coordinate of a specified second atom
  /// @param [in] i Particle index
  inline const double& ry2(const uint i) const { return m_ry2[i];}

  /// @brief Get a pointer to the z coordinate of the second atoms
  inline double* rz2() { return m_rz2.data();}
  /// @brief Get a constant pointer to the z coordinate of the second atoms
  inline const double* rz2() const { return m_rz2.data();}
  /// @brief Get a reference to the z coordinate of a specified second atom
  /// @param [in] i Particle index
  inline double& rz2(const uint i) { return m_rz2[i];}
  /// @brief Get a constant reference to the z coordinate of a specified second atom
  /// @param [in] i Particle index
  inline const double& rz2(const uint i) const { return m_rz2[i];}

  /// @brief Get a pointer to the x coordinate of the third atoms
  inline double* rx3() { return m_rx3.data();}
  /// @brief Get a constant pointer to the x coordinate of the third atoms
  inline const double* rx3() const { return m_rx3.data();}
  /// @brief Get a reference to the x coordinate of a specified third atom
  /// @param [in] i Particle index
  inline double& rx3(const uint i) { return m_rx3[i];}
  /// @brief Get a constant reference to the x coordinate of a specified third atom
  /// @param [in] i Particle index
  inline const double& rx3(const uint i) const { return m_rx3[i];}

  /// @brief Get a pointer to the y coordinate of the third atoms
  inline double* ry3() { return m_ry3.data();}
  /// @brief Get a constant pointer to the y coordinate of the third atoms
  inline const double* ry3() const { return m_ry3.data();}
  /// @brief Get a reference to the y coordinate of a specified third atom
  /// @param [in] i Particle index
  inline double& ry3(const uint i) { return m_ry3[i];}
  /// @brief Get a constant reference to the y coordinate of a specified third atom
  /// @param [in] i Particle index
  inline const double& ry3(const uint i) const { return m_ry3[i];}

  /// @brief Get a pointer to the z coordinate of the third atoms
  inline double* rz3() { return m_rz3.data();}
  /// @brief Get a constant pointer to the z coordinate of the third atoms
  inline const double* rz3() const { return m_rz3.data();}
  /// @brief Get a reference to the z coordinate of a specified third atom
  /// @param [in] i Particle index
  inline double& rz3(const uint i) { return m_rz3[i];}
  /// @brief Get a constant reference to the z coordinate of a specified third atom
  /// @param [in] i Particle index
  inline const double& rz3(const uint i) const { return m_rz3[i];}

  /// @brief Get a pointer to the x components of the first forces
  inline double* fx1() { return m_fx1.data();}
  /// @brief Get a constant pointer to the x components of the first forces
  inline const double* fx1() const { return m_fx1.data();}
  /// @brief Get a reference to the x components of the first force for a specified dof
  /// @param [in] i Particle index
  inline double& fx1(const uint i) { return m_fx1[i];}
  /// @brief Get a constant reference to the x components of the first force for a specified dof
  /// @param [in] i Particle index
  inline const double& fx1(const uint i) const { return m_fx1[i];}

  /// @brief Get a pointer to the y components of the first forces
  inline double* fy1() { return m_fy1.data();}
  /// @brief Get a constant pointer to the y components of the first forces
  inline const double* fy1() const { return m_fy1.data();}
  /// @brief Get a reference to the y components of the first force for a specified dof
  /// @param [in] i Particle index
  inline double& fy1(const uint i) { return m_fy1[i];}
  /// @brief Get a constant reference to the y components of the first force for a specified dof
  /// @param [in] i Particle index
  inline const double& fy1(const uint i) const { return m_fy1[i];}

  /// @brief Get a pointer to the z components of the first forces
  inline double* fz1() { return m_fz1.data();}
  /// @brief Get a constant pointer to the z components of the first forces
  inline const double* fz1() const { return m_fz1.data();}
  /// @brief Get a reference to the z components of the first force for a specified dof
  /// @param [in] i Particle index
  inline double& fz1(const uint i) { return m_fz1[i];}
  /// @brief Get a constant reference to the z components of the first force for a specified dof
  /// @param [in] i Particle index
  inline const double& fz1(const uint i) const { return m_fz1[i];}

  /// @brief Get a pointer to the x components of the second forces
  inline double* fx2() { return m_fx2.data();}
  /// @brief Get a constant pointer to the x components of the second forces
  inline const double* fx2() const { return m_fx2.data();}
  /// @brief Get a reference to the x components of the second force for a specified dof
  /// @param [in] i Particle index
  inline double& fx2(const uint i) { return m_fx2[i];}
  /// @brief Get a constant reference to the x components of the second force for a specified dof
  /// @param [in] i Particle index
  inline const double& fx2(const uint i) const { return m_fx2[i];}

  /// @brief Get a pointer to the y components of the second forces
  inline double* fy2() { return m_fy2.data();}
  /// @brief Get a constant pointer to the y components of the second forces
  inline const double* fy2() const { return m_fy2.data();}
  /// @brief Get a reference to the y components of the second force for a specified dof
  /// @param [in] i Particle index
  inline double& fy2(const uint i) { return m_fy2[i];}
  /// @brief Get a constant reference to the y components of the second force for a specified dof
  /// @param [in] i Particle index
  inline const double& fy2(const uint i) const { return m_fy2[i];}

  /// @brief Get a pointer to the z components of the second forces
  inline double* fz2() { return m_fz2.data();}
  /// @brief Get a constant pointer to the z components of the second forces
  inline const double* fz2() const { return m_fz2.data();}
  /// @brief Get a reference to the z components of the second force for a specified dof
  /// @param [in] i Particle index
  inline double& fz2(const uint i) { return m_fz2[i];}
  /// @brief Get a constant reference to the z components of the second force for a specified dof
  /// @param [in] i Particle index
  inline const double& fz2(const uint i) const { return m_fz2[i];}

  /// @brief Get a pointer to the x components of the third forces
  inline double* fx3() { return m_fx3.data();}
  /// @brief Get a constant pointer to the x components of the third forces
  inline const double* fx3() const { return m_fx3.data();}
  /// @brief Get a reference to the x components of the third force for a specified dof
  /// @param [in] i Particle index
  inline double& fx3(const uint i) { return m_fx3[i];}
  /// @brief Get a constant reference to the x components of the third force for a specified dof
  /// @param [in] i Particle index
  inline const double& fx3(const uint i) const { return m_fx3[i];}

  /// @brief Get a pointer to the y components of the third forces
  inline double* fy3() { return m_fy3.data();}
  /// @brief Get a constant pointer to the y components of the third forces
  inline const double* fy3() const { return m_fy3.data();}
  /// @brief Get a reference to the y components of the third force for a specified dof
  /// @param [in] i Particle index
  inline double& fy3(const uint i) { return m_fy3[i];}
  /// @brief Get a constant reference to the y components of the third force for a specified dof
  /// @param [in] i Particle index
  inline const double& fy3(const uint i) const { return m_fy3[i];}

  /// @brief Get a pointer to the z components of the third forces
  inline double* fz3() { return m_fz3.data();}
  /// @brief Get a constant pointer to the z components of the third forces
  inline const double* fz3() const { return m_fz3.data();}
  /// @brief Get a reference to the z components of the third force for a specified dof
  /// @param [in] i Particle index
  inline double& fz3(const uint i) { return m_fz3[i];}
  /// @brief Get a constant reference to the z components of the third force for a specified dof
  /// @param [in] i Particle index
  inline const double& fz3(const uint i) const { return m_fz3[i];}

  /// @brief Get a pointer to the potential energies
  inline double* ep() { return m_ep.data();}
  /// @brief Get a constant pointer to the potential energies
  inline const double* ep() const { return m_ep.data();}
  /// @brief Get a reference to the potential energy for a specified dof
  /// @param [in] i Particle index
  inline double& ep(const uint i) { return m_ep[i];}
  /// @brief Get a constant reference to the potential energy for a specified dof
  /// @param [in] i Particle index
  inline const double& ep(const uint i) const { return m_ep[i];}

  /// @brief Get a pointer to the first parameters
  inline double* p1() { return m_p1.data();}
  /// @brief Get a constant pointer to the first parameters
  inline const double* p1() const { return m_p1.data();}
  /// @brief Get a reference to the first parameter for a specified dof
  /// @param [in] i Particle index
  inline double& p1(const uint i) { return m_p1[i];}
  /// @brief Get a constant reference to the first parameter for a specified dof
  /// @param [in] i Particle index
  inline const double& p1(const uint i) const { return m_p1[i];}

  /// @brief Get a pointer to the second parameters
  inline double* p2() { return m_p2.data();}
  /// @brief Get a constant pointer to the second parameters
  inline const double* p2() const { return m_p2.data();}
  /// @brief Get a reference to the second parameter for a specified dof
  /// @param [in] i Particle index
  inline double& p2(const uint i) { return m_p2[i];}
  /// @brief Get a constant reference to the second parameter for a specified dof
  /// @param [in] i Particle index
  inline const double& p2(const uint i) const { return m_p2[i];}

  /// @brief Get a pointer to the third parameters
  inline double* p3() { return m_p3.data();}
  /// @brief Get a constant pointer to the third parameters
  inline const double* p3() const { return m_p3.data();}
  /// @brief Get a reference to the third parameter for a specified dof
  /// @param [in] i Particle index
  inline double& p3(const uint i) { return m_p3[i];}
  /// @brief Get a constant reference to the third parameter for a specified dof
  /// @param [in] i Particle index
  inline const double& p3(const uint i) const { return m_p3[i];}

  /// @brief Get a pointer to the forth parameters
  inline double* p4() { return m_p4.data();}
  /// @brief Get a constant pointer to the forth parameters
  inline const double* p4() const { return m_p4.data();}
  /// @brief Get a reference to the forth parameter for a specified dof
  /// @param [in] i Particle index
  inline double& p4(const uint i) { return m_p4[i];}
  /// @brief Get a constant reference to the forth parameter for a specified dof
  /// @param [in] i Particle index
  inline const double& p4(const uint i) const { return m_p4[i];}

  /// @brief Get a array of the pointer to the parameters
  /// @param [in] num Number of parameters
  inline double** params(const uint num) {
  	m_pAccess.resize(num);
  	switch(num){
  	case 4 :
  		m_pAccess[3]=p4();
  	case 3 :
  		m_pAccess[2]=p3();
  	case 2 :
  		m_pAccess[1]=p2();
  	case 1 :
  		m_pAccess[0]=p1();
  	case 0 :
  		break;
  	}
  	return m_pAccess.data();
  }

  /// @brief Resize all arrays for a specified dofs force calculation
  /// @tparam kind Specify the dofs (0 for bonds, 1 for angles, 2 for dihedrals, 3 for impropers)
  /// @param [in] n New size
  template <uint8_t kind> void resize(const uint n) {
  	switch(kind) {
  	case 0 :
  		resizeBond(n);
  		break;
  	case 1 :
  		resizeAngle(n);
  		break;
  	case 2 :
  		resizeDihedral(n);
  		break;
  	case 3 :
  		resizeImproper(n);
  		break;
  	}
  }

  /// @brief Resize parameters arrays
  /// @param [in] n New size
  /// @param [in] num Number of parameters
  inline void resizeParam(const uint n, const uint num) {
  	switch(num) {
  	case 4 :
  		m_p4.resize(n);
  	case 3 :
  		m_p3.resize(n);
  	case 2 :
  		m_p2.resize(n);
  	case 1 :
  		m_p1.resize(n);
  	default :
  		break;
  	}
  }

  /// @brief Resize all arrays used for bond force calculations
  /// @param [in] n New size
  inline void resizeBond(const uint n) {

    // Only if size increase
    if (n > m_rx1.size()) {

      m_rx1.resize(n);
      m_ry1.resize(n);
      m_rz1.resize(n);

      m_fx1.resize(n);
      m_fy1.resize(n);
      m_fz1.resize(n);

      m_ep.resize(n);

    }

    resizeParam(n,std::get<0>(ForceField::paramBuild.at(Global::ffield->type())));

  }

  /// @brief Resize all arrays used for angle force calculations
  /// @param [in] n New size
  inline void resizeAngle(const uint n) {
    // Only if size increase
    if (n > m_rx2.size()) {

      m_rx1.resize(n);
      m_ry1.resize(n);
      m_rz1.resize(n);

      m_rx2.resize(n);
      m_ry2.resize(n);
      m_rz2.resize(n);

      m_fx1.resize(n);
      m_fy1.resize(n);
      m_fz1.resize(n);

      m_fx2.resize(n);
      m_fy2.resize(n);
      m_fz2.resize(n);

      m_ep.resize(n);

    }

    resizeParam(n,std::get<1>(ForceField::paramBuild.at(Global::ffield->type())));

  }

  /// @brief Resize all arrays used for dihedral angle force calculations
  /// @param [in] n New size
  inline void resizeDihedral(const uint n) {

    // Only if size increase
    if (n > m_rx3.size()) {

      m_rx1.resize(n);
      m_ry1.resize(n);
      m_rz1.resize(n);

      m_rx2.resize(n);
      m_ry2.resize(n);
      m_rz2.resize(n);

      m_rx3.resize(n);
      m_ry3.resize(n);
      m_rz3.resize(n);

      m_fx1.resize(n);
      m_fy1.resize(n);
      m_fz1.resize(n);

      m_fx2.resize(n);
      m_fy2.resize(n);
      m_fz2.resize(n);

      m_fx3.resize(n);
      m_fy3.resize(n);
      m_fz3.resize(n);

      m_ep.resize(n);

    }

    resizeParam(n,std::get<2>(ForceField::paramBuild.at(Global::ffield->type())));

  }

  /// @brief Resize all arrays used for improper torsion force calculations
  /// @param [in] n New size
  inline void resizeImproper(const uint n) {

    // Only if size increase
    if (n > m_rx4.size()) {

      m_rx1.resize(n);
      m_ry1.resize(n);
      m_rz1.resize(n);

      m_rx2.resize(n);
      m_ry2.resize(n);
      m_rz2.resize(n);

      m_rx3.resize(n);
      m_ry3.resize(n);
      m_rz3.resize(n);

      m_fx1.resize(n);
      m_fy1.resize(n);
      m_fz1.resize(n);

      m_fx2.resize(n);
      m_fy2.resize(n);
      m_fz2.resize(n);

      m_fx3.resize(n);
      m_fy3.resize(n);
      m_fz3.resize(n);

      m_ep.resize(n);

    }

    resizeParam(n,std::get<3>(ForceField::paramBuild.at(Global::ffield->type())));

  }

private :

  static constexpr uint8_t chunk = 16; ///< Chunk for the arrays

  ExtArray<double, chunk, align> m_rx1; ///< Storage for the x coordinate of the first atoms
  ExtArray<double, chunk, align> m_ry1; ///< Storage for the y coordinate of the first atoms
  ExtArray<double, chunk, align> m_rz1; ///< Storage for the z coordinate of the first atoms

  ExtArray<double, chunk, align> m_rx2; ///< Storage for the x coordinate of the second atoms
  ExtArray<double, chunk, align> m_ry2; ///< Storage for the y coordinate of the second atoms
  ExtArray<double, chunk, align> m_rz2; ///< Storage for the z coordinate of the second atoms

  ExtArray<double, chunk, align> m_rx3; ///< Storage for the x coordinate of the third atoms
  ExtArray<double, chunk, align> m_ry3; ///< Storage for the y coordinate of the third atoms
  ExtArray<double, chunk, align> m_rz3; ///< Storage for the z coordinate of the third atoms

  ExtArray<double, chunk, align> m_rx4; ///< Storage for the x coordinate of the forth atoms
  ExtArray<double, chunk, align> m_ry4; ///< Storage for the y coordinate of the forth atoms
  ExtArray<double, chunk, align> m_rz4; ///< Storage for the z coordinate of the forth atoms

  ExtArray<double, chunk, align> m_fx1; ///< Storage for the x components of the first forces
  ExtArray<double, chunk, align> m_fy1; ///< Storage for the y components of the first forces
  ExtArray<double, chunk, align> m_fz1; ///< Storage for the z components of the first forces

  ExtArray<double, chunk, align> m_fx2; ///< Storage for the x components of the second forces
  ExtArray<double, chunk, align> m_fy2; ///< Storage for the y components of the second forces
  ExtArray<double, chunk, align> m_fz2; ///< Storage for the z components of the second forces

  ExtArray<double, chunk, align> m_fx3; ///< Storage for the x components of the third forces
  ExtArray<double, chunk, align> m_fy3; ///< Storage for the y components of the third forces
  ExtArray<double, chunk, align> m_fz3; ///< Storage for the z components of the third forces

  ExtArray<double, chunk, align> m_ep; ///< Storage for the potential energy

  ExtArray<double, chunk, align> m_p1; ///< Storage for the first parameters
  ExtArray<double, chunk, align> m_p2; ///< Storage for the second parameters
  ExtArray<double, chunk, align> m_p3; ///< Storage for the third parameters
  ExtArray<double, chunk, align> m_p4; ///< Storage for the third parameters
  ExtArray<double*, chunk, align> m_pAccess; ///< Store pointers to the parameters

};

#endif /* VECTORIZATIONMBUFFER_HPP_ */
