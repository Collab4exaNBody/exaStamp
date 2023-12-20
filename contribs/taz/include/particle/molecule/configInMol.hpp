/*
 * configInMol.hpp
 *
 *  Created on: Sep 16, 2016
 *      Author: giarda
 */

///@file
/// @brief Specialization of the Configuration for the molecules data


#ifndef CONFIGINMOL_HPP_
#define CONFIGINMOL_HPP_


#include <vector>
#include <string>

#include "forceField/forceField.hpp"

#include "io/input.hpp"

#include "parallel/types/MPI_particle.hpp"


template <class T> class Configuration;
class MPI__InMol;
class MPI__Particle;


/// @brief Temporary structure gathering all data to initialize the molecules
template <> class Configuration<MPI__InMol> {

public :

	Configuration();
	Configuration(const Input& input);
	/// @brief Destructor (nothing to do)
	~Configuration() {}
	inline void clear();

	// Init file data
	FF::Type m_forceField; ///< Type of the force field (add to init file)
	FF::ChargeMethod m_charges; ///< Method to set the partial charges (add to init file)
	std::string m_fileName; ///< Name of the molecular data file (allready in the domain config)
	std::string m_fileDir; ///< Directory of the molecular data file (allready in the domain config)
	std::string m_format; ///< Format of the molecular data file (add to init file and get a default)
	bool m_add; ///< Flag indicating if hydrogens must be added to the molecular data file (add to init file and make a default)
	int m_conformer; ///< Number of the conformer to consider in the molecular data file (add to init file and make a default)
  vec3<double> m_margins; ///< Empty border around the molecule (on the upper side only)
  double m_initialTemperature; ///< Initial temperature of the system, used to initialize the velocities
  double m_maxBond; ///< Maximum bond length
  vec3<int> m_nCells; ///< Number of cells of the lattice in each dimension
  // Molecular input data
  uint64_t m_numMol; ///< Number of molecules in the base cell
  vec3<double> m_cellSize; ///< Size of the base cell
  std::vector<MPI__Particle*> m_particles; ///< Point particles of the molecules in the base cell
  std::vector<std::string> m_listOfSubtypes; ///< List of the all subtypes of the atoms


};


/// @brief Clear the particle vectors
///
///
void Configuration<MPI__InMol>::clear() {
	for(uint i(0); i<m_particles.size();++i) {
		delete m_particles[i];
	}
	m_particles.clear();
}

#endif /* CONFIGINMOL_HPP_ */
