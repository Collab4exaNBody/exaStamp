/// @file
/// @brief Definition of class Stat

#ifndef __DIAG_HPP_INCLUDED
#define __DIAG_HPP_INCLUDED


#include <iostream>
#include <iomanip>
#include <string>

#include "utils/auxMath.hpp"


/// @brief Tool to print basic stats about a set of objects (not used)
/// @param [in] T Type of the set
template <class T> struct Stat {

	/// @brief Default constructor
  Stat() {}

  /// @brief Constructor from a set
	/// @param [in] data Data set
  /// @param [in] n Size of the data set
  Stat(const T* data, const uint n) : size(n) {
    quickAnalysis(data, n, mean, std, min, max);
  }

  /// @brief Destructor (nothing to do)
  ~Stat() {}

  uint size; ///< Size of the data set
  double mean ///< Mean
  double std; ///< Standart deviation
  T min; ///< Minimum
  T max; ///< Maximum

  /// @brief Print the stats
  /// @param [in] str Introductory comment
  void print(const std::string& str) {

    std::cout<< str
	     << "mean " << std::fixed << std::setprecision(1) << mean << " "
	     << "std "  << std::fixed << std::setprecision(1) << std  << " "
	     << "min "  << std::fixed << std::setprecision(1) << min  << " "
	     << "max "  << std::fixed << std::setprecision(1) << max  << " "
	     << std::endl;

  }

};


#endif // __DIAG_HPP_INCLUDED
