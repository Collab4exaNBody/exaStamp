/// @file 
/// @brief Definition of Zone and ZoneTwo classes, used to build the traversals

#ifndef __ZONE_HPP_INCLUDED
#define __ZONE_HPP_INCLUDED


#include "utils/array/array.hpp"
#include "utils/vec3/vec3.hpp"


/// @brief Class ZoneTow : defines a rectangular box in the cells grid
class ZoneTwo {

public:

	/// @brief Default constructor
	///
	/// Fixe limits size at 2 (upper limit and lower limit)
  ZoneTwo() : limits(Array< vec3<int> >(2)) {}
  /// @brief Constructor from limits
	/// @param [in] a Lower limit
  /// @param [in] b Upper limit
  ZoneTwo(const vec3<int>& a, const vec3<int>& b)
    : limits(Array< vec3<int> >(2)) {
    limits[0] = a;
    limits[1] = b;
  }

  /// @brief Destructor (nothing to do)
  ~ZoneTwo() {}

  /// @brief Calculate the size of the rectangular box
  /// @return Size of the zone
  uint estimatedSize() {
    vec3<int> size;
    for (uint dim=0; dim<VEC3_NDIMS; ++dim) size[dim]  = limits[1][dim] - limits[0][dim];
    return (uint) product(size);
  }

  /// @brief Find if specified cell is in the zone
  /// @param [in] coords Cell coordinates
  /// @return True if the cell is in the zone
  bool isIn(const vec3<int>& coords) {
    for (uint dim=0; dim<VEC3_NDIMS; ++dim) {
      if ( coords[dim]<limits[0][dim] || coords[dim]>=limits[1][dim] ) return false;
    }
    return true;
  }

private: 
  
  Array< vec3<int> > limits; ///< Limits of the box

};


/// @brief Defines a rectangular box from which a smaller rectangular box has been extracted
class Zone {

public:

	/// @brief Default constructor
	///
	///
  Zone() : big(), small() {}

  /// @brief Constructor
	/// @param [in] a Lower limit of the external box
  /// @param [in] b Lower limit of the internal box
  /// @param [in] c Upper limit of the internal box
  /// @param [in] d Upper limit of the external box
  Zone(const vec3<int>& a, const vec3<int>& b, const vec3<int>& c, const vec3<int>& d)
    : big(a, d), small(b, c) {}

  /// @brief Destructor (nothing to do)
  ~Zone() {}

  /// @brief Calculate the size of the zone
  /// @return Size of the zone
  int estimatedSize() {
    return big.estimatedSize()-small.estimatedSize();
  }

  /// @brief Find if specified cell is in the zone
  /// @param [in] coords Cell coordinates
  /// @return True if the cell is in the zone
  bool isIn(const vec3<int>& coords) {
    return big.isIn(coords) && (!small.isIn(coords));
  }

private:

  ZoneTwo big; ///< External box
  ZoneTwo small; ///< Internal box

};

#endif //  __ZONE_HPP_INCLUDED
