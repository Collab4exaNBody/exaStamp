/// @file 
/// @brief Definition and implementation of distance modulator : tool to modulate a distance and set it back between to bounds

#ifndef __MODULER_HPP_INCLUDED
#define __MODULER_HPP_INCLUDED


#include "utils/vec3/vec3.hpp"


/// @brief Distance modulator : tool to modulate a distance and set it back between two bounds
/// @tparam T Type of the distance to modulate (double or int)
template <class T> 
class Moduler {

public:

  /// @brief Default constructor
  Moduler() {}

  /// @brief Constructor with data (not used)
  /// @tparam T Type of the distance to modulate (double or int)
  /// @param [in] setup Dimensions in which distance modulation is necessary
  /// @param [in] minBounds Lower bounds in each dimension
  /// @param [in] maxBounds Upper bounds in each dimension
  Moduler(const vec3<bool>& setup, const vec3<T>& minBounds, const vec3<T>& maxBounds)
    : dim(setup), min(minBounds), max(maxBounds), dif(max-min) {
  }

  /// @brief Destructor (nothing to do)
  ~Moduler() {}

  /// @brief Set Moduler from parameters
  /// @tparam T Type of the distance to modulate (double or int)
  /// @param [in] setup Dimensions in which distance modulation is necessary
  /// @param [in] minBounds Lower bound in each dimension
  /// @param [in] maxBounds Upper bound in each dimension
  void set(const vec3<bool>& setup, const vec3<T>& minBounds, const vec3<T>& maxBounds) {
    dim = setup;
    min = minBounds;
    max = maxBounds;
    dif = max-min;
  }

  /// @brief Modulate a 3 dimensional distance
  /// @tparam T Type of the distance to modulate (double or int)
  /// @param [in,out] distance 3D distance to modulate
  inline void module(vec3<T>& distance) const {
    if (dim.x) module<0>(distance.x);
    if (dim.y) module<1>(distance.y);
    if (dim.z) module<2>(distance.z);
  }

  /// @brief Modulate 3 one dimensional distances
  /// @tparam T Type of the distance to modulate (double or int)
  /// @param [in,out] dx X component to modulate
  /// @param [in,out] dy Y component to modulate
  /// @param [in,out] dz Z component to modulate
  inline void module(T& dx, T& dy, T& dz) const {
    if (dim.x) module<0>(dx);
    if (dim.y) module<1>(dy);
    if (dim.z) module<2>(dz);
  }

  /// @brief Modulate a component to a distance
  /// @tparam T Type of the distance to modulate (double or int)
  /// @tparam DIM Identifier of the dimension (1 for x, 2 for y and 3 for z)
  /// @param [in,out] distance Component to modulate
  template <uint8_t DIM> inline void module(T& distance) const {
    if      (distance <  min[DIM]) distance += dif[DIM];
    else if (distance >= max[DIM]) distance -= dif[DIM];
  }

private:

  vec3<bool> dim; ///< Identify dimensions in which distances must be modulated

  vec3<T> min; ///< Lower bound in each dimension
  vec3<T> max; ///< Upper bound in each dimension
  vec3<T> dif; ///< Gap between the bounds in each dimension

};

#endif // __MODULER_HPP_INCLUDED
