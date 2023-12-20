/// @file 
/// @brief Classes to map 1D to 3D indexes

#ifndef __CURVE_HPP_INCLUDED
#define __CURVE_HPP_INCLUDED


#include "utils/vec3/vec3.hpp"


/// @brief Class to map a rectangular box with a 1D index
/// @tparam T Type of the index and coordinates (must be an integer type)
template <class T=int>
class LCurve {

	/// @brief Shortcut the index type
  typedef T index_t;
  /// @brief Shortcut for the coordinates type
  typedef vec3<T> coord_t;

public:

  /// @brief Default constructor
  ///
	///
  LCurve() : m_size(0), m_yz(0) {}

  /// @brief Destructor (nothing to do)
  virtual ~LCurve() {}

  /// @brief Constructor from a size
	/// @param [in] size Size of the box to map
  LCurve(const coord_t& size) : m_size(size), m_yz(m_size.y*m_size.z) {}

  /// @brief Copy constructor (not used)
	/// @param [in] converter LCurve to copy
  LCurve(const LCurve<T>& converter) : m_size(converter.m_size), m_yz(m_size.y*m_size.z) {}

  /// @brief Operator = (not used)
  /// @param [in] converter Copied LCurve
  LCurve<T>& operator = (const LCurve<T>& converter) {
    set(converter.m_size);
    return *this;
  }

  /// @brief Set LCurve from a box size
  /// @brief size Box size
  void set(const coord_t& size) {
    m_size = size;
    m_yz   = m_size.y*m_size.z;
  }

  /// @brief Convert an index to coordinates
  /// @param [in] index Index to convert
  /// @return Coordinates
  inline coord_t convert(const index_t& index) const {
    return coord_t( index/m_yz, (index/m_size.z) % m_size.y, index % m_size.z );
  }
  
  /// @brief Convert coordinates to an index
  /// @param [in] coord Coordinates to convert
  /// @return Index
  inline index_t convert(const coord_t& coord) const {
    return coord.x*m_yz + coord.y*m_size.z + coord.z;
  }
  
private:
  
  coord_t m_size; ///< Size of the box to map
  index_t m_yz; ///< Product of y and z components of the size
  
};


/// @brief A complex bits operation used in ZCurve to convert coordinates to indexes
/// @param [in] v_ Integer before
/// @return Integer after
static inline uint64_t spreadBits(uint64_t v_) {
  uint64_t v = v_;
  uint64_t m1 = 0x30c30c30c30c30c3UL;
  uint64_t m2 = 0x9249249249249249UL;
  v = (v | (v << 32)) & 0xffff00000000ffff; 
  v = (v | (v << 16)) & 0x00ff0000ff0000ff;
  v = (v | (v <<  8)) & 0xf00f00f00f00f00f; 
  v = (v | (v <<  4)) & m1;
  v = (v | (v <<  2)) & m2;
  return v;
}


/// @brief A complex bits operation used in ZCurve to convert indexes to coordinates
/// @param [in] v_ Integer before
/// @return Integer aft
static inline uint64_t unspreadBits(uint64_t v_) {
  uint64_t v = v_;

  v = (v )            & 0x9249249249249249UL;
  v = (v ^ (v >>  2)) & 0x30c30c30c30c30c3UL;
  v = (v ^ (v >>  4)) & 0xf00f00f00f00f00fUL; 
  v = (v ^ (v >>  8)) & 0x00ff0000ff0000ffUL;
  v = (v ^ (v >> 16)) & 0xffff00000000ffffUL; 
  v = (v ^ (v >> 32)) & 0xffffffffUL; 

  return v;
}


/// @brief Another class to convert indexes to coordinates (not used)
class ZCurve {

public:

	/// @brief Default constructor
  ZCurve() {}
  /// @brief Destructor (nothing to do)
  ~ZCurve() {}

  /// @brief Convert an index to coordinates
  /// @param [in] index Index to convert
  /// @return Coordinates
  inline vec3<uint64_t> convert(uint64_t index) {

    vec3<uint64_t> coords;
    
    uint64_t m = index;
    coords.x = unspreadBits(m);
    m = index >> 1;
    coords.y = unspreadBits(m);
    m = index >> 2;
    coords.z = unspreadBits(m);

    return coords;

  }

  /// @brief Convert coordinates to index
  /// @param [in] coords Coordinates to convert
  /// @return Index
  inline uint64_t convert(const vec3<uint64_t>& coords) {

    uint64_t r1 = spreadBits(coords.x);
    uint64_t r2 = spreadBits(coords.y);
    uint64_t r3 = r2 << 1;
    uint64_t r4 = r1 | r3;
    uint64_t r5 = spreadBits(coords.z);
    uint64_t r6 = r5 << 2;
    uint64_t r7 = r4 | r6;
    
    return r7;

  }

};

#endif // __CURVE_HPP_INCLUDED
