/// @file 
/// @brief Definition of a basic templated three dimension vector

#ifndef __VEC3_HPP_INCLUDED
#define __VEC3_HPP_INCLUDED


#include <ostream>

#include "utils/auxMath.hpp"


/// @brief Number of dimensions in a vec3
#define VEC3_NDIMS 3


/// @brief A basic templated three dimension vector
/// @tparam T Class inside the vector
template<class T> struct vec3 {

  /// @brief Default constructor / constructor from a value
	/// @tparam T Class inside the vector
	/// @param [in] a Value for all the components, default=0
  vec3(const T& a = 0) : x(a), y(a), z(a) {}

  /// @brief Constructor from three values
  /// @tparam T Class inside the vector
  /// @param [in] a Value for the x component
  /// @param [in] b Value for the y component
  /// @param [in] c Value for the z component
  vec3(const T& a, const T& b, const T& c) : x(a), y(b), z(c) {}

  /// @brief Copy constructor
  /// @tparam T Class inside the vector
  /// @param [in] v Vector to copy
  vec3(const vec3<T> &v) : x(v.x), y(v.y), z(v.z) {}

  /// @brief Copy constructor
  /// @tparam T Class inside the vector
  /// @tparam U Class inside the vector to copy
  /// @param [in] v Vector to copy
  template <class U>
  vec3(const vec3<U>& v) : x(static_cast<T>(v.x)), y(static_cast<T>(v.y)), z(static_cast<T>(v.z)) {}

  /// @brief Destructor (nothing to do)
  ~vec3() {}

  T& operator [] (const uint8_t dim);
  const T& operator [] (const uint8_t dim) const;
  vec3<T>& operator = (const vec3<T>& v);

  bool operator == (const vec3<T>& v) const;
  bool operator != (const vec3<T>& v) const;

  vec3<T>& operator += (const vec3<T>& v);
  vec3<T>& operator -= (const vec3<T>& v);
  vec3<T>& operator *= (const vec3<T>& v);
  vec3<T>& operator /= (const vec3<T>& v);

  vec3<T>& operator += (const T& a);
  vec3<T>& operator -= (const T& a);
  vec3<T>& operator *= (const T& a);
  vec3<T>& operator /= (const T& a);

  vec3<T> operator + (const vec3<T>& v) const;
  vec3<T> operator - (const vec3<T>& v) const;
  vec3<T> operator * (const vec3<T>& v) const;
  vec3<T> operator / (const vec3<T>& v) const;

  vec3<T> operator + (const T& a) const;
  vec3<T> operator - (const T& a) const;
  vec3<T> operator * (const T& a) const;
  vec3<T> operator / (const T& a) const;

  vec3<bool> operator >  (const vec3<T>& v) const;
  vec3<bool> operator <  (const vec3<T>& v) const;
  vec3<bool> operator >= (const vec3<T>& v) const;
  vec3<bool> operator <= (const vec3<T>& v) const;
  
  vec3<bool> operator >  (const T& a) const;
  vec3<bool> operator <  (const T& a) const;
  vec3<bool> operator >= (const T& a) const;
  vec3<bool> operator <= (const T& a) const;

  template<class Tt>
  friend std::ostream& operator << (std::ostream& out, const vec3<Tt>& v);

  T x; ///< X component
  T y; ///< Y component
  T z; ///< Z component

};


// =============================================================================
// ==== ARITHMETIC OPERATORS ===================================================
// =============================================================================


template<class T> vec3<T> operator + (const T& a, const vec3<T>& u);
template<class T> vec3<T> operator - (const T& a, const vec3<T>& u);
template<class T> vec3<T> operator * (const T& a, const vec3<T>& u);
/// No reverse operator for division (don't divide a constant by a vector ...)


#include "utils/vec3/vec3_arithmeticOperators.hxx"


// =============================================================================
// ==== LOGICAL OPERATORS ======================================================
// =============================================================================


#include "vec3_logicalOperators.hxx"


/// @brief Logical AND operator
///
/// Do a AND for each pair of components
/// @param [in] a Boolean vector
/// @param [in] b Boolean vector
inline vec3<bool> operator && (const vec3<bool>& a, const vec3<bool>& b) { 
  return vec3<bool>(a.x && b.x, a.y && b.y, a.z && b.z);
}


/// @brief Logical OR operator
///
/// Do a OR for each pair of components
/// @param [in] a Boolean vector
/// @param [in] b Boolean vector
inline vec3<bool> operator || (const vec3<bool>& a, const vec3<bool>& b) { 
  return vec3<bool>(a.x || b.x, a.y || b.y, a.z || b.z);
}


// =============================================================================
// ==== SOME UTILITY FUNCTIONS 
// =============================================================================


template<class T> T dot(const vec3<T>& u, const vec3<T>& v);
template<class T> vec3<T> cross(const vec3<T>& u, const vec3<T>& v);
template<class T> T norm2(const vec3<T>& v);


template<class T> T sum(const vec3<T>& v);
template<class T> T product(const vec3<T>& v);


template<class T> int countNumberOf(T val, const vec3<T>& v);


#include "utils/vec3/vec3_utils.hxx"


// =============================================================================
// ==== AUX MATH PACKAGE =======================================================
// =============================================================================


#include "utils/vec3/vec3_auxMath.hxx"


// =============================================================================
// ==== IMPLEMENTATION OF OTHER MEMBERS ========================================
// =============================================================================


/// @brief Component access operator
/// @tparam T Class inside the vector
/// @param [in] dim Component to access
/// @return Reference to the component
template<class T> inline T& vec3<T>::operator [] (const uint8_t dim) {
    return *(reinterpret_cast<T*>(this) + dim);
}


/// @brief Constant component access operator
/// @tparam T Class inside the vector
/// @param [in] dim Component to access
/// @return Constant reference to the component
template<class T> inline const T& vec3<T>::operator [] (const uint8_t dim) const {
  return *(reinterpret_cast<const T*>(this) + dim);
}


/// @brief Assignment operator
/// @tparam T Class inside the vectors
/// @param [in] v Vector to copy
template<class T> inline vec3<T>& vec3<T>::operator = (const vec3<T>& v) {
  x = v.x;
  y = v.y;
  z = v.z;
  return *this;
}


/// @brief Stream insertion operator for a vector
/// @tparam Tt Class inside the vector
/// @param [in,out] out Stream
/// @param [in] v Vector to print to print
template<class Tt> inline std::ostream& operator << (std::ostream& out, const vec3<Tt>& v) {
  out<< v.x << " " << v.y << " " << v.z;
  return out;
}


/// @brief Get a vector with zeros
/// @return Vector with zeros
inline vec3<int> zeros() {
  return vec3<int>(0);
}


/// @brief Get a vector with ones
/// @return Vector with ones
inline vec3<int> ones() {
  return vec3<int>(1);
}

#endif // __VEC3_HPP_INCLUDED
