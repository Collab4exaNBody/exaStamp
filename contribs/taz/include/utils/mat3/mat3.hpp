/// @file
/// @brief A basic templated 3x3 matrix structure

#ifndef MAT3_HPP
#define MAT3_HPP


#include <iostream>
#include <string>

#include "utils/vec3/vec3.hpp"


/// @brief Number of dimensions in a mat3
#define MAT3_NDIMS 3


/// @brief A basic templated three dimension matrix
/// @tparam T Class inside the matrix
template<class T> struct mat3 {

  /// @brief Default constructor / constructor from a value
	/// @tparam T Class inside the matrix
	/// @param [in] a Value for all the components, default=0
  mat3(T a = 0) : m11(a),  m12(a), m13(a), m21(a), m22(a), m23(a), m31(a), m32(a), m33(a){}

  /// @brief Constructor with nine values
  /// @tparam T Class inside the matrix
  /// @param [in] a Value for the 11 component
  /// @param [in] b Value for the 12 component
  /// @param [in] c Value for the 13 component
  /// @param [in] d Value for the 21 component
  /// @param [in] e Value for the 22 component
  /// @param [in] f Value for the 23 component
  /// @param [in] g Value for the 31 component
  /// @param [in] h Value for the 32 component
  /// @param [in] i Value for the 33 component
  mat3(T a, T b, T c, T d, T e, T f,T g, T h, T i) : m11(a),  m12(b), m13(c), m21(d), m22(e), m23(f), m31(g), m32(h), m33(i) {}

  /// @brief Copy constructor
  /// @tparam T Class inside the matrix
  /// @param [in] m Matrix to copy
  mat3(const mat3<T> &m) : m11(m.m11),  m12(m.m12), m13(m.m13), m21(m.m21), m22(m.m22), m23(m.m23), m31(m.m31), m32(m.m32), m33(m.m33) {}

  /// @brief Destructor (nothing to do)
  ~mat3() {}


  T& operator () (int dim1, int dim2);
  const T& operator () (int dim1, int dim2) const;
  mat3<T>& operator = (const mat3<T>& m);

  bool operator == (const mat3<T>& m) const;
  bool operator != (const mat3<T>& m) const;


  mat3<T>& operator += (const mat3<T>& m);
  mat3<T>& operator -= (const mat3<T>& m);
  mat3<T>& operator *= (const mat3<T>& m);
  mat3<T>& operator /= (const mat3<T>& m);

  mat3<T>& operator += (const T& a);
  mat3<T>& operator -= (const T& a);
  mat3<T>& operator *= (const T& a);
  mat3<T>& operator /= (const T& a);


  mat3<T> operator + (const mat3<T>& m) const;
  mat3<T> operator - (const mat3<T>& m) const;
  mat3<T> operator * (const mat3<T>& m) const;
  mat3<T> operator / (const mat3<T>& m) const;

  vec3<T> operator * (const vec3<T>& u) const;

  mat3<T> operator + (const T& a) const;
  mat3<T> operator - (const T& a) const;
  mat3<T> operator * (const T& a) const;
  mat3<T> operator / (const T& a) const;



  mat3<bool> operator >  (const mat3<T>& m) const;
  mat3<bool> operator <  (const mat3<T>& m) const;
  mat3<bool> operator >= (const mat3<T>& m) const;
  mat3<bool> operator <= (const mat3<T>& m) const;
  
  mat3<bool> operator >  (const T& a) const;
  mat3<bool> operator <  (const T& a) const;
  mat3<bool> operator >= (const T& a) const;
  mat3<bool> operator <= (const T& a) const;


  template<class Tt>
  friend std::ostream& operator << (std::ostream& out, const mat3<Tt>& m);

  T m11; ///< 11 matrix component
  T m12; ///< 12 matrix component
  T m13; ///< 12 matrix component
  T m21; ///< 21 matrix component
  T m22; ///< 22 matrix component
  T m23; ///< 23 matrix component
  T m31; ///< 31 matrix component
  T m32; ///< 32 matrix component
  T m33; ///< 33 matrix component
};



// =============================================================================
// ==== ARITHMETIC OPERATORS ===================================================
// =============================================================================

template<class T> mat3<T> operator + (const T& a, const mat3<T>& m);
template<class T> mat3<T> operator - (const T& a, const mat3<T>& m);
template<class T> mat3<T> operator * (const T& a, const mat3<T>& m);
/// No reverse operator for division (don't divide a constant by a matrice ...)

#include "mat3_arithmeticOperators.hxx"



// =============================================================================
// ==== LOGICAL OPERATORS ======================================================
// =============================================================================

#include "mat3_logicalOperators.hxx"



// =============================================================================
// ==== SOME UTILITY FUNCTIONS 
// =============================================================================


template<class T> T dot(const mat3<T>& m, const mat3<T>& n);
template<class T> vec3<T> dot(const mat3<T>& m, const vec3<T>& v);
template<class T> mat3<T> multiply(const mat3<T>& m, const mat3<T>& n);
template<class T> T norm2(const mat3<T>& m);
template<class T> mat3<T> tensor(const vec3<T>& u, const vec3<T>& v);

template<class T> T sum(const mat3<T>& m);
template<class T> T product(const mat3<T>& m);
template<class T> mat3<T> transpose(const mat3<T>& m);


template<class T> void print(const mat3<T>& m, std::string separator=" ");
template<class T> int countNumberOf(T val, const mat3<T>& m);

#include "mat3_utils.hxx"



// =============================================================================
// ==== AUX MATH PACKAGE =======================================================
// =============================================================================

#include "mat3_auxMath.hxx"



// =============================================================================
// ==== IMPLEMENTATION OF OTHER MEMBERS ========================================
// =============================================================================



/// @brief Access with () operator
/// @tparam T Class inside the matrix
/// @param [in] dim1 Line to access
/// @param [in] dim2 Column to access
/// @return Component
template<class T> T& mat3<T>::operator () (int dim1, int dim2) {
  switch (dim1) {
  case 0 : switch (dim2) {
    case 0 : return m11;
    case 1 : return m12;
    case 2 : return m13;
    }
  case 1 :   switch (dim2) {
    case 0 : return m21;
    case 1 : return m22;
    case 2 : return m23;
    }
  case 2 :   switch (dim2) {
    case 0 : return m31;
    case 1 : return m32;
    case 2 : return m33;
    }
  }
  std::cout<< "WARNING : in mat3<T>::operator() (int,int) : "
	   << "accessing a component with strange index value" << std::endl;
  return (*this) (dim1>0 ? dim1%MAT3_NDIMS : dim1+MAT3_NDIMS, dim2>0 ? dim2%MAT3_NDIMS : dim2+MAT3_NDIMS);
}


/// @brief Constant access with () operator
/// @tparam T Class inside the matrix
/// @param [in] dim1 Line to access
/// @param [in] dim2 Column to access
/// @return Component
template<class T> const T& mat3<T>::operator () (int dim1, int dim2) const {
  switch (dim1) {
  case 0 : switch (dim2) {
    case 0 : return m11;
    case 1 : return m12;
    case 2 : return m13;
    }
  case 1 :   switch (dim2) {
    case 0 : return m21;
    case 1 : return m22;
    case 2 : return m23;
    }
  case 2 :   switch (dim2) {
    case 0 : return m31;
    case 1 : return m32;
    case 2 : return m33;
    }
  }
  std::cout<< "WARNING : in mat3<T>::operator() (int,int) : "
	   << "accessing a component with strange index value" << std::endl;
  return (*this) (dim1>0 ? dim1%MAT3_NDIMS : dim1+MAT3_NDIMS, dim2>0 ? dim2%MAT3_NDIMS : dim2+MAT3_NDIMS);
}



/// @brief Assignment operator
/// @tparam T Class inside the matrix
/// @param [in] m Matrix to copy
template<class T> mat3<T>& mat3<T>::operator = (const mat3<T>& m) {
  m11 = m.m11;
  m12 = m.m12;
  m13 = m.m13;
  m21 = m.m21;
  m22 = m.m22;
  m23 = m.m23;
  m31 = m.m31;
  m32 = m.m32;
  m33 = m.m33;
  return *this;
}



/// @brief Output stream operator
/// @param [in,out] out Output stream
/// @param [in] m Matrix
template<class Tt> std::ostream& operator << (std::ostream& out, const mat3<Tt>& m) {
  out<< m.m11 << " " << m.m12 << " " << m.m13 << std::endl
     << m.m21 << " " << m.m22 << " " << m.m23 << std::endl
     << m.m31 << " " << m.m32 << " " << m.m33 ;
  return out;
}

/// @brief Create a matrix full of 0
inline mat3<int> zerosMat() {
  return mat3<int>(0,0,0,0,0,0,0,0,0);
}

/// @brief Create a matrix full of 1
inline mat3<int> onesMat() {
  return mat3<int>(1,1,1,1,1,1,1,1,1);
}

/// @brief Create a diagonal matrix
/// @param [in] val Value on the diagonal
template<class T> inline mat3<T> diag(const T& val) {
  return mat3<T>(val,0,0,0,val,0,0,0,val);
}


#endif // MAT3_HPP
