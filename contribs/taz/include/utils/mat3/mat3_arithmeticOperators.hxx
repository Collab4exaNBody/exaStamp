/// @file
/// @brief Arithmetic operators for mat3
/// (and specialization for int/double cases)


/// @brief Operator += with another mat3
/// @tparam T Class inside the matrices
/// @param [in] m Mat3 to add
template<class T> mat3<T>& mat3<T>::operator += (const mat3<T>& m) {
  m11 += m.m11;
  m12 += m.m12;
  m13 += m.m13;
  m21 += m.m21;
  m22 += m.m22;
  m23 += m.m23;
  m31 += m.m31;
  m32 += m.m32;
  m33 += m.m33;
  return *this;
}


/// @brief Operator -= with another mat3
/// @tparam T Class inside the matrices
/// @param [in] m Mat3 to subtract
template<class T> mat3<T>& mat3<T>::operator -= (const mat3<T>& m) {
  m11 -= m.m11;
  m12 -= m.m12;
  m13 -= m.m13;
  m21 -= m.m21;
  m22 -= m.m22;
  m23 -= m.m23;
  m31 -= m.m31;
  m32 -= m.m32;
  m33 -= m.m33;
  return *this;
}


/// @brief Operator *= with another mat3
///
/// Component to component multiplication
/// @tparam T Class inside the matrices
/// @param [in] m Mat3 to multiply
template<class T> mat3<T>& mat3<T>::operator *= (const mat3<T>& m) {
  m11 *= m.m11;
  m12 *= m.m12;
  m13 *= m.m13;
  m21 *= m.m21;
  m22 *= m.m22;
  m23 *= m.m23;
  m31 *= m.m31;
  m32 *= m.m32;
  m33 *= m.m33;
  return *this;
}


/// @brief Operator /= with another mat3
///
/// Component to component division
/// @tparam T Class inside the matrices
/// @param [in] m Mat3 to divide by
template<class T> mat3<T>& mat3<T>::operator /= (const mat3<T>& m) {
  m11 /= m.m11;
  m12 /= m.m12;
  m13 /= m.m13;
  m21 /= m.m21;
  m22 /= m.m22;
  m23 /= m.m23;
  m31 /= m.m31;
  m32 /= m.m32;
  m33 /= m.m33; 
  return *this;
}


/// @brief Operator += with a constant
///
/// Add a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to add
template<class T> mat3<T>& mat3<T>::operator += (const T& a) {
  return *this+=mat3<T>(a);
}


/// @brief Operator -= with a constant
///
/// Subtract a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to subtract
template<class T> mat3<T>& mat3<T>::operator -= (const T& a) {
  return *this-=mat3<T>(a);
}


///  @brief Operator *= with a constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to multiply
template<class T> mat3<T>& mat3<T>::operator *= (const T& a) {
  return *this*=mat3<T>(a);
}


///  @brief Operator /= with a constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to divide by
template<class T> mat3<T>& mat3<T>::operator /= (const T& a) {
  return *this/=mat3<T>(a);
}


/// @brief Sum of two mat3
/// @tparam T Class inside the matrices
/// @param [in] m Mat3 to add
template<class T> mat3<T> mat3<T>::operator + (const mat3<T>& m) const {
  return mat3<T>(*this)+=m;
}


/// @brief Difference of two mat3
/// @tparam T Class inside the matrices
/// @param [in] m Mat3 to subtract
template<class T> mat3<T> mat3<T>::operator - (const mat3<T>& m) const {
  return mat3<T>(*this)-=m;
}


/// @brief Member to member multiplication
///
/// Component to component multiplication
/// @tparam T Class inside the matrices
/// @param [in] m Mat3 to multiply
template<class T> mat3<T> mat3<T>::operator * (const mat3<T>& m) const {
  return mat3<T>(*this)*=m;
}


/// @brief Member to member division
///
/// Component to component division
/// @tparam T Class inside the matrices
/// @param [in] m Mat3 to divide by
template<class T> mat3<T> mat3<T>::operator / (const mat3<T>& m) const {
  return mat3<T>(*this)/=m;
}


/// @brief Multiplication with a vector
/// @tparam T Class inside
/// @param [in] u Vector to multiply
template<class T> vec3<T> mat3<T>::operator * (const vec3<T>& u) const {
  return vec3<T>(m11*u.x + m12*u.y + m13*u.z, m21*u.x + m22*u.y + m23*u.z, m31*u.x + m32*u.y + m33*u.z);
}


/// @brief Sum with a constant
///
/// Add a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to add
template<class T> mat3<T> mat3<T>::operator + (const T& a) const {
  return mat3<T>(*this)+=a;
}


/// @brief Difference with a constant
///
/// Subtract a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to subtract
template<class T> mat3<T> mat3<T>::operator - (const T& a) const {
  return mat3<T>(*this)-=a;
}


/// @brief Multiplication with a constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to multiply
template<class T> mat3<T> mat3<T>::operator * (const T& a) const {
  return mat3<T>(*this)*=a;
}


/// @brief Division with a constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to divide by
template<class T> mat3<T> mat3<T>::operator / (const T& a) const {
  return mat3<T>(*this)/=a;
}



/// @brief Sum with a constant (reverse)
///
/// Add a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to add
/// @param [in] u Mat3 to add
template<class T> mat3<T> operator + (const T& a, const mat3<T>& u) {
  return u+a;
}


/// @brief Difference with a constant (reverse)
///
/// Subtract to a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant
/// @param [in] u Mat3 to subtract
template<class T> mat3<T> operator - (const T& a, const mat3<T>& u) {
  return u-a;
}


/// @brief Multiplication with a constant (reverse)
/// @tparam T Class inside the matrix
/// @param [in] a Constant to multiply
/// @param [in] u Mat3 to multiply
template<class T> mat3<T> operator * (const T& a, const mat3<T>& u) {
  return u*a;
}


// =============================================================================
// ==== SOME SPECIALIZATIONS MAT3<T> ===========================================
// =============================================================================


/// @brief Sum of an integer and a double matrix
/// @param [in] c Integer constant
/// @param [in] m Matrix of doubles to add
/// @return Sum as a matrix of doubles
inline mat3<double> operator + (const int c, const mat3<double> &m) {
  return mat3<double>( c+m.m11, c+m.m12, c+m.m13, c+m.m21, c+m.m22, c+m.m23, c+m.m31, c+m.m32, c+m.m33 );
}


/// @brief Sum of a double and an integer matrix
/// @param [in] c Double constant
/// @param [in] m Matrix of integers to add
/// @return Sum as a matrix of doubles
inline mat3<double> operator + (const double c, const mat3<int> &m) {
  return mat3<double>( c+m.m11, c+m.m12, c+m.m13, c+m.m21, c+m.m22, c+m.m23, c+m.m31, c+m.m32, c+m.m33 );
}


/// @brief Sum of an integer and a double matrix (reverse)
/// @param [in] u Matrix of integers
/// @param [in] c Double constant to add
/// @return Sum as a matrix of doubles
inline mat3<double> operator + (const mat3<int> &u, const double c) {
  return c+u;
}


/// @brief Sum of a double and an integer matrix (reverse)
/// @param [in] u Matrix of doubles
/// @param [in] c Integer constant to add
/// @return Sum as a matrix of doubles
inline mat3<double> operator + (const mat3<double> &u, const int c) {
  return c+u;
}


/// @brief Sum between integer and double matrices
/// @param [in] u Matrix of integers
/// @param [in] m Matrix of doubles
/// @return Sum as a matrix of doubles
inline mat3<double> operator + (const mat3<int> &u, const mat3<double> &m) {
  return mat3<double>( u.m11+m.m11, u.m12+m.m12, u.m13+m.m13, u.m21+m.m21, u.m22+m.m22, u.m23+m.m23, u.m31+m.m31, u.m32+m.m32, u.m33+m.m33 );
}


/// @brief Sum between integer and double matrices (reverse)
/// @param [in] u Matrix of doubles
/// @param [in] m Matrix of integers
/// @return Sum as a matrix of doubles
inline mat3<double> operator + (const mat3<double> &u, const mat3<int> &m) {
  return m+u;
}


/// @brief Difference between an integer and a double matrix
/// @param [in] c Integer constant
/// @param [in] m Matrix of doubles to subtract
/// @return Difference as a matrix of doubles
inline mat3<double> operator - (const int c, const mat3<double> &m) {
  return mat3<double>( c-m.m11, c-m.m12, c-m.m13, c-m.m21, c-m.m22, c-m.m23, c-m.m31, c-m.m32, c-m.m33 );
}


/// @brief Difference between a double and an integer matrix
/// @param [in] c Double constant
/// @param [in] m Matrix of integers to subtract
/// @return Difference as a matrix of doubles
inline mat3<double> operator - (const double c, const mat3<int> &m) {
  return mat3<double>( c-m.m11, c-m.m12, c-m.m13, c-m.m21, c-m.m22, c-m.m23, c-m.m31, c-m.m32, c-m.m33 );
}


/// @brief Difference between an integer and a double matrix (reverse)
/// @param [in] u Matrix of integers
/// @param [in] c Double constant to subtract
/// @return Difference as a matrix of doubles
inline mat3<double> operator - (const mat3<int> &u, const double c) {
  return mat3<double>( u.m11-c, u.m12-c, u.m13-c, u.m21-c, u.m22-c, u.m23-c, u.m31-c, u.m32-c, u.m33-c );
}


/// @brief Difference betwwen a double and an integer matrix (reverse)
/// @param [in] u Matrix of doubles
/// @param [in] c Integer constant to subtract
/// @return Difference as a matrix of doubles
inline mat3<double> operator - (const mat3<double> &u, const int c) {
  return mat3<double>( u.m11-c, u.m12-c, u.m13-c, u.m21-c, u.m22-c, u.m23-c, u.m31-c, u.m32-c, u.m33-c );
}


/// @brief Difference between integer and double matrices
/// @param [in] u Matrix of integers
/// @param [in] m Matrix of doubles
/// @return Difference as a matrix of doubles
inline mat3<double> operator - (const mat3<int> &u, const mat3<double> &m) {
  return mat3<double>( u.m11-m.m11, u.m12-m.m12, u.m13-m.m13, u.m21-m.m21, u.m22-m.m22, u.m23-m.m23, u.m31-m.m31, u.m32-m.m32, u.m33-m.m33 );
}


/// @brief Difference between integer and double matrices (reverse)
/// @param [in] u Matrix of doubles
/// @param [in] m Matrix of integers
/// @return Difference as a matrix of doubles
inline mat3<double> operator - (const mat3<double> &u, const mat3<int> &m) {
  return mat3<double>( u.m11-m.m11, u.m12-m.m12, u.m13-m.m13, u.m21-m.m21, u.m22-m.m22, u.m23-m.m23, u.m31-m.m31, u.m32-m.m32, u.m33-m.m33 );
}


/// @brief Product between an integer and a double matrix
/// @param [in] c Integer constant
/// @param [in] m Matrix of doubles
/// @return Product as a matrix of doubles
inline mat3<double> operator * (const int c, const mat3<double> &m) {
  return mat3<double>( c*m.m11, c*m.m12, c*m.m13, c*m.m21, c*m.m22, c*m.m23, c*m.m31, c*m.m32, c*m.m33 );
}


/// @brief Product between a double and an integer matrix
/// @param [in] c Double constant
/// @param [in] m Matrix of integers
/// @return Product as a matrix of doubles
inline mat3<double> operator * (const double c, const mat3<int> &m) {
  return mat3<double>( c*m.m11, c*m.m12, c*m.m13, c*m.m21, c*m.m22, c*m.m23, c*m.m31, c*m.m32, c*m.m33 );
}


/// @brief Product between an integer and a double matrix (reverse)
/// @param [in] u Matrix of integers
/// @param [in] c Double constant
/// @return Product as a matrix of doubles
inline mat3<double> operator * (const mat3<int> &u, const double c) {
  return c*u;
}


/// @brief Product betwwen a double and an integer matrix (reverse)
/// @param [in] u Matrix of doubles
/// @param [in] c Integer constant
/// @return Product as a matrix of doubles
inline mat3<double> operator * (const mat3<double> &u, const int c) {
  return c*u;
}


/// @brief Product between integer and double matrices
///
/// Component to component multiplication
/// @param [in] u Matrix of doubles
/// @param [in] m Matrix of integers
/// @return Product as a matrix of doubles
inline mat3<double> operator * (const mat3<double> &u, const mat3<int> &m) {
  return mat3<double>( u.m11*m.m11, u.m12*m.m12, u.m13*m.m13, u.m21*m.m21, u.m22*m.m22, u.m23*m.m23, u.m31*m.m31, u.m32*m.m32, u.m33*m.m33 );
}


/// @brief Product between integer and double matrices (reverse)
///
/// Component to component multiplication
/// @param [in] u Matrix of integers
/// @param [in] m Matrix of doubles
/// @return Product as a matrix of doubles
inline mat3<double> operator * (const mat3<int> &u, const mat3<double> &m) {
  return m*u;
}


// =============================================================================
// ==== DIVISION OPERATOR (ONLY A FEW) =========================================
// =============================================================================


/// @brief Division of an integer matrix by a double
/// @param [in] u Integer matrix
/// @param [in] c Double constant
/// @return Result as a matrix of doubles
inline mat3<double> operator / (const mat3<int>& u, const double c) {
  return (1.*u)/c;
}


/// @brief Division of a double matrix by an integer
/// @param [in] u Double matrix
/// @param [in] c Integer constant
/// @return Result as a matrix of doubles
inline mat3<double> operator / (const mat3<double>& u, const int c) {
  return u/(1.*c);
}


/// @brief Division of an integer matrix by a double matrix
///
/// Component to component division
/// @param [in] u Integer matrix
/// @param [in] m Double matrix
/// @return Result as a matrix of doubles
inline mat3<double> operator / (const mat3<int>& u, const mat3<double>& m) {
  return (1.*u)/m;
}


/// @brief Division of a double matrix by an integer matrix
///
/// Component to component division
/// @param [in] u Double matrix
/// @param [in] m Integer matrix
/// @return Result as a matrix of doubles
inline mat3<double> operator / (const mat3<double>& u, const mat3<int>& m) {
  return u/(1.*m);
}
