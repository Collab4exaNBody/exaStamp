/// @file 
/// @brief Arithmetic operators for vec3
/// (and specialization for int/double cases)

/// @brief Sum assignment operator with another vector
/// @tparam T Class inside the vectors
/// @param [in] v Vector to add
template<class T> inline vec3<T>& vec3<T>::operator += (const vec3<T>& v) {
  x+=v.x;
  y+=v.y;
  z+=v.z;
  return *this;
}


/// @brief Difference assignment operator with another vector
/// @tparam T Class inside the vectors
/// @param [in] v Vector to subtract
template<class T> inline vec3<T>& vec3<T>::operator -= (const vec3<T>& v) {
  x-=v.x;
  y-=v.y;
  z-=v.z;
  return *this;
}


/// @brief Product assignment operator with another vector
/// @tparam T Class inside the vectors
/// @param [in] v Vector to multiply
template<class T> inline vec3<T>& vec3<T>::operator *= (const vec3<T>& v) {
  x*=v.x;
  y*=v.y;
  z*=v.z;
  return *this;
}


/// @brief Quotient assignment operator with another vector
/// @tparam T Class inside the vectors
/// @param [in] v Vector to divide by
template<class T> inline vec3<T>& vec3<T>::operator /= (const vec3<T>& v) {
  x/=v.x;
  y/=v.y;
  z/=v.z;
  return *this;
}


/// @brief Sum assignment operator with a constant
///
/// Add a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant to add
template<class T> inline vec3<T>& vec3<T>::operator += (const T& a) {
  return *this+=vec3<T>(a);
}


/// @brief Difference assignment operator with a constant
///
/// Subtract a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant to subtract
template<class T> inline vec3<T>& vec3<T>::operator -= (const T& a) {
  return *this-=vec3<T>(a);
}


/// @brief Product assignment operator with a constant
/// @tparam T Class inside the vector
/// @param [in] a Constant to multiply
template<class T> inline vec3<T>& vec3<T>::operator *= (const T& a) {
  return *this*=vec3<T>(a);
}


/// @brief Quotient assignment operator with a constant
/// @tparam T Class inside the vector
/// @param [in] a Constant to divide by
template<class T> inline vec3<T>& vec3<T>::operator /= (const T& a) {
  return *this/=vec3<T>(a);
}


/// @brief Sum operator for two vectors
/// @tparam T Class inside the vectors
/// @param [in] v Vector to add
template<class T> inline vec3<T> vec3<T>::operator + (const vec3<T>& v) const {
  return vec3<T>(*this)+=v;
}


/// @brief Difference operator for two vectors
/// @tparam T Class inside the vectors
/// @param [in] v Vector to subtract
template<class T> inline vec3<T> vec3<T>::operator - (const vec3<T>& v) const {
  return vec3<T>(*this)-=v;
}


/// @brief Member to member multiplication operator for two vectors
/// @tparam T Class inside the vectors
/// @param [in] v Vector to multiply
template<class T> inline vec3<T> vec3<T>::operator * (const vec3<T>& v) const {
  return vec3<T>(*this)*=v;
}


/// @brief Member to member division operator for two vectors
/// @tparam T Class inside the vectors
/// @param [in] v Vector to multiply
template<class T> inline vec3<T> vec3<T>::operator / (const vec3<T>& v) const {
  return vec3<T>(*this)/=v;
}


/// @brief Sum operator for a vector and a constant
///
/// Add a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant to add
template<class T> inline vec3<T> vec3<T>::operator + (const T& a) const {
  return vec3<T>(*this)+=a;
}


/// @brief Difference operator for a vector and a constant
///
/// Subtract a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant to subtract
template<class T> inline vec3<T> vec3<T>::operator - (const T& a) const {
  return vec3<T>(*this)-=a;
}


/// @brief Multiplication operator for a vector and a constant
/// @tparam T Class inside the vector
/// @param [in] a Constant to multiply
template<class T> inline vec3<T> vec3<T>::operator * (const T& a) const {
  return vec3<T>(*this)*=a;
}


/// @brief Quotient operator for a vector and a constant
/// @tparam T Class inside the vector
/// @param [in] a Constant to divide by
template<class T> inline vec3<T> vec3<T>::operator / (const T& a) const {
  return vec3<T>(*this)/=a;
}


/// @brief Sum operator for a vector and a constant (reverse)
///
/// Add the vector to a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant
/// @param [in] u Vector to add
template<class T> inline vec3<T> operator + (const T& a, const vec3<T>& u) {
  return u+a;
}


/// @brief Difference operator for a vector and a constant (reverse)
///
/// Subtract the vector from a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant
/// @param [in] u Vector to subtract
template<class T> inline vec3<T> operator - (const T& a, const vec3<T>& u) {
  return vec3<T>(a) - u;
}


/// @brief Product operator for a vector and a constant (reverse)
/// @tparam T Class inside the vector
/// @param [in] a Constant
/// @param [in] u Vector to multiply
template<class T> inline vec3<T> operator * (const T& a, const vec3<T>& u) {
  return u*a;
}


// =============================================================================
// ==== SOME SPECIALIZATIONS VEC3<T> ===========================================
// =============================================================================


/// @brief Sum operator for an integer constant and a vector of doubles
///
/// Get a vector whose three components are the constant
/// @param [in] c Integer constant
/// @param [in] v Vector of doubles to add
/// @return Sum as a vector of doubles
inline vec3<double> operator + (const int c, const vec3<double> &v) {
  return vec3<double>( c+v.x, c+v.y, c+v.z );
}


/// @brief Sum operator for a double constant and a vector of integers
///
/// Get a vector whose three components are the constant
/// @param [in] c Double constant
/// @param [in] v Vector of integers to add
/// @return Sum as a vector of doubles
inline vec3<double> operator + (const double c, const vec3<int> &v) {
  return vec3<double>( c+v.x, c+v.y, c+v.z );
}


/// @brief Sum operator for a vector of integers and a double constant
///
/// Get a vector whose three components are the constant
/// @param [in] u Vector of integers
/// @param [in] c Double constant to add
/// @return Sum as a vector of doubles
inline vec3<double> operator + (const vec3<int> &u, const double c) {
  return c+u;
}


/// @brief Sum operator for a vector of doubles and an integer constant
///
/// Get a vector whose three components are the constant
/// @param [in] u Vector of doubles
/// @param [in] c Integer constant to add
/// @return Sum as a vector of doubles
inline vec3<double> operator + (const vec3<double> &u, const int c) {
  return c+u;
}


/// @brief Sum operator for a vector of integers and a vector of doubles
/// @param [in] u Vector of integers
/// @param [in] v Vector of doubles to add
/// @return Sum as a vector of doubles
inline vec3<double> operator + (const vec3<int> &u, const vec3<double> &v) {
  return vec3<double>( u.x+v.x, u.y+v.y, u.z+v.z );
}


/// @brief Sum operator for a vector of doubles and a vector of integers
/// @param [in] u Vector of doubles
/// @param [in] v Vector of integers to add
/// @return Sum as a vector of doubles
inline vec3<double> operator + (const vec3<double> &u, const vec3<int> &v) {
  return v+u;
}


/// @brief Difference operator between an integer constant and a vector of doubles
///
/// Get a vector whose three components are the constant
/// @param [in] c Integer constant
/// @param [in] v Vector of doubles to subtract
/// @return Difference as a vector of doubles
inline vec3<double> operator - (const int c, const vec3<double> &v) {
  return vec3<double>( c-v.x, c-v.y, c-v.z );
}


/// @brief Difference operator between a double constant and a vector of integers
///
/// Get a vector whose three components are the constant
/// @param [in] c Double constant
/// @param [in] v Vector of integers to subtract
/// @return Difference as a vector of doubles
inline vec3<double> operator - (const double c, const vec3<int> &v) {
  return vec3<double>( c-v.x, c-v.y, c-v.z );
}


/// @brief Difference operator between a vector of integers and a double constant
///
/// Get a vector whose three components are the constant
/// @param [in] u Vector of integers
/// @param [in] c Double constant to subtract
/// @return Difference as a vector of doubles
inline vec3<double> operator - (const vec3<int> &u, const double c) {
  return vec3<double>( u.x-c, u.y-c, u.z-c );
}


/// @brief Difference operator between a vector of doubles and an integer constant
///
/// Get a vector whose three components are the constant
/// @param [in] u Vector of doubles
/// @param [in] c Integer constant to subtract
/// @return Difference as a vector of doubles
inline vec3<double> operator - (const vec3<double> &u, const int c) {
  return vec3<double>( u.x-c, u.y-c, u.z-c );
}


/// @brief Difference operator between a vector of integers and a vector of doubles
/// @param [in] u Vector of integers
/// @param [in] v Vector of doubles to subtract
/// @return Difference as a vector of doubles
inline vec3<double> operator - (const vec3<int> &u, const vec3<double> &v) {
  return vec3<double>( u.x-v.x, u.y-v.y, u.z-v.z );
}


/// @brief Difference operator between a vector of doubles and a vector of integers
/// @param [in] u Vector of doubles
/// @param [in] v Vector of integers to subtract
/// @return Difference as a vector of doubles
inline vec3<double> operator - (const vec3<double> &u, const vec3<int> &v) {
  return vec3<double>( u.x-v.x, u.y-v.y, u.z-v.z );
}


/// @brief Product operator for an integer constant and a vector of doubles
/// @param [in] c Integer constant
/// @param [in] v Vector of doubles to multiply
/// @return Product as a vector of doubles
inline vec3<double> operator * (const int c, const vec3<double> &v) {
  return vec3<double>( c*v.x, c*v.y, c*v.z );
}


/// @brief Product operator for a double constant and a vector of integers
/// @param [in] c Double constant
/// @param [in] v Vector of integers to multiply
/// @return Product as a vector of doubles
inline vec3<double> operator * (const double c, const vec3<int> &v) {
  return vec3<double>( c*v.x, c*v.y, c*v.z );
}


/// @brief Product operator for a vector of integers and a double constant
/// @param [in] u Vector of integers
/// @param [in] c Double constant to multiply
/// @return Product as a vector of doubles
inline vec3<double> operator * (const vec3<int> &u, const double c) {
  return c*u;
}


/// @brief Product operator for a vector of doubles and an integer constant
/// @param [in] u Vector of doubles
/// @param [in] c Integer constant to multiply
/// @return Product as a vector of doubles
inline vec3<double> operator * (const vec3<double> &u, const int c) {
  return c*u;
}


/// @brief Member to member multiplication operator for a vector of doubles and a vector of integers
/// @param [in] u Vector of doubles
/// @param [in] v Vector of integers to multiply
/// @return Product as a vector of doubles
inline vec3<double> operator * (const vec3<double> &u, const vec3<int> &v) {
  return vec3<double>( u.x*v.x, u.y*v.y, u.z*v.z );
}


/// @brief Member to member multiplication operator for a vector of integers and a vector of doubles
/// @param [in] u Vector of integers
/// @param [in] v Vector of doubles to multiply
/// @return Product as a vector of doubles
inline vec3<double> operator * (const vec3<int> &u, const vec3<double> &v) {
  return v*u;
}


// =============================================================================
// ==== DIVISION OPERATOR (ONLY A FEW) =========================================
// =============================================================================


/// @brief Quotient operator of a vector of integers by a double constant
/// @param [in] u Vector of integers
/// @param [in] c Double constant to divide by
/// @return Quotient as a vector of doubles
inline vec3<double> operator / (const vec3<int>& u, const double c) {
  return (1.*u)/c;
}


/// @brief Quotient operator of a vector of doubles by an integer constant
/// @param [in] u Vector of doubles
/// @param [in] c Integer constant to divide by
/// @return Quotient as a vector of doubles
inline vec3<double> operator / (const vec3<double>& u, const int c) {
  return u/(1.*c);
}


/// @brief Member to member division operator of a vector of doubles by a vector of integers
/// @param [in] u Vector of doubles
/// @param [in] v Vector of integers to divide by
/// @return Quotient as a vector of doubles
inline vec3<double> operator / (const vec3<int>& u, const vec3<double>& v) {
  return (1.*u)/v;
}


/// @brief Member to member division operator of a vector of integers by a vector of doubles
/// @param [in] u Vector of integers
/// @param [in] v Vector of doubles to divide by
/// @return Quotient as a vector of doubles
inline vec3<double> operator / (const vec3<double>& u, const vec3<int>& v) {
  return u/(1.*v);
}
