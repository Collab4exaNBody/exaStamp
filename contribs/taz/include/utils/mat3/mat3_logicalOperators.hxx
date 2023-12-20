/// @file
/// @brief Logical operators for mat3


/// @brief Equality operator between matrices
///
/// The nine components are equal
/// @tparam T Class inside the matrices
/// @param [in] m Matrix to compare
template<class T> bool mat3<T>::operator == (const mat3<T>& m) const {
  return m11==m.m11 && m12==m.m12 && m13==m.m13 && m21==m.m21 && m22==m.m22 && m23==m.m23 && m31==m.m31 && m32==m.m32 && m33==m.m33 ;
}


/// @brief Inequality operator between matrices
/// @tparam T Class inside the matrices
/// @param [in] m Matrix to compare
template<class T> bool mat3<T>::operator != (const mat3<T>& m) const {
  return !(*this==m);
}


/// @brief Operator > (with another mat3)
///
/// The nine components are strictly superior
/// @tparam T Class inside the matrices
/// @param [in] m Matrix to compare
template<class T> mat3<bool> mat3<T>::operator >  (const mat3<T>& m) const { 
  return mat3<bool>(m11>m.m11, m12>m.m12, m13>m.m13, m21>m.m21, m22>m.m22, m23>m.m23, m31>m.m31, m32>m.m32, m33>m.m33); 
}


/// @brief Operator < (with another mat3)
///
/// The nine components are strictly inferior
/// @tparam T Class inside the matrices
/// @param [in] m Matrix to compare
template<class T> mat3<bool> mat3<T>::operator <  (const mat3<T>& m) const { 
  return mat3<bool>(m11<m.m11, m12<m.m12, m13<m.m13, m21<m.m21, m22<m.m22, m23<m.m23, m31<m.m31, m32<m.m32, m33<m.m33); 
}


/// @brief Operator >= (with another mat3)
///
/// The nine components are superior
/// @tparam T Class inside the matrices
/// @param [in] m Matrix to compare
template<class T> mat3<bool> mat3<T>::operator >= (const mat3<T>& m) const { 
  return mat3<bool>(m11>=m.m11, m12>=m.m12, m13>=m.m13, m21>=m.m21, m22>=m.m22, m23>=m.m23, m31>=m.m31, m32>=m.m32, m33>=m.m33); 
}


/// @brief Operator <= (with another mat3)
///
/// The nine components are inferior
/// @tparam T Class inside the matrices
/// @param [in] m Matrix to compare
template<class T> mat3<bool> mat3<T>::operator <= (const mat3<T>& m) const { 
  return mat3<bool>(m11<=m.m11, m12<=m.m12, m13<=m.m13, m21<=m.m21, m22<=m.m22, m23<=m.m23, m31<=m.m31, m32<=m.m32, m33<=m.m33); 
}
  

/// @brief Operator > (with a constant)
///
/// Compare to a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to compare
template<class T> mat3<bool> mat3<T>::operator >  (const T& a) const { 
  return (*this)>mat3<T>(a); 
}


/// @brief Operator < (with a constant)
///
/// Compare to a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to compare
template<class T> mat3<bool> mat3<T>::operator <  (const T& a) const { 
  return (*this)<mat3<T>(a); 
}


/// @brief Operator >= (with a constant)
///
/// Compare to a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to compare
template<class T> mat3<bool> mat3<T>::operator >= (const T& a) const { 
  return (*this)>=mat3<T>(a); 
}


/// @brief Operator <= (with a constant)
///
/// Compare to a matrix whose nine components are the constant
/// @tparam T Class inside the matrix
/// @param [in] a Constant to compare
template<class T> mat3<bool> mat3<T>::operator <= (const T& a) const { 
  return (*this)<=mat3<T>(a); 
}
