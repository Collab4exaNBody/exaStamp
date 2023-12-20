/// @file 
/// @brief Comparison operators for vec3

/// @brief Equality operator between vectors
///
/// The three components are equal
/// @tparam T Class inside the vectors
/// @param [in] v Vector to compare
template<class T> inline bool vec3<T>::operator == (const vec3<T>& v) const {
  return x==v.x && y==v.y && z==v.z;
}


/// @brief Inequality operator between vectors
/// @tparam T Class inside the vectors
/// @param [in] v Vector to compare
template<class T> inline bool vec3<T>::operator != (const vec3<T>& v) const {
  return !(*this==v);
}


/// @brief Comparison (>) operator between vectors
///
/// The three components are strictly superior
/// @tparam T Class inside the vectors
/// @param [in] v Vector to compare
template<class T> inline vec3<bool> vec3<T>::operator >  (const vec3<T>& v) const { 
  return vec3<bool>(x>v.x, y>v.y, z>v.z); 
}


/// @brief Comparison (<) operator between vectors
///
/// The three components are strictly inferior
/// @tparam T Class inside the vectors
/// @param [in] v Vector to compare
template<class T> inline vec3<bool> vec3<T>::operator <  (const vec3<T>& v) const { 
  return vec3<bool>(x<v.x, y<v.y, z<v.z); 
}


/// @brief Comparison (>=) operator between vectors
///
/// The three components are superior
/// @tparam T Class inside the vectors
/// @param [in] v Vector to compare
template<class T> inline vec3<bool> vec3<T>::operator >= (const vec3<T>& v) const { 
  return vec3<bool>(x>=v.x, y>=v.y, z>=v.z); 
}


/// @brief Comparison (<=) operator between vectors
///
/// The three components are inferior
/// @tparam T Class inside the vectors
/// @param [in] v Vector to compare
template<class T> inline vec3<bool> vec3<T>::operator <= (const vec3<T>& v) const { 
  return vec3<bool>(x<=v.x, y<=v.y, z<=v.z); 
}
  

/// @brief Comparison (>) operator to a constant
///
/// Compare to a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant to compare
template<class T> inline vec3<bool> vec3<T>::operator >  (const T& a) const { 
  return (*this)>vec3<T>(a); 
}


/// @brief Comparison (<) operator to a constant
///
/// Compare to a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant to compare
template<class T> inline vec3<bool> vec3<T>::operator <  (const T& a) const { 
  return (*this)<vec3<T>(a); 
}


/// @brief Comparison (>=) operator to a constant
///
/// Compare to a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant to compare
template<class T> inline vec3<bool> vec3<T>::operator >= (const T& a) const { 
  return (*this)>=vec3<T>(a); 
}


/// @brief Comparison (<=) operator to a constant
///
/// Compare to a vector whose three components are the constant
/// @tparam T Class inside the vectors
/// @param [in] a Constant to compare
template<class T> inline vec3<bool> vec3<T>::operator <= (const T& a) const { 
  return (*this)<=vec3<T>(a); 
}
