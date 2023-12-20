/// @file 
/// @brief Utility functions for vec3

/// @brief Euclidian dot product
/// @tparam T Type inside the vectors
/// @param [in] u First input vector
/// @param [in] v Second input vector
/// @return Dot product of the vectors
template <typename T> 
inline T dot(const vec3<T> &u, const vec3<T> &v) {
  return u.x*v.x + u.y*v.y + u.z*v.z;
}


/// @brief Cross product
/// @tparam T Type inside the vectors
/// @param [in] u First input vector
/// @param [in] v Second input vector
/// @return Cross product of the vectors
template <typename T> 
inline vec3<T> cross(const vec3<T> &u, const vec3<T> &v) {
  return vec3<T>( u.y*v.z - u.z*v.y,
		  u.z*v.x - u.x*v.z,
		  u.x*v.y - u.y*v.x);    
}


/// @brief Squared euclidian norm
/// @tparam T Type inside the vector
/// @param [in] v Input vector
/// @return Squared norm
template <typename T> 
inline T norm2(const vec3<T> &v) {
  return v.x*v.x + v.y*v.y + v.z*v.z;
}


/// @brief Count number of occurrence of a value in the components of a vector
/// @tparam T Type inside the vector
/// @param [in] val Wanted value
/// @param [in] v Vector to study
/// @return Number of occurrence of the value
template <typename T>
inline int countNumberOf(T val, const vec3<T> &v) {
  int res=0;
  if (v.x==val) ++res;
  if (v.y==val) ++res;
  if (v.z==val) ++res;
  return res;
}


/// @brief Components sum
/// @tparam T Class inside the vector
/// @param [in] v Input vector
/// @return Sum of all the components of the vector
template<class T> inline T sum(const vec3<T>& v) {
  return v.x + v.y + v.z;
}


/// @brief Components product
/// @tparam T Class inside the vector
/// @param [in] v Input vector
/// @return Product of all the components of the vector
template<class T> inline T product(const vec3<T>& v) {
  return v.x * v.y * v.z;
}
