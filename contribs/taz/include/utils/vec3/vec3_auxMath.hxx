/// @file 
/// @brief Extending functions from auxMath.hpp to vec3 class


/// @brief Absolute value
/// @tparam Type inside the vector
/// @param [in] v Input vector
/// @return Absolute value on each component
template <typename T>
inline vec3<T> auxAbs(const vec3<T>& v) {
  return vec3<T>(auxAbs(v.x), auxAbs(v.y), auxAbs(v.z));
}


/// @brief Exponential
/// @tparam Type inside the vector
/// @param [in] v Input vector
/// @return Exponential on each component
inline vec3<double> auxExp(const vec3<double>& v) {
  return vec3<double>(auxExp(v.x), auxExp(v.y), auxExp(v.z));
}


/// @brief Floor function
/// @tparam Type inside the vector
/// @param [in] v Input vector
/// @return Floor on each component
template <typename T=double>
inline vec3<T> auxFloor(const vec3<double>& v) {
  return vec3<T>(auxFloor<T>(v.x), auxFloor<T>(v.y), auxFloor<T>(v.z));
}


/// @brief Max function
/// @tparam Type inside the vectors
/// @param [in] a First input vector
/// @param [in] b Second input vector
/// @return Component to component maximum
template <typename T>
inline vec3<T> auxMax(const vec3<T>& a, const vec3<T>& b) {
  return vec3<T>(auxMax(a.x, b.x), auxMax(a.y, b.y), auxMax(a.z, b.z));
}


/// @brief Min function
/// @tparam Type inside the vectors
/// @param [in] a First input vector
/// @param [in] b Second input vector
/// @return Component to component minimum
template <typename T>
inline vec3<T> auxMin(const vec3<T>& a, const vec3<T>& b) {
  return vec3<T>(auxMin(a.x, b.x), auxMin(a.y, b.y), auxMin(a.z, b.z));
}


/// @brief Modulo function
/// @param [in] a Input vector
/// @param [in] m Vector to modulate by
/// @return Component to component modulo
inline vec3<int> auxMod(const vec3<int>& a, const vec3<int>& m) {
  return vec3<int>( auxMod(a.x, m.x), auxMod(a.y, m.y), auxMod(a.z, m.z) );
}


/// @brief Square function
/// @tparam Type inside the vector
/// @param [in] v Input vector
/// @return Square on each component
template <typename T>
inline vec3<T> auxSq(const vec3<T>& v) {
	return vec3<T>(auxSq(v.x), auxSq(v.y), auxSq(v.z));
}
