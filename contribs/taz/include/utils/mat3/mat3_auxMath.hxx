/// @file
/// @brief Extending functions from auxMath.hpp to mat3 class


/// @brief Absolute value
/// @tparam Type inside the matrix
/// @param [in] m Input matrix
/// @return Absolute value on each component
template <typename T>
inline mat3<T> auxAbs(const mat3<T>& m) {
  return mat3<T>(auxAbs(m.m11), auxAbs(m.m12), auxAbs(m.m13), auxAbs(m.m21), auxAbs(m.m22), auxAbs(m.m23), auxAbs(m.m31), auxAbs(m.m32), auxAbs(m.m33));
}


/// @brief Exponential
/// @tparam Type inside the matrix
/// @param [in] m Input matrix
/// @return Exponential on each component
inline mat3<double> auxExp(const mat3<double>& m) {
  return mat3<double>(auxExp(m.m11), auxExp(m.m12), auxExp(m.m13), auxExp(m.m21), auxExp(m.m22), auxExp(m.m23), auxExp(m.m31), auxExp(m.m32), auxExp(m.m33));
}


/// @brief Floor function
/// @tparam Type inside the matrix
/// @param [in] m Input matrix
/// @return Floor on each component
template <typename T=double>
inline mat3<T> auxFloor(const mat3<double>& m) {
  return mat3<T>(auxFloor<T>(m.m11), auxFloor<T>(m.m12), auxFloor<T>(m.m13), auxFloor<T>(m.m21), auxFloor<T>(m.m22), auxFloor<T>(m.m23), auxFloor<T>(m.m31), auxFloor<T>(m.m32), auxFloor<T>(m.m33));
}


/// @brief Max function
/// @tparam Type inside the matrices
/// @param [in] a First input matrix
/// @param [in] b Second input matrix
/// @return Component to component maximum
template <typename T>
inline mat3<T> auxMax(const mat3<T>& a, const mat3<T>& b) {
  return mat3<T>(auxMax(a.m11, b.m11), auxMax(a.m12, b.m12), auxMax(a.m13, b.m13), auxMax(a.m21, b.m21), auxMax(a.m22, b.m22), auxMax(a.m23, b.m23), auxMax(a.m31, b.m31), auxMax(a.m32, b.m32), auxMax(a.m33, b.m33));
}


/// @brief Min function
/// @tparam Type inside the matrices
/// @param [in] a First input matrix
/// @param [in] b Second input matrix
/// @return Component to component minimum
template <typename T>
inline mat3<T> auxMin(const mat3<T>& a, const mat3<T>& b) {
  return mat3<T>(auxMin(a.m11, b.m11), auxMin(a.m12, b.m12), auxMin(a.m13, b.m13), auxMin(a.m21, b.m21), auxMin(a.m22, b.m22), auxMin(a.m23, b.m23), auxMin(a.m31, b.m31), auxMin(a.m32, b.m32), auxMin(a.m33, b.m33));
}


/// @brief Modulo function
/// @tparam Type inside the matrices
/// @param [in] a First input matrix
/// @param [in] b Second input matrix
/// @return Component to component modulo
inline mat3<int> auxMod(const mat3<int>& a, const mat3<int>& b) {
  return mat3<int>(auxMod(a.m11, b.m11), auxMod(a.m12, b.m12), auxMod(a.m13, b.m13), auxMod(a.m21, b.m21), auxMod(a.m22, b.m22), auxMod(a.m23, b.m23), auxMod(a.m31, b.m31), auxMod(a.m32, b.m32), auxMod(a.m33, b.m33));
}


/// @brief Square root function
/// @tparam Type inside the matrix
/// @param [in] m Input matrix
/// @return Square on each component
template <typename T>
inline mat3<T> auxSq(const mat3<T>& m) {
  return mat3<T>(auxSq(m.m11), auxSq(m.m12), auxSq(m.m13), auxSq(m.m21), auxSq(m.m22), auxSq(m.m23), auxSq(m.m31), auxSq(m.m32), auxSq(m.m33));
}
