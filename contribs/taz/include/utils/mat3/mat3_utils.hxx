/// @file
/// @brief Utility functions for mat3


/// @brief Euclidian dot product
/// @tparam T Type inside the matrices
/// @param [in] m First input matrix
/// @param [in] n Second input matrix
/// @return Dot product of the matrices
template <typename T> 
T dot(const mat3<T> &m, const mat3<T> &n) {
  return m.m11*n.m11 + m.m12*n.m12 + m.m13*n.m13 + m.m21*n.m21 + m.m22*n.m22 + m.m23*n.m23 + m.m31*n.m31 + m.m32*n.m32 + m.m33*n.m33 ;
}

/// @brief Euclidian dot product beween a matrix and a vector
/// @tparam T Type inside
/// @param [in] m Input matrix
/// @param [in] v Input vector
/// @return Dot product
template<typename T> 
vec3<T> dot(const mat3<T>& m, const vec3<T>& v) {
  return vec3<T>(m.m11*v.x + m.m12*v.y + m.m13*v.z,
		 m.m21*v.x + m.m22*v.y + m.m23*v.z,
		 m.m31*v.x + m.m32*v.y + m.m33*v.z);
}

/// @brief Matrix product
/// @tparam T Type inside the matrices
/// @param [in] m First input matrix
/// @param [in] n Second input matrix
/// @return Product of the matrices
template <typename T> 
mat3<T> multiply(const mat3<T> &m, const mat3<T> &n) {
  return mat3<T>( m.m11*n.m11 + m.m12*n.m21 + m.m13*n.m31,
		  m.m11*n.m12 + m.m12*n.m22 + m.m13*n.m32,
		  m.m11*n.m13 + m.m12*n.m23 + m.m13*n.m33,
		  m.m21*n.m11 + m.m22*n.m21 + m.m23*n.m31,
		  m.m21*n.m12 + m.m22*n.m22 + m.m23*n.m32,
		  m.m21*n.m13 + m.m22*n.m23 + m.m23*n.m33,
		  m.m31*n.m11 + m.m32*n.m21 + m.m33*n.m31,
		  m.m31*n.m12 + m.m32*n.m22 + m.m33*n.m32,
		  m.m31*n.m13 + m.m32*n.m23 + m.m33*n.m33);    
}


/// @brief Square of the euclidian norm
/// @tparam T Type inside
/// @param [in] m Input matrix
/// @return Squared norm
template <typename T> 
T norm2(const mat3<T> &m) {
  return dot(m,m);
}


/// @brief Matrix product 3x1 to 1x3 to produce 3x3 matrix
/// @tparam T Type inside
/// @param [in] u 3x1 matrix
/// @param [in] v 1x3 matrix
/// @return Result matrix
template <typename T> 
mat3<T> tensor(const vec3<T>& u, const vec3<T>& v) {
  return mat3<T>(u.x*v.x, u.x*v.y, u.x*v.z, u.y*v.x, u.y*v.y, u.y*v.z, u.z*v.x, u.z*v.y, u.z*v.z);
}


/// @brief Display on standard output (for debugging purposes)
/// @tparam Type inside the matrix
/// @param [in] m Input matrix
/// @param [in] separator String separator between the components
template <typename T> 
void print(const mat3<T> &m, std::string separator) {
  std::cout<< m.m11 << separator << m.m12 << separator << m.m13 << std::endl
	   << m.m21 << separator << m.m22 << separator << m.m23 << std::endl
	   << m.m31 << separator << m.m32 << separator << m.m33 << std::endl;
}


/// @brief Count number of occurrence of a value in the components of a matrix
/// @tparam T Type inside the matrix
/// @param [in] val Wanted value
/// @param [in] m Matrix to study
/// @return Number of occurrence of the value
template <typename T>
int countNumberOf(T val, const mat3<T> &m) {
  int res=0;
  if (m.m11==val) ++res;
  if (m.m12==val) ++res;
  if (m.m13==val) ++res;
  if (m.m21==val) ++res;
  if (m.m22==val) ++res;
  if (m.m23==val) ++res;
  if (m.m31==val) ++res;
  if (m.m32==val) ++res;
  if (m.m33==val) ++res;
  return res;
}


/// @brief Components sum
/// @tparam Class inside the matrix
/// @param [in] m Matrix input
/// @return Sum of all the components of the matrix
template<class T> T sum(const mat3<T>& m) {
  return m.m11+m.m12+m.m13+m.m21+m.m22+m.m23+m.m31+m.m32+m.m33;
}


/// @brief Components product
/// @tparam Class inside the matrix
/// @param [in] m Input matrix
/// @return Product of all the components of the matrix
template<class T> T product(const mat3<T>& m) {
  return m.m11*m.m12*m.m13*m.m21*m.m22*m.m23*m.m31*m.m32*m.m33;
}

/// @brief Matrix transposition
/// @tparam Class inside the matrix
/// @param [in] m Input matrix
/// @return Transposed matrix
template<class T> mat3<T> transpose(const mat3<T>& m){
  return mat3<T>(m.m11, m.m21, m.m31, m.m12, m.m22, m.m32, m.m13, m.m23, m.m33);
}
