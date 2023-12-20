/// @file 
/// @brief A set of basic math functions
///
/// Right now it is just an encapsulation of the \<math.h\> or \<cmath\> libraries,
/// but we do it in case we have to change those later (we could need a fast sqrt
/// or something similar).

#ifndef __AUX_MATH_HPP_INCLUDED
#define __AUX_MATH_HPP_INCLUDED


#include <cmath>
#include <iostream>


/// @brief Absolute value function
/// @param [in] a Argument integer
/// @return \f$ |a| \f$
inline int auxAbs (int a) {
  return a>0 ? a : -a;
}


/// @brief Absolute value function
/// @param [in] a Argument float
/// @return \f$ |a| \f$
inline float auxAbs (float a) {
  return fabs(a);
}


/// @brief Absolute value function
/// @param [in] a Argument double
/// @return \f$ |a| \f$
inline double auxAbs (double a) {
  return fabs(a);
}


/// @brief Cosine function
/// @param [in] a Argument double
/// @return cos(a)
inline double auxCos (double a) {
  return cos(a);
}


/// @brief Cubic root
/// @param [in] a Argument double
/// @return \f$ a^{1/3} \f$
inline double auxCbrt(double a) {
  return cbrt(a);
}


/// @brief Exponential function
/// @param [in] a Argument double
/// @return \f$ \exp \left( a \right) \f$
inline double auxExp (double a) {
  return exp(a);
}


/// @brief Floor function
/// @tparam T Return type
/// @param [in] a Argument double
/// @return \f$ \max\{ m \in \mathtt{Z} \,|\, m\leq a\}\f$
template <typename T>
inline T auxFloor (double a) {
  return static_cast<T>(floor(a));
}


/// @brief Floor function
///
/// Specialization for a double return
/// @param [in] a Argument double
/// @return \f$ \max\{ m \in \mathtt{Z} \,|\, m\leq a\}\f$
template <>
inline double auxFloor (double a) {
  return floor(a);
}


/// @brief logarithme fuction
/// @param [in] a Argument double
/// @return log(a)
inline double auxLog (double a) {
  return log(a);
}


/// @brief Max function
/// @tparam Type of the arguments and return
/// @param [in] a Argument T
/// @param [in] b Argument T
/// @return \f$ \max(a, b) \f$
template <typename T>
inline const T& auxMax (const T& a, const T& b) {
  return (a>b ? a:b);
}


/// @brief Min function
/// @tparam Type of the arguments and return
/// @param [in] a Argument T
/// @param [in] b Argument T
/// @return \f$ \min(a, b) \f$
template <typename T>
inline const T& auxMin (const T& a, const T& b) {
  return (a<b ? a:b);
}


/// @brief Max function
/// @tparam Type of the arguments and return
/// @param [in] a Argument T
/// @param [in] b Argument T
/// @param [in] c Argument T
/// @return \f$ \max(a, b, c) \f$
template <typename T>
inline const T& auxMax (const T& a, const T& b, const T&c) {
  return auxMax(auxMax(a, b), c);
}


/// @brief Min function
/// @tparam Type of the arguments and return
/// @param [in] a Argument T
/// @param [in] b Argument T
/// @param [in] c Argument T
/// @return \f$ \min(a, b, c) \f$
template <typename T>
inline const T& auxMin (const T& a, const T& b, const T&c) {
  return auxMin(auxMin(a, b), c);
}


/// @brief Modulo function
///
/// An integer modulo that handles negatives
/// @param [in] a Argument integer
/// @param [in] m Modulo integer
/// @return a modulo m
inline int auxMod (int a, int m) {
  return ( (a<0) ? ( (a%m)+m ) : a ) % m;
}


/// @brief Power function
/// @tparam T Type of the argument and return
/// @tparam U Type of the power
/// @param [in] t Argument
/// @param [in] u Power
/// @return \f$ t^u \f$
template <typename T, typename U>
inline T auxPow(const T& t, const U& u) {
  return pow(t, u);
}


/// @brief Sine function
/// @param [in] a Argument double
/// @return sin(a)
inline double auxSin (double a) {
  return sin(a);
}


/// @brief Square function
/// @tparam  T Type of the argument
/// @param [in] a Argument T
/// @return \f$ a^2 \f$
template <typename T>
inline const T auxSq (const T &a) {
  return a*a;
}


/// @brief Square root
/// @param [in] a Argument double
/// @return \f$ a^{1/2} \f$
inline double auxSqrt(double a) {
  return sqrt(a);
}


/// @brief Tangent function
/// @param [in] a Argument double
/// @return tan(a)
inline double auxTan (double a) {
  return tan(a);
}


/// @brief Round an integer to a multiple of a chunk
/// @tparam N Chunk
/// @param [in] a Argument integer
/// @return The smallest integer larger than the original integer and multiple of the chunk
template <uint N> uint roundToChunk(uint a) {
  return (1+a/N)*N;
}


/// @brief Compute the mean, the standard deviation, the min and the max of a
/// set of data (used only in debug, poor performances)
/// @tparam T Type of the data set
/// @param [in] val Pointer to the data set
/// @param [in] size Size of the data set
/// @param [out] mean Mean
/// @param [out] std Standard deviation
/// @param [out] min Minimum
/// @param [out] max Maximum
template <class T> void quickAnalysis(const T* val, const int size, double& mean, double& std, T& min ,T& max) {

  mean = 0.;
  std  = 0.;
  min  = val[0];
  max  = val[0];

  for (int i=0; i<size; ++i) {
    mean += val[i];
    min   = auxMin(min, val[i]);
    max   = auxMax(max, val[i]);
  }

  mean /= size;

  for (int i=0; i<size; ++i) std += auxSq(val[i]-mean);

  std /= size;
  std  = auxSqrt(std);

}

#endif // __AUX_MATH_HPP_INCLUDED

