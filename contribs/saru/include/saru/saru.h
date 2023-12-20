/// @file
/// @brief Implementation of functions that generate random numbers

#ifndef __SARU_HPP_INCLUDED
#define __SARU_HPP_INCLUDED


#include <cmath>
#include <limits.h>
#include <random>
#include <stdlib.h>


/// @brief Namespace for the random number generation
namespace Saru {

  /// @brief 2*RAND_MAX into a double
  static constexpr double RAND_MAX_PLS_TWO         = 2.0 + (double) RAND_MAX;
  /// @brief Inverse of 2*RAND_MAX+3
  static constexpr double INV_2_RAND_MAX_PLS_THREE = 1.0 / (2.0*((double) RAND_MAX) + 3.0);
  /// @brief Two pi
  static constexpr double TWO_PI                   = 2.0 * M_PI;

  /// @brief RAND_MAX into an int32_t
  static constexpr int32_t int32_max = 2147483647;

  int32_t saru (const int32_t seed1, const int32_t seed2, const int32_t seed3);

  double randU (const int32_t i, const int32_t j, const int32_t nit);
  double randN (const int32_t i, const int32_t j, const int32_t nit);

  double randU (const int32_t i, const int32_t nit);
  double randN (const int32_t i, const int32_t nit);

}


/// @brief Generat a random integer from the Saru algorithme
///
/// Generate a random integer between -RAND_MAX et RAND_MAX (cf
/// Afshar, Schmid, Pishevar, Worley in Comp. Phy. Comm. , 2012) @cite afshar_2012
/// @param [in] seed1 First seed
/// @param [in] seed2 Second seed
/// @param [in] seed3 Third seed
/// @warning Must have seed2 > seed1
inline int32_t Saru::saru (int32_t seed1, int32_t seed2, int32_t seed3) {

  const int32_t oWeylOffset = 0x8009d14b;
  const int32_t oWeylPeriod = 0xda879add;

  int32_t state, wstate;

  // Seeding Saru
  seed3 ^= (seed1<<7) ^ (seed2>>6);
  seed2 += (seed1>>4) ^ (seed3>>15);
  seed1 ^= (seed2<<9) + (seed3<<8 );
  seed3 ^= 0xA5366B4D * ( (seed2>>11) ^ (seed1>>15) );
  seed2 += 0x72BE1579 * ( (seed1<<4 ) ^ (seed3>>15) );
  seed1 ^= 0X3F38A6ED * ( (seed3>>5 ) ^ (seed2>>22) );
  seed2 += seed1 * seed3;
  seed1 += seed3 ^ (seed2>>2);
  seed2 ^= seed2 >> 17;

  // Transforming the seeds
  state  = 0x79dedea3 * ( seed1 ^ (seed1>>14) );
  wstate = (state + seed2) ^ (state>> 8);
  state  = state + ( wstate * (wstate^0xdddf97f5) );
  wstate = 0xABCB96F7 + (wstate >> 1);

  // Advancing Saru's state
  state  = 0x4beb5d59 * state + 0x2600e1f7;                       // LCG
  wstate = wstate + oWeylOffset + ( (wstate>>31) & oWeylPeriod ); // OWS

  uint32_t v = (state ^ (state>>26)) + wstate;
  return (v ^ (v>>20)) * 0x6957f5a7;

}


/// @brief Generate a random integer between 0 and 1 from an uniform distribution
/// @param [in] i First seed
/// @param [in] j Second seed
/// @param [in] nit Third seed
/// @warning Must have j > i
/// @version Tree seeds
inline double Saru::randU (const int32_t i, const int32_t j, const int32_t nit) {
  return ( (double) saru(i, j, nit) + RAND_MAX_PLS_TWO ) * INV_2_RAND_MAX_PLS_THREE;
}


/// @brief Generate a random integer between 0 and 1 from a normal distribution
/// @param [in] i First seed
/// @param [in] j Second seed
/// @param [in] nit Third seed
/// @warning Must have j > i
/// @version Tree seeds
inline double Saru::randN (const int32_t i, const int32_t j, const int32_t nit) {
  double tmp1 = randU((i<<1),    (j<<1),    nit);
  double tmp2 = randU((i<<1)+ 1, (j<<1)+ 1, nit);
  return sqrt( -2.0*log(tmp1) ) * cos( TWO_PI*tmp2 );
}


/// @brief Generate a random integer between 0 and 1 from an uniform distribution
/// @param [in] i First seed
/// @param [in] nit Third seed
/// @version Two seeds and RAND_MAX
inline double Saru::randU (const int32_t i, const int32_t nit) {
  return randU(i, int32_max, nit);
}


/// @brief Generate a random integer between 0 and 1 from a normal distribution
/// @param [in] i First seed
/// @param [in] nit Third seed
/// @version Two seeds and RAND_MAX
inline double Saru::randN (const int32_t i, const int32_t nit) {

  double tmp1 = randU((i<<1),   nit);
  double tmp2 = randU((i<<1)+1, nit);

  // if (tmp1<=0.) {
  //   std::cerr<< "\n" << "Fatal error in randN(int32_t i, int32_t nit) : a NaN is on the way ... \n"
  // 	     << "  randN(i, nit)     called with i=" << std::setw(11) << i      
  // 	     << ",              "                  << ", nit=" << std::setw(11) << nit << "\n"
  // 	     << "  randU(i, nit)     called with i=" << std::setw(11) << (i<<1) 
  // 	     << ",              "                  << ", nit=" << std::setw(11) << nit << " and returned " << tmp1 << "\n"
  // 	     << "  Saru (i, j, nit)  called with i=" << std::setw(11) << (i<<1) 
  // 	     << ", j=" << std::setw(11) << int32_max << ", nit=" << std::setw(11) << nit << " and returned " << saru((i<<1), int32_max, nit)
  // 	     << std::endl;
  //   std::abort();
  // }

  return sqrt( -2.0*log(tmp1) ) * cos( TWO_PI*tmp2 );

}

#endif // __SARU_HPP_INCLUDED
