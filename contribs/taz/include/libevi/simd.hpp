/// @file
/// @brief Main Header

#ifndef __SIMD_HPP_INCLUDED
#define __SIMD_HPP_INCLUDED

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <assert.h>

#ifndef __evi_ser
#include <emmintrin.h>
#include <immintrin.h>
#endif

#ifndef NDEBUG
#include <iostream>
#endif

#if defined __INTEL_COMPILER



#define __INTEL_CXX_VERSION (__INTEL_COMPILER*100 + __INTEL_COMPILER_UPDATE)

#if __INTEL_CXX_VERSION < 150001
#error This version of the code requires icc 15.0.1 or higher
#endif

#define __USE_SVML



#elif defined __GNUG__


#ifdef __clang__
#define __GNU_CXX_VERSION (8*10000 + 3*100 )
#else
#define __GNU_CXX_VERSION (__GNUC__*10000 + __GNUC_MINOR__*100 + __GNUC_PATCHLEVEL__)
#endif

#if __GNU_CXX_VERSION < 40900
#error This version of the code requires g++ 4.9.0 or higher
#else
#define __evi_use_svml_gcc
#endif



#endif





/// @brief Flag to check the use of SVML
#ifdef __INTEL_COMPILER
#define __evi_use_svml
#endif





namespace simd {

  /// @brief All vectorization types 
  enum instr_t { SER, SSE, AVX, MIC };



  namespace ser {
    constexpr uint16_t size  = 0;
    constexpr uint16_t align = 8;
  }
  namespace sse {
    constexpr uint16_t size  = 128;
    constexpr uint16_t align = 16;
  }
  namespace avx {
    constexpr uint16_t size  = 256;
    constexpr uint16_t align = 32;
  }
  namespace mic {
    constexpr uint16_t size  = 512;
    constexpr uint16_t align = 64;
  }



  template <instr_t instr, class base_t> 
  class vector;

  template <> class vector<SER, float >;
  template <> class vector<SER, double>;

#define TMPL template <instr_t instr, class base_t> 
#define TMPL_vector vector<instr, base_t>

  TMPL TMPL_vector ceil    (const TMPL_vector& a);
  TMPL TMPL_vector floor   (const TMPL_vector& a);
  TMPL TMPL_vector max     (const TMPL_vector& a, const TMPL_vector& b);
  TMPL TMPL_vector min     (const TMPL_vector& a, const TMPL_vector& b);

  TMPL base_t reduce_add   (const TMPL_vector& a, const uint16_t n);
  TMPL void print          (const TMPL_vector& a, const uint16_t n);
  TMPL int  movemask       (const TMPL_vector& a);

  TMPL TMPL_vector sqrt    (const TMPL_vector& a);
  TMPL TMPL_vector cbrt    (const TMPL_vector& a);
  TMPL TMPL_vector inv_sqrt(const TMPL_vector& a);
  TMPL TMPL_vector inv_cbrt(const TMPL_vector& a);

  TMPL TMPL_vector exp     (const TMPL_vector& a);
  TMPL TMPL_vector exp2    (const TMPL_vector& a);
  TMPL TMPL_vector exp10   (const TMPL_vector& a);
  TMPL TMPL_vector log     (const TMPL_vector& a);
  TMPL TMPL_vector log2    (const TMPL_vector& a);
  TMPL TMPL_vector log10   (const TMPL_vector& a);
  TMPL TMPL_vector pow     (const TMPL_vector& a, const TMPL_vector& b);

  TMPL TMPL_vector cos     (const TMPL_vector& a);
  TMPL TMPL_vector sin     (const TMPL_vector& a);
  TMPL TMPL_vector tan     (const TMPL_vector& a);
  TMPL TMPL_vector cosh    (const TMPL_vector& a);
  TMPL TMPL_vector sinh    (const TMPL_vector& a);
  TMPL TMPL_vector tanh    (const TMPL_vector& a);

  TMPL TMPL_vector hypot   (const TMPL_vector& a, const TMPL_vector& b);

  TMPL TMPL_vector acos    (const TMPL_vector& a);
  TMPL TMPL_vector asin    (const TMPL_vector& a);
  TMPL TMPL_vector atan    (const TMPL_vector& a);
  TMPL TMPL_vector acosh   (const TMPL_vector& a);
  TMPL TMPL_vector asinh   (const TMPL_vector& a);
  TMPL TMPL_vector atanh   (const TMPL_vector& a);

  TMPL TMPL_vector blendv  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& mask);

  TMPL TMPL_vector fmadd   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c);
  TMPL TMPL_vector fmsub   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c);
  TMPL TMPL_vector fnmadd  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c);
  TMPL TMPL_vector fnmsub  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c);


#ifndef __evi_use_svml
  namespace no_vec {

  TMPL base_t min_reduce   (const TMPL_vector& a, const uint16_t n);
  TMPL base_t mul_reduce   (const TMPL_vector& a, const uint16_t n);
  TMPL base_t add_reduce   (const TMPL_vector& a, const uint16_t n);


  TMPL TMPL_vector cbrt    (const TMPL_vector& a);


  TMPL TMPL_vector exp     (const TMPL_vector& a);
  TMPL TMPL_vector exp2    (const TMPL_vector& a);
  TMPL TMPL_vector exp10   (const TMPL_vector& a);
  TMPL TMPL_vector log     (const TMPL_vector& a);
  TMPL TMPL_vector log2    (const TMPL_vector& a);
  TMPL TMPL_vector log10   (const TMPL_vector& a);
  TMPL TMPL_vector pow     (const TMPL_vector& a, const TMPL_vector& b);

  TMPL TMPL_vector cos     (const TMPL_vector& a);
  TMPL TMPL_vector sin     (const TMPL_vector& a);
  TMPL TMPL_vector tan     (const TMPL_vector& a);
  TMPL TMPL_vector cosh    (const TMPL_vector& a);
  TMPL TMPL_vector sinh    (const TMPL_vector& a);
  TMPL TMPL_vector tanh    (const TMPL_vector& a);

  TMPL TMPL_vector acos    (const TMPL_vector& a);
  TMPL TMPL_vector asin    (const TMPL_vector& a);
  TMPL TMPL_vector atan    (const TMPL_vector& a);
  TMPL TMPL_vector acosh   (const TMPL_vector& a);
  TMPL TMPL_vector asinh   (const TMPL_vector& a);
  TMPL TMPL_vector atanh   (const TMPL_vector& a);



  }
  #endif
#undef TMPL_vector
#undef TMPL

  
  
  template <class T, uint16_t align> 
  static inline T* aligned_malloc(size_t n) 
  {
    void* memptr = nullptr;
    size_t alignment = align;
    if( alignment < sizeof(void*) ) { alignment = sizeof(void*); } // required by posix_memalign
    assert( posix_memalign( &memptr, alignment, n*sizeof(T) ) == 0 );
    assert( memptr != nullptr );
    return (T*) memptr;
  }

  template <class T> 
  static inline void aligned_free(T* ptr) {
    if (ptr!=nullptr) free(ptr);
  }

  template <class T, uint16_t size> 
  uint16_t constexpr get_chunk_size() {
    return size==0 ? 1 : size/(8*sizeof(T));
  }



}




#include "libevi/no_vec.hxx"
#include "libevi/ser_float.hxx"
#include "libevi/ser_double.hxx"






#ifdef __evi_sse

namespace simd {
  template <> class vector<SSE, float >;
  template <> class vector<SSE, double>;
}

#include "libevi/sse_float.hxx"
#include "libevi/sse_double.hxx"

#endif

#ifdef __evi_avx

namespace simd {
  template <> class vector<AVX, float >;
  template <> class vector<AVX, double>;
}

#include "libevi/avx_float.hxx"
#include "libevi/avx_double.hxx"

#endif

#ifdef __evi_mic

namespace simd {
  template <> class vector<MIC, float >;
  template <> class vector<MIC, double>;
}

#include "libevi/mic_float.hxx"
#include "libevi/mic_double.hxx"

#endif





namespace simd {

#if defined __evi_sse
  template <class T> using vector_t = vector<SSE, T>;
#elif defined __evi_avx
  template <class T> using vector_t = vector<AVX, T>;
#elif defined __evi_mic
  template <class T> using vector_t = vector<MIC, T>;
#else
  template <class T> using vector_t = vector<SER, T>;
#endif


//  static inline vector_t s_zero(0.0);
//  static inline vector_t s_one(1.0);
}





#endif // __SIMD_HPP_INCLUDED
