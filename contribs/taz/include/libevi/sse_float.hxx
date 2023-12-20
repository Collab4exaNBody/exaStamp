/// @file
/// @brief Vectorization using intrinsics

#include <iostream>

namespace simd {

  /// @brief Specialization of simd::vector for SSE floats
  template <> 
  class vector<SSE, float> {

  public:



    typedef vector<SSE, float> self_type;
    typedef float  base_t;
    typedef __m128 chunk_t;


    static constexpr uint16_t size  = sse::size;
    static constexpr uint16_t align = sse::align;

    static constexpr uint16_t chunk_size = get_chunk_size<base_t, size>();

    static constexpr self_type& zero() { return s_zero; };
    static constexpr self_type& one () { return s_one ; };





    vector() {}
    ~vector() {}

    vector(const base_t& in) : m_data(_mm_set1_ps(in)) {}
    vector(const base_t* in) : m_data(_mm_load_ps(in)) {}

    vector(const chunk_t& c) : m_data(c) {}

    vector(const self_type& a) : m_data(a.m_data) {}

    self_type& operator = (const self_type& a) {
      m_data = a.m_data;
      return *this;
    }

    void load (const base_t* in) { m_data = _mm_load_ps(in); }

    void store(base_t* out) const { _mm_store_ps(out, m_data); }

    self_type operator + (const self_type& a) const { return self_type( _mm_add_ps(m_data, a.m_data) ); }
    self_type operator - (const self_type& a) const { return self_type( _mm_sub_ps(m_data, a.m_data) ); }
    self_type operator * (const self_type& a) const { return self_type( _mm_mul_ps(m_data, a.m_data) ); }
    self_type operator / (const self_type& a) const { return self_type( _mm_div_ps(m_data, a.m_data) ); }

    self_type operator && (const self_type& a) const { return self_type( _mm_and_ps(m_data, a.m_data) ); }
    self_type operator || (const self_type& a) const { return self_type( _mm_or_ps(m_data, a.m_data) ); }

    self_type operator == (const self_type& a) const { return self_type( _mm_cmpeq_ps (m_data, a.m_data) ); }
    self_type operator != (const self_type& a) const { return self_type( _mm_cmpneq_ps(m_data, a.m_data) ); }
    self_type operator >= (const self_type& a) const { return self_type( _mm_cmpge_ps (m_data, a.m_data) ); }
    self_type operator >  (const self_type& a) const { return self_type( _mm_cmpgt_ps (m_data, a.m_data) ); }
    self_type operator <= (const self_type& a) const { return self_type( _mm_cmple_ps (m_data, a.m_data) ); }
    self_type operator <  (const self_type& a) const { return self_type( _mm_cmplt_ps (m_data, a.m_data) ); }

    chunk_t& data() { return m_data; }
    const chunk_t& data() const { return m_data; }

  private:

    static self_type s_zero;
    static self_type s_one ;

    chunk_t m_data;

  };



#define TMPL template <>
#define TMPL_vector vector<SSE, float>

//For MEAM
// No vector function.
  TMPL inline TMPL_vector::base_t reduce_add (const TMPL_vector& a, const uint16_t n) {
    return no_vec::reduce_add(a,n);
  }
    
#ifndef NDEBUG
  TMPL   inline void print(const TMPL_vector& a, const uint16_t n) {
    TMPL_vector::base_t *tmp=(TMPL_vector::base_t *)(&a.data());
    for (uint16_t i=0; i<n ; i++) std::cout<<tmp[i]<<std::endl;
  }
#endif

  TMPL   inline int movemask (const TMPL_vector& a) {
  #if defined(__evi_use_svml) || defined(__evi_use_svml_gcc)
    return _mm_movemask_ps(a.data());
  #else
    return no_vec::movemask(a);
  #endif
  }

  TMPL inline TMPL_vector ceil    (const TMPL_vector& a) {
    return TMPL_vector( _mm_ceil_ps(a.data()) );
  }

  TMPL inline TMPL_vector floor   (const TMPL_vector& a) {
    return TMPL_vector( _mm_floor_ps(a.data()) );
  }

  TMPL inline TMPL_vector max     (const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm_max_ps(a.data(), b.data()) );
  }

  TMPL inline TMPL_vector min     (const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm_min_ps(a.data(), b.data()) );
  }



  TMPL inline TMPL_vector sqrt    (const TMPL_vector& a) {
    return TMPL_vector( _mm_sqrt_ps(a.data()) );
  }

  TMPL inline TMPL_vector cbrt    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_cbrt_ps(a.data()) );
#else
    return no_vec::cbrt(a);
#endif
  }

  TMPL inline TMPL_vector inv_sqrt(const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_invsqrt_ps(a.data()) );
#else
    return TMPL_vector::one() / sqrt(a);
#endif
  }

  TMPL inline TMPL_vector inv_cbrt(const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_invcbrt_ps(a.data()) );
#else
    return TMPL_vector::one() / cbrt(a);
#endif
  }



  TMPL inline TMPL_vector exp     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_exp_ps(a.data()) );
#else
    return no_vec::exp(a);
#endif
  }

  TMPL inline TMPL_vector exp2    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_exp2_ps(a.data()) );
#else
    return no_vec::exp2(a);
#endif
  }

  TMPL inline TMPL_vector exp10   (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_exp10_ps(a.data()) );
#else
    return no_vec::exp10(a);
#endif
  }

  TMPL inline TMPL_vector log     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_log_ps(a.data()) );
#else
    return no_vec::log(a);
#endif
  }

  TMPL inline TMPL_vector log2    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_log2_ps(a.data()) );
#else
    return no_vec::log2(a);
#endif
  }

  TMPL inline TMPL_vector log10   (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_log10_ps(a.data()) );
#else
    return no_vec::log10(a);
#endif
  }

  TMPL inline TMPL_vector pow     (const TMPL_vector& a, const TMPL_vector& b) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_pow_ps(a.data(), b.data()) );
#else
    return no_vec::pow(a, b);
#endif
  }



  TMPL inline TMPL_vector cos     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_cos_ps(a.data()) );
#else
    return no_vec::cos(a);
#endif
  }
  
  TMPL inline TMPL_vector sin     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_sin_ps(a.data()) );
#else
    return no_vec::sin(a);    
#endif
  }
  
  TMPL inline TMPL_vector tan     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_tan_ps(a.data()) );
#else
    return no_vec::tan(a);
#endif
    
  }

  TMPL inline TMPL_vector cosh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_cosh_ps(a.data()) );
#else
    return no_vec::cosh(a);
#endif
  }
  
  TMPL inline TMPL_vector sinh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_sinh_ps(a.data()) );
#else
    return no_vec::sinh(a);    
#endif
  }
  
  TMPL inline TMPL_vector tanh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_tanh_ps(a.data()) );
#else
    return no_vec::tanh(a);
#endif
  }



  TMPL inline TMPL_vector hypot   (const TMPL_vector& a, const TMPL_vector& b) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_hypot_ps(a.data(), b.data()) );
#else
    return sqrt(a*a+b*b);
#endif
  }



  TMPL inline TMPL_vector acos     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_acos_ps(a.data()) );
#else
    return no_vec::acos(a);
#endif
  }
  
  TMPL inline TMPL_vector asin     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_asin_ps(a.data()) );
#else
    return no_vec::asin(a);    
#endif
  }
  
  TMPL inline TMPL_vector atan     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_atan_ps(a.data()) );
#else
    return no_vec::atan(a);
#endif
    
  }

  TMPL inline TMPL_vector acosh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_acosh_ps(a.data()) );
#else
    return no_vec::acosh(a);
#endif
  }
  
  TMPL inline TMPL_vector asinh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_asinh_ps(a.data()) );
#else
    return no_vec::asinh(a);    
#endif
  }
  
  TMPL inline TMPL_vector atanh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_atanh_ps(a.data()) );
#else
    return no_vec::atanh(a);
#endif
  }



  TMPL inline TMPL_vector blendv  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& mask) {
    return TMPL_vector( _mm_blendv_ps(a.data(), b.data(), mask.data()) );
  }



  TMPL inline TMPL_vector fmadd   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return a*b + c;
#else
    return TMPL_vector( _mm_fmadd_ps(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fmsub   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return a*b - c;
#else
    return TMPL_vector( _mm_fmsub_ps(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fnmadd  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return c - a*b;
#else
    return TMPL_vector( _mm_fnmadd_ps(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fnmsub  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return TMPL_vector( _mm_xor_ps((a*b+c).data(), _mm_set1_ps(-0.f)) );
#else
    return TMPL_vector( _mm_fnmadd_ps(a.data(), b.data(), c.data()) );
#endif
  }
  

  
#undef TMPL_vector
#undef TMPL



}
