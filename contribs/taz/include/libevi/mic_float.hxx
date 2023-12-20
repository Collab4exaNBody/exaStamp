/// @file
/// @brief Vectorization using intrinsics



namespace simd {

  /// @brief Specialization of simd::vector for floats on mics
  template <> 
  class vector<MIC, float> {

  public:



    typedef vector<MIC, float> self_type;
    typedef float base_t ;
    typedef __m512 chunk_t;
    typedef __mmask16 mask_t;

    static constexpr uint16_t size  = mic::size;
    static constexpr uint16_t align = mic::align;

    static constexpr uint16_t chunk_size = get_chunk_size<base_t, size>();

    static inline const self_type& zero() { return s_zero; };
    static inline const self_type& one () { return s_one ; };





    vector() {}
    ~vector() {}

    vector(const base_t& in) : m_data(_mm512_set1_ps(in)) {}
    vector(const base_t* in) : m_data(_mm512_load_ps(in)) {}

    vector(const chunk_t& c) : m_data(c) {}

    vector(const self_type& a) : m_data(a.m_data) {}

    self_type& operator = (const self_type& a) {
      m_data = a.m_data;
      return *this;
    }

    void load (const base_t* in) { m_data = _mm512_load_ps(in); }

    void store(base_t* out) const { _mm512_store_ps(out, m_data); }

    self_type operator + (const self_type& a) const { return self_type( _mm512_add_ps(m_data, a.m_data) ); }
    self_type operator - (const self_type& a) const { return self_type( _mm512_sub_ps(m_data, a.m_data) ); }
    self_type operator * (const self_type& a) const { return self_type( _mm512_mul_ps(m_data, a.m_data) ); }
    self_type operator / (const self_type& a) const { return self_type( _mm512_div_ps(m_data, a.m_data) ); }



    mask_t operator == (const self_type& a) const { return _mm512_cmpeq_ps_mask (m_data, a.m_data); }
    mask_t operator != (const self_type& a) const { return _mm512_cmpneq_ps_mask(m_data, a.m_data); }
    mask_t operator >= (const self_type& a) const { return _mm512_cmpnlt_ps_mask(m_data, a.m_data); }
    mask_t operator >  (const self_type& a) const { return _mm512_cmpnle_ps_mask(m_data, a.m_data); }
    mask_t operator <= (const self_type& a) const { return _mm512_cmple_ps_mask (m_data, a.m_data); }
    mask_t operator <  (const self_type& a) const { return _mm512_cmplt_ps_mask (m_data, a.m_data); }

    chunk_t& data() { return m_data; }
    const chunk_t& data() const { return m_data; }

  private:

    static self_type s_zero;
    static self_type s_one ;

    chunk_t m_data;

  };



#define TMPL template <>
#define TMPL_vector vector<MIC, float>

//For MEAM
// No vector function.
#ifndef NDEBUG
  TMPL   inline void print(const TMPL_vector& a, const uint16_t n) {
    TMPL_vector::base_t *tmp=(TMPL_vector::base_t *)(&a.data());
    for (uint16_t i=0; i<n ; i++) std::cout<<tmp[i]<<std::endl;
  }
#endif


  TMPL inline TMPL_vector ceil    (const TMPL_vector& a) {
    return TMPL_vector( _mm512_ceil_ps(a.data()) );
  }

  TMPL inline TMPL_vector floor   (const TMPL_vector& a) {
    return TMPL_vector( _mm512_floor_ps(a.data()) );
  }

  TMPL inline TMPL_vector max     (const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm512_max_ps(a.data(), b.data()) );
  }

  TMPL inline TMPL_vector min     (const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm512_min_ps(a.data(), b.data()) );
  }

  TMPL inline TMPL_vector sqrt    (const TMPL_vector& a) {
    return TMPL_vector( _mm512_sqrt_ps(a.data()) );
  }

  TMPL inline TMPL_vector cbrt    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_cbrt_ps(a.data()) );
#else
    return no_vec::cbrt(a);
#endif
  }

  TMPL inline TMPL_vector inv_sqrt(const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_invsqrt_ps(a.data()) );
#else
    return TMPL_vector::one() / sqrt(a);
#endif
  }

  TMPL inline TMPL_vector inv_cbrt(const TMPL_vector& a) {



    return TMPL_vector::one() / cbrt(a);

  }



  TMPL inline TMPL_vector exp     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_exp_ps(a.data()) );
#else
    return no_vec::exp(a);
#endif
  }

  TMPL inline TMPL_vector exp2    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_exp2_ps(a.data()) );
#else
    return no_vec::exp2(a);
#endif
  }

  TMPL inline TMPL_vector exp10   (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_exp10_ps(a.data()) );
#else
    return no_vec::exp10(a);
#endif
  }

  TMPL inline TMPL_vector log     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_log_ps(a.data()) );
#else
    return no_vec::log(a);
#endif
  }

  TMPL inline TMPL_vector log2    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_log2_ps(a.data()) );
#else
    return no_vec::log2(a);
#endif
  }

  TMPL inline TMPL_vector log10   (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_log10_ps(a.data()) );
#else
    return no_vec::log10(a);
#endif
  }

  TMPL inline TMPL_vector pow     (const TMPL_vector& a, const TMPL_vector& b) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_pow_ps(a.data(), b.data()) );
#else
    return no_vec::pow(a, b);
#endif
  }



  TMPL inline TMPL_vector cos     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_cos_ps(a.data()) );
#else
    return no_vec::cos(a);
#endif
  }
  
  TMPL inline TMPL_vector sin     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_sin_ps(a.data()) );
#else
    return no_vec::sin(a);    
#endif
  }
  
  TMPL inline TMPL_vector tan     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_tan_ps(a.data()) );
#else
    return no_vec::tan(a);
#endif
    
  }

  TMPL inline TMPL_vector cosh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_cosh_ps(a.data()) );
#else
    return no_vec::cosh(a);
#endif
  }
  
  TMPL inline TMPL_vector sinh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_sinh_ps(a.data()) );
#else
    return no_vec::sinh(a);    
#endif
  }
  
  TMPL inline TMPL_vector tanh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_tanh_ps(a.data()) );
#else
    return no_vec::tanh(a);
#endif
  }



  TMPL inline TMPL_vector hypot   (const TMPL_vector& a, const TMPL_vector& b) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_hypot_ps(a.data(), b.data()) );
#else
    return sqrt(a*a+b*b);
#endif
  }



  TMPL inline TMPL_vector acos     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_acos_ps(a.data()) );
#else
    return no_vec::acos(a);
#endif
  }
  
  TMPL inline TMPL_vector asin     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_asin_ps(a.data()) );
#else
    return no_vec::asin(a);    
#endif
  }
  
  TMPL inline TMPL_vector atan     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_atan_ps(a.data()) );
#else
    return no_vec::atan(a);
#endif
    
  }

  TMPL inline TMPL_vector acosh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_acosh_ps(a.data()) );
#else
    return no_vec::acosh(a);
#endif
  }
  
  TMPL inline TMPL_vector asinh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_asinh_ps(a.data()) );
#else
    return no_vec::asinh(a);    
#endif
  }
  
  TMPL inline TMPL_vector atanh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm512_atanh_ps(a.data()) );
#else
    return no_vec::atanh(a);
#endif
  }



  TMPL inline TMPL_vector fmadd   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return a*b + c;
#else
    return TMPL_vector( _mm512_fmadd_ps(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fmsub   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return a*b - c;
#else
    return TMPL_vector( _mm512_fmsub_ps(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fnmadd  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return c - a*b;
#else
    return TMPL_vector( _mm512_fnmadd_ps(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fnmsub  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return TMPL_vector( _mm512_xor_ps((a*b+c).data(), _mm512_set1_ps(-0.f)) );
#else
    return TMPL_vector( _mm512_fnmadd_ps(a.data(), b.data(), c.data()) );
#endif
  }
  

  
#undef TMPL_vector
#undef TMPL



// inline TMPL_vector::mask_t TMPL_vector::operator == (const self_type& a) const {
//   return _mm512_cmp_ps_mask(m_data, a.m_data, _CMP_EQ_OQ);
// }

// inline TMPL_vector::mask_t TMPL_vector::operator != (const self_type& a) const {
//   return _mm512_cmp_ps_mask(m_data, a.m_data, _CMP_NEQ_OQ);
// }

// inline TMPL_vector::mask_t TMPL_vector::operator >= (const self_type& a) const {
//   return _mm512_cmp_ps_mask(m_data, a.m_data, _CMP_GE_OQ);
// }

// inline TMPL_vector::mask_t TMPL_vector::operator > (const self_type& a) const {
//   return _mm512_cmp_ps_mask(m_data, a.m_data, _CMP_GT_OQ);
// }

// inline TMPL_vector::mask_t TMPL_vector::operator <= (const self_type& a) const {
//   return _mm512_cmp_ps_mask(m_data, a.m_data, _CMP_LE_OQ);
// }

// inline TMPL_vector::mask_t TMPL_vector::operator < (const self_type& a) const {
//   return _mm512_cmp_ps_mask(m_data, a.m_data, _CMP_LT_OQ);
// }



#define TMPL_vector vector<MIC, float>
#define TMPL_mask vector<MIC, float>::mask_t



  inline TMPL_vector add (const TMPL_mask& mask, const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm512_maskz_add_ps(mask, a.data(), b.data()) );
  }
  
  inline TMPL_vector sub (const TMPL_mask& mask, const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm512_maskz_sub_ps(mask, a.data(), b.data()) );
  }
  
  inline TMPL_vector mul (const TMPL_mask& mask, const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm512_maskz_mul_ps(mask, a.data(), b.data()) );
  }
  
  inline TMPL_vector div (const TMPL_mask& mask, const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm512_maskz_div_ps(mask, a.data(), b.data()) );
  }
  
  inline TMPL_vector blend (const TMPL_mask& mask, const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm512_mask_blend_ps(mask, a.data(), b.data()) );
  }

  inline TMPL_vector blendv (const TMPL_vector& a, const TMPL_vector& b, const TMPL_mask& mask)  {
    return TMPL_vector( _mm512_mask_blend_ps( mask, a.data(), b.data()));
  }

 inline TMPL_vector::base_t reduce_add (const TMPL_vector& a, const TMPL_mask& mask) {
  #ifdef __evi_use_svml
	return _mm512_mask_reduce_add_ps(mask, a.data());
  #else
    return no_vec::reduce_add(a,mask);
  #endif
   }

   inline int movemask (const TMPL_mask mask) {
  #ifdef __evi_use_svml
    return _mm512_mask2int(mask);
  #else
    return no_vec::movemask(a);
  #endif
  }
  
#undef TMPL_vector
#undef TMPL_mask



}
