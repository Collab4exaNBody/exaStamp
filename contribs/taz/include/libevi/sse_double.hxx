/// @file
/// @brief Vectorization using intrinsics


namespace simd {

  /// @brief Specialization of simd::vector for SSE doubles
  template <> 
  class vector<SSE, double> {

  public:



    typedef vector<SSE, double> self_type;
    typedef double  base_t;
    typedef __m128d chunk_t;


    static constexpr uint16_t size  = sse::size;
    static constexpr uint16_t align = sse::align;

    static constexpr uint16_t chunk_size = get_chunk_size<base_t, size>();

    static constexpr self_type& zero() { return s_zero; };
    static constexpr self_type& one () { return s_one ; };





    vector() {}
    ~vector() {}

    vector(const base_t& in) : m_data(_mm_set1_pd(in)) {}
    vector(const base_t* in) : m_data(_mm_load_pd(in)) {}

    vector(const chunk_t& c) : m_data(c) {}

    vector(const self_type& a) : m_data(a.m_data) {}

    self_type& operator = (const self_type& a) {
      m_data = a.m_data;
      return *this;
    }

    void load (const base_t* in) { m_data = _mm_load_pd(in); }

    void store(base_t* out) const { _mm_store_pd(out, m_data); }

    self_type operator + (const self_type& a) const { return self_type( _mm_add_pd(m_data, a.m_data) ); }
    self_type operator - (const self_type& a) const { return self_type( _mm_sub_pd(m_data, a.m_data) ); }
    self_type operator * (const self_type& a) const { return self_type( _mm_mul_pd(m_data, a.m_data) ); }
    self_type operator / (const self_type& a) const { return self_type( _mm_div_pd(m_data, a.m_data) ); }

    self_type operator && (const self_type& a) const { return self_type( _mm_and_pd(m_data, a.m_data) ); }
    self_type operator || (const self_type& a) const { return self_type( _mm_or_pd(m_data, a.m_data) ); }

    self_type operator == (const self_type& a) const { return self_type( _mm_cmpeq_pd (m_data, a.m_data) ); }
    self_type operator != (const self_type& a) const { return self_type( _mm_cmpneq_pd(m_data, a.m_data) ); }
    self_type operator >= (const self_type& a) const { return self_type( _mm_cmpge_pd (m_data, a.m_data) ); }
    self_type operator >  (const self_type& a) const { return self_type( _mm_cmpgt_pd (m_data, a.m_data) ); }
    self_type operator <= (const self_type& a) const { return self_type( _mm_cmple_pd (m_data, a.m_data) ); }
    self_type operator <  (const self_type& a) const { return self_type( _mm_cmplt_pd (m_data, a.m_data) ); }

    chunk_t& data() { return m_data; }
    const chunk_t& data() const { return m_data; }

  private:

    static self_type s_zero;
    static self_type s_one ;

    chunk_t m_data;

  };



#define TMPL template <>
#define TMPL_vector vector<SSE, double>

// No vector function.
  TMPL inline TMPL_vector::base_t reduce_add (const TMPL_vector& a, const uint16_t n) {
    return no_vec::reduce_add(a,n);
  }
    
// To debug, print a vector
#ifndef NDEBUG
  TMPL   inline void print(const TMPL_vector& a, const uint16_t n) {
    TMPL_vector::base_t *tmp=(TMPL_vector::base_t *)(&a.data());
    for (uint16_t i=0; i<n ; i++) std::cout<<tmp[i]<<std::endl;
  }
#endif

  TMPL   inline int movemask (const TMPL_vector& a) {
  #if defined(__evi_use_svml) || defined(__evi_use_svml_gcc)
    return _mm_movemask_pd(a.data());
  #else
    return no_vec::movemask(a);
  #endif
  }

  TMPL inline TMPL_vector ceil    (const TMPL_vector& a) {
    return TMPL_vector( _mm_ceil_pd(a.data()) );
  }

  TMPL inline TMPL_vector floor   (const TMPL_vector& a) {
    return TMPL_vector( _mm_floor_pd(a.data()) );
  }

  TMPL inline TMPL_vector max     (const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm_max_pd(a.data(), b.data()) );
  }

  TMPL inline TMPL_vector min     (const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector( _mm_min_pd(a.data(), b.data()) );
  }



  TMPL inline TMPL_vector sqrt    (const TMPL_vector& a) {
    return TMPL_vector( _mm_sqrt_pd(a.data()) );
  }

  TMPL inline TMPL_vector cbrt    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_cbrt_pd(a.data()) );
#else
    return no_vec::cbrt(a);
#endif
  }

  TMPL inline TMPL_vector inv_sqrt(const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_invsqrt_pd(a.data()) );
#else
    return TMPL_vector::one() / sqrt(a);
#endif
  }

  TMPL inline TMPL_vector inv_cbrt(const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_invcbrt_pd(a.data()) );
#else
    return TMPL_vector::one() / cbrt(a);
#endif
  }



  TMPL inline TMPL_vector exp     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_exp_pd(a.data()) );
#else
    return no_vec::exp(a);
#endif
  }

  TMPL inline TMPL_vector exp2    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_exp2_pd(a.data()) );
#else
    return no_vec::exp2(a);
#endif
  }

  TMPL inline TMPL_vector exp10   (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_exp10_pd(a.data()) );
#else
    return no_vec::exp10(a);
#endif
  }

  TMPL inline TMPL_vector log     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_log_pd(a.data()) );
#else
    return no_vec::log(a);
#endif
  }

  TMPL inline TMPL_vector log2    (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_log2_pd(a.data()) );
#else
    return no_vec::log2(a);
#endif
  }

  TMPL inline TMPL_vector log10   (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_log10_pd(a.data()) );
#else
    return no_vec::log10(a);
#endif
  }

  TMPL inline TMPL_vector pow     (const TMPL_vector& a, const TMPL_vector& b) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_pow_pd(a.data(), b.data()) );
#else
    return no_vec::pow(a, b);
#endif
  }



  TMPL inline TMPL_vector cos     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_cos_pd(a.data()) );
#else
    return no_vec::cos(a);
#endif
  }
  
  TMPL inline TMPL_vector sin     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_sin_pd(a.data()) );
#else
    return no_vec::sin(a);    
#endif
  }
  
  TMPL inline TMPL_vector tan     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_tan_pd(a.data()) );
#else
    return no_vec::tan(a);
#endif
    
  }

  TMPL inline TMPL_vector cosh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_cosh_pd(a.data()) );
#else
    return no_vec::cosh(a);
#endif
  }
  
  TMPL inline TMPL_vector sinh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_sinh_pd(a.data()) );
#else
    return no_vec::sinh(a);    
#endif
  }
  
  TMPL inline TMPL_vector tanh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_tanh_pd(a.data()) );
#else
    return no_vec::tanh(a);
#endif
  }



  TMPL inline TMPL_vector hypot   (const TMPL_vector& a, const TMPL_vector& b) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_hypot_pd(a.data(), b.data()) );
#else
    return sqrt(a*a+b*b);
#endif
  }



  TMPL inline TMPL_vector acos     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_acos_pd(a.data()) );
#else
    return no_vec::acos(a);
#endif
  }
  
  TMPL inline TMPL_vector asin     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_asin_pd(a.data()) );
#else
    return no_vec::asin(a);    
#endif
  }
  
  TMPL inline TMPL_vector atan     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_atan_pd(a.data()) );
#else
    return no_vec::atan(a);
#endif
    
  }

  TMPL inline TMPL_vector acosh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_acosh_pd(a.data()) );
#else
    return no_vec::acosh(a);
#endif
  }
  
  TMPL inline TMPL_vector asinh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_asinh_pd(a.data()) );
#else
    return no_vec::asinh(a);    
#endif
  }
  
  TMPL inline TMPL_vector atanh     (const TMPL_vector& a) {
#ifdef __evi_use_svml
    return TMPL_vector( _mm_atanh_pd(a.data()) );
#else
    return no_vec::atanh(a);
#endif
  }

  TMPL inline TMPL_vector blendv  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& mask) {
  #if defined(__evi_use_svml) || defined(__evi_use_svml_gcc)
    return TMPL_vector( _mm_blendv_pd(a.data(), b.data(), mask.data()) );
  #else
    return no_vec::blendv(a,b,mask);
  #endif
  }

  TMPL inline TMPL_vector fmadd   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return a*b + c;
#else
    return TMPL_vector( _mm_fmadd_pd(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fmsub   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return a*b - c;
#else
    return TMPL_vector( _mm_fmsub_pd(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fnmadd  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return c - a*b;
#else
    return TMPL_vector( _mm_fnmadd_pd(a.data(), b.data(), c.data()) );
#endif
  }

  TMPL inline TMPL_vector fnmsub  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
#ifdef __evi_no_fma
    return TMPL_vector( _mm_xor_pd((a*b+c).data(), _mm_set1_pd(-0.d)) );
#else
    return TMPL_vector( _mm_fnmadd_pd(a.data(), b.data(), c.data()) );
#endif
  }
  

  
#undef TMPL_vector
#undef TMPL



}
