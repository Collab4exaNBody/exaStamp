/// @file
/// @brief Vectorization using intrinsics



namespace simd {

  /// @brief Specialization of simd::vector for SER floats
  template <> 
  class vector<SER, float> {

  public:



    typedef vector<SER, float> self_type;
    typedef float base_t;
    typedef float chunk_t;


    static constexpr uint16_t size  = ser::size;
    static constexpr uint16_t align = ser::align;

    static constexpr uint16_t chunk_size = get_chunk_size<base_t, size>();

    static constexpr self_type& zero() { return s_zero; };
    static constexpr self_type& one () { return s_one ; };





    vector() {}
    ~vector() {}

    vector(const base_t& in) : m_data( in) {}
    vector(const base_t* in) : m_data(*in) {}



    vector(const self_type& a) : m_data(a.m_data) {}

    self_type& operator = (const self_type& a) {
      m_data = a.m_data;
      return *this;
    }

    void load (const base_t* in) { m_data = *in; }

    void store(base_t* out) const { *out = m_data; }

    self_type operator + (const self_type& a) const { return self_type(m_data + a.m_data); }
    self_type operator - (const self_type& a) const { return self_type(m_data - a.m_data); }
    self_type operator * (const self_type& a) const { return self_type(m_data * a.m_data); }
    self_type operator / (const self_type& a) const { return self_type(m_data / a.m_data); }

    self_type operator && (const self_type& a) const { return self_type(m_data && a.m_data); }
    self_type operator || (const self_type& a) const { return self_type(m_data || a.m_data); }

    self_type operator == (const self_type& a) const { return m_data == a.m_data ? self_type(0xFFFFFFFF) : s_zero; }
    self_type operator != (const self_type& a) const { return m_data != a.m_data ? self_type(0xFFFFFFFF) : s_zero; }
    self_type operator >= (const self_type& a) const { return m_data >= a.m_data ? self_type(0xFFFFFFFF) : s_zero; }
    self_type operator >  (const self_type& a) const { return m_data >  a.m_data ? self_type(0xFFFFFFFF) : s_zero; }
    self_type operator <= (const self_type& a) const { return m_data <= a.m_data ? self_type(0xFFFFFFFF) : s_zero; }
    self_type operator <  (const self_type& a) const { return m_data <  a.m_data ? self_type(0xFFFFFFFF) : s_zero; }

    chunk_t& data() { return m_data; }
    const chunk_t& data() const { return m_data; }

  private:

    static self_type s_zero;
    static self_type s_one ;

    chunk_t m_data;

  };



#define TMPL template <>
#define TMPL_vector vector<SER, float>

//For MEAM
// No vector function.
  TMPL inline TMPL_vector::base_t reduce_add (const TMPL_vector& a, const uint16_t n) {
    TMPL_vector::base_t *tmp=(TMPL_vector::base_t *)(&a.data());
    TMPL_vector::base_t tmp1=0;
    for (uint16_t i=0; i<n ; i++) tmp1 += tmp[i];
    return tmp1;
  }
    
//for debug, print a vector
#ifndef NDEBUG
  TMPL   inline void print(const TMPL_vector& a, const uint16_t n) {
    TMPL_vector::base_t *tmp=(TMPL_vector::base_t *)(&a.data());
    for (uint16_t i=0; i<n ; i++) std::cout<<tmp[i]<<std::endl;
  }
#endif

  TMPL   inline int movemask (const TMPL_vector& a) {
    if((TMPL_vector::base_t)(a.data()) == 0xFFFFFFFF) return 1;
    else return 0;
  }


  TMPL inline TMPL_vector ceil    (const TMPL_vector& a) { 
    return TMPL_vector(::ceil(a.data()));
  }

  TMPL inline TMPL_vector floor   (const TMPL_vector& a) { 
    return TMPL_vector(::floor(a.data()));
  }

  TMPL inline TMPL_vector max     (const TMPL_vector& a, const TMPL_vector& b) { 
    return TMPL_vector(::fmax(a.data(), b.data()));
  }

  TMPL inline TMPL_vector min     (const TMPL_vector& a, const TMPL_vector& b) { 
    return TMPL_vector(::fmin(a.data(), b.data()));
  }



  TMPL inline TMPL_vector sqrt    (const TMPL_vector& a) { 
    return TMPL_vector(::sqrt(a.data()));
  }

  TMPL inline TMPL_vector cbrt    (const TMPL_vector& a) { 
    return TMPL_vector(::cbrt(a.data()));
  }

  TMPL inline TMPL_vector inv_sqrt(const TMPL_vector& a) { 
    return TMPL_vector::one() / sqrt(a);
  }

  TMPL inline TMPL_vector inv_cbrt(const TMPL_vector& a) { 
    return TMPL_vector::one() / cbrt(a);
  }



  TMPL inline TMPL_vector exp     (const TMPL_vector& a) { 
    return TMPL_vector(::exp (a.data()));
  }

  TMPL inline TMPL_vector exp2    (const TMPL_vector& a) { 
    return TMPL_vector(::exp2(a.data()));
  }

  TMPL inline TMPL_vector exp10   (const TMPL_vector& a) { 
    return TMPL_vector(::pow (10.f, a.data()));
  }

  TMPL inline TMPL_vector log     (const TMPL_vector& a) { 
    return TMPL_vector(::log  (a.data()));
  }

  TMPL inline TMPL_vector log2    (const TMPL_vector& a) { 
    return TMPL_vector(::log2 (a.data()));
  }

  TMPL inline TMPL_vector log10   (const TMPL_vector& a) { 
    return TMPL_vector(::log10(a.data()));
  }

  TMPL inline TMPL_vector pow     (const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector(::pow(a.data(), b.data()));
  }



  TMPL inline TMPL_vector cos     (const TMPL_vector& a) { 
    return TMPL_vector(::cos (a.data()));
  }

  TMPL inline TMPL_vector sin     (const TMPL_vector& a) { 
    return TMPL_vector(::sin (a.data()));
  }

  TMPL inline TMPL_vector tan     (const TMPL_vector& a) { 
    return TMPL_vector(::tan (a.data()));
  }

  TMPL inline TMPL_vector cosh    (const TMPL_vector& a) { 
    return TMPL_vector(::cosh(a.data()));
  }

  TMPL inline TMPL_vector sinh    (const TMPL_vector& a) { 
    return TMPL_vector(::sinh(a.data()));
  }

  TMPL inline TMPL_vector tanh    (const TMPL_vector& a) { 
    return TMPL_vector(::tanh(a.data()));
  }



  TMPL inline TMPL_vector hypot   (const TMPL_vector& a, const TMPL_vector& b) {
    return TMPL_vector(::hypot(a.data(), b.data())); 
  }



  TMPL inline TMPL_vector acos    (const TMPL_vector& a) { 
    return TMPL_vector(::acosh(a.data()));
  }

  TMPL inline TMPL_vector asin    (const TMPL_vector& a) { 
    return TMPL_vector(::asinh(a.data()));
  }

  TMPL inline TMPL_vector atan    (const TMPL_vector& a) { 
    return TMPL_vector(::atanh(a.data()));
  }

  TMPL inline TMPL_vector acosh   (const TMPL_vector& a) { 
    return TMPL_vector(::acosh(a.data()));
  }

  TMPL inline TMPL_vector asinh   (const TMPL_vector& a) { 
    return TMPL_vector(::asinh(a.data()));
  }

  TMPL inline TMPL_vector atanh   (const TMPL_vector& a) { 
    return TMPL_vector(::atanh(a.data()));
  }



  TMPL inline TMPL_vector blendv  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& mask) { 
    return TMPL_vector( mask.data() ? b.data() : a.data() );
  }



  TMPL inline TMPL_vector fmadd   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
    return TMPL_vector(::fma(a.data(), b.data(), c.data()));
  }

  TMPL inline TMPL_vector fmsub   (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
    return TMPL_vector(::fma(a.data(), b.data(), -c.data()));
  }

  TMPL inline TMPL_vector fnmadd  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
    return TMPL_vector(::fma(-a.data(), b.data(), c.data()));
  }

  TMPL inline TMPL_vector fnmsub  (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& c) {
    return TMPL_vector(::fma(-a.data(), b.data(), -c.data()));
  }



#undef TMPL_vector
#undef TMPL



}
