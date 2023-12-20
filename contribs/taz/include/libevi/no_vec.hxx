/// @file
/// @brief

#ifdef __evi_use_svml
namespace simd {

  namespace no_vec {

#define TMPL template <instr_t instr, class base_t> 
#define TMPL_vector vector<instr, base_t>

  TMPL inline typename TMPL_vector::base_t reduce_add (const TMPL_vector& a, const uint16_t n) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    typename TMPL_vector::base_t tmp1=0;
    for (uint16_t i=0; i<n ; i++){tmp1 += tmp[i];}
    return tmp1;
  }

#undef TMPL_vector
#undef TMPL

  }

}

#endif

#ifndef __evi_use_svml



namespace simd {

  namespace no_vec {

#define TMPL template <instr_t instr, class base_t> 
#define TMPL_vector vector<instr, base_t>

  TMPL inline typename TMPL_vector::base_t reduce_add (const TMPL_vector& a, const uint16_t n) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    typename TMPL_vector::base_t tmp1=0;
    for (uint16_t i=0; i<n ; i++){tmp1 += tmp[i];}
    return tmp1;
  }

#ifndef NDEBUG
  TMPL inline void print(const TMPL_vector& a, const uint16_t n) {
    typename TMPL_vector::base_t *tmp=(typename TMPL_vector::base_t *)(&a.data());
    typename TMPL_vector::base_t tmp1=0;
    for (uint16_t i=0; i<n ; i++){std::cout<<tmp[i]<<std::endl;}
  }
#endif

  TMPL   inline int movemask (const TMPL_vector& a) {;
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size ; i++){if(tmp[i]!=0) return 1;}
    return 0;
  }

  TMPL inline TMPL_vector blendv (const TMPL_vector& a, const TMPL_vector& b, const TMPL_vector& mask) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmpa [TMPL_vector::size];
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmpb [TMPL_vector::size];
    alignas(TMPL_vector::align) typename TMPL_vector::base_t val [TMPL_vector::size];

    a.store(tmpa);
    b.store(tmpb);
    mask.store(val);

    for (uint16_t i=0; i<TMPL_vector::size; ++i){ if(val[i]==0) tmp[i] = tmpa[i]; else tmp[i] = tmpb[i] ;}
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector cbrt (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::cbrt(tmp[i]);
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector exp (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::exp(tmp[i]);
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector exp2 (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::exp2(tmp[i]);
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector exp10 (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::pow (10.f, tmp[i]);
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector log (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::log(tmp[i]);
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector log2 (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::log2(tmp[i]);
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector log10 (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::log10(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector pow (const TMPL_vector& a, const TMPL_vector& b) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmpA [TMPL_vector::size];
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmpB [TMPL_vector::size];
    a.store(tmpA);
    b.store(tmpB);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmpA[i] = ::pow(tmpA[i], tmpB[i]);
    return TMPL_vector(tmpA);
  }
  
  TMPL inline TMPL_vector cos (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::cos(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector sin (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::sin(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector tan (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::tan(tmp[i]);
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector cosh (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::cosh(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector sinh (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::sinh(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector tanh (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::tanh(tmp[i]);
    return TMPL_vector(tmp);
  }


  
  TMPL inline TMPL_vector acos (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::acos(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector asin (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::asin(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector atan (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::atan(tmp[i]);
    return TMPL_vector(tmp);
  }

  TMPL inline TMPL_vector acosh (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::acosh(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector asinh (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::asinh(tmp[i]);
    return TMPL_vector(tmp);
  }
  
  TMPL inline TMPL_vector atanh (const TMPL_vector& a) {
    alignas(TMPL_vector::align) typename TMPL_vector::base_t tmp [TMPL_vector::size];
    a.store(tmp);
    for (uint16_t i=0; i<TMPL_vector::size; ++i) tmp[i] = ::atanh(tmp[i]);
    return TMPL_vector(tmp);
  }

#undef TMPL_vector
#undef TMPL

  }

}



#endif
