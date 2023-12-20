#pragma once

#include <complex>

#include <onika/dac/slices.h>

namespace exaStamp
{

  struct my_complex3d
  {
    std::complex<double> x;
    std::complex<double> y;
    std::complex<double> z;
  };

  // slice identifiers can be decalred anywhere, keeing in mind it should be both easily usable
  // without polluating the namespace
  struct c3d_x_real_t {}; static inline constexpr c3d_x_real_t c3d_x_real;
  struct c3d_x_imag_t {}; static inline constexpr c3d_x_imag_t c3d_x_imag;
  struct c3d_y_real_t {}; static inline constexpr c3d_y_real_t c3d_y_real;
  struct c3d_y_imag_t {}; static inline constexpr c3d_y_imag_t c3d_y_imag;
  struct c3d_z_real_t {}; static inline constexpr c3d_z_real_t c3d_z_real;
  struct c3d_z_imag_t {}; static inline constexpr c3d_z_imag_t c3d_z_imag;    

}

namespace onika
{
  namespace dac
  {

    template<class T>
    struct ComplexRealReference
    {
      std::complex<T>& m_ref;
      inline operator T() const { return m_ref.real(); }
      inline ComplexRealReference operator = (const T& x) { m_ref.real(x); return *this; }
      inline ComplexRealReference operator += (const T& x) { m_ref.real( m_ref.real() + x ); return *this; }
    };

    template<class T>
    struct ComplexImagReference
    {
      std::complex<T>& m_ref;
      inline operator T() const { return m_ref.imag(); }
      inline ComplexImagReference operator = (const T& x) { m_ref.imag(x); return *this; }
      inline ComplexImagReference operator += (const T& x) { m_ref.imag( m_ref.imag() + x ); return *this; }
    };
    
    // 
    template<>
    struct DataSlicing<::exanb::my_complex3d>
    {
      using c3d_x_real_t = ::exanb::c3d_x_real_t;
      using c3d_x_imag_t = ::exanb::c3d_x_imag_t;
      using c3d_y_real_t = ::exanb::c3d_y_real_t;
      using c3d_y_imag_t = ::exanb::c3d_y_imag_t;
      using c3d_z_real_t = ::exanb::c3d_z_real_t;
      using c3d_z_imag_t = ::exanb::c3d_z_imag_t;
      
      using value_t = ::exanb::my_complex3d;
      using slices_t = DataSlices< c3d_x_real_t, c3d_x_imag_t, c3d_y_real_t, c3d_y_imag_t, c3d_z_real_t, c3d_z_imag_t >;

      static inline ComplexRealReference<double> get_slice_rw(value_t& v , c3d_x_real_t ) { return ComplexRealReference<double>{v.x}; }
      static inline ComplexImagReference<double> get_slice_rw(value_t& v , c3d_x_imag_t ) { return ComplexImagReference<double>{v.x}; }

      static inline ComplexRealReference<double> get_slice_rw(value_t& v , c3d_y_real_t ) { return ComplexRealReference<double>{v.y}; }
      static inline ComplexImagReference<double> get_slice_rw(value_t& v , c3d_y_imag_t ) { return ComplexImagReference<double>{v.y}; }

      static inline ComplexRealReference<double> get_slice_rw(value_t& v , c3d_z_real_t ) { return ComplexRealReference<double>{v.z}; }
      static inline ComplexImagReference<double> get_slice_rw(value_t& v , c3d_z_imag_t ) { return ComplexImagReference<double>{v.z}; }

      static inline auto get_slice_ro(const value_t& v , c3d_x_real_t ) { return v.x.real(); }
      static inline auto get_slice_ro(const value_t& v , c3d_x_imag_t ) { return v.x.imag(); }

      static inline auto get_slice_ro(const value_t& v , c3d_y_real_t ) { return v.y.real(); }
      static inline auto get_slice_ro(const value_t& v , c3d_y_imag_t ) { return v.y.imag(); }

      static inline auto get_slice_ro(const value_t& v , c3d_z_real_t ) { return v.z.real(); }
      static inline auto get_slice_ro(const value_t& v , c3d_z_imag_t ) { return v.z.imag(); }
    };

  }
}



