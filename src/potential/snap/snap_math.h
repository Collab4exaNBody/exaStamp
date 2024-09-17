#pragma once

#include <onika/cuda/cuda.h>

namespace SnapMath
{

  struct Complexd
  {
    double r;
    double i;

    ONIKA_HOST_DEVICE_FUNC inline double real() const { return r; }
    ONIKA_HOST_DEVICE_FUNC inline double imag() const { return i; }

    ONIKA_HOST_DEVICE_FUNC inline bool operator == (const Complexd& o) const { return r==o.r && i==o.i; }

    ONIKA_HOST_DEVICE_FUNC inline Complexd& operator += (const Complexd& b) { r+=b.r; i+=b.i; return *this; }
    ONIKA_HOST_DEVICE_FUNC inline Complexd& operator -= (const Complexd& b) { r-=b.r; i-=b.i; return *this; }

    ONIKA_HOST_DEVICE_FUNC inline Complexd& operator *= (const Complexd& b)
    {
      const double _r = r*b.r - i*b.i;
      const double _i = r*b.i + i*b.r;
      r=_r; i=_i;
      return *this;
    }

    ONIKA_HOST_DEVICE_FUNC inline Complexd& operator *= (double b)
    {
      r *= b;
      i *= b;
      return *this;
    }

    ONIKA_HOST_DEVICE_FUNC inline Complexd& operator /= (const Complexd& o)
    {
      const double a=r, b=i, c=o.r, d=o.i;
      const double q = c*c + d*d;
      r = (a*c+b*d) / q;
      i = (b*c-a*d) / q;
      return *this;
    }

    ONIKA_HOST_DEVICE_FUNC inline Complexd operator / (const Complexd& o) const
    {
      const double a=r, b=i, c=o.r, d=o.i;
      const double q = c*c + d*d;
      return { (a*c+b*d) / q , (b*c-a*d) / q };
    }

    ONIKA_HOST_DEVICE_FUNC inline Complexd operator - () const { return {-r,-i}; }

    ONIKA_HOST_DEVICE_FUNC inline Complexd operator * (const Complexd& b) const { return { r*b.r - i*b.i , r*b.i + i*b.r }; }
    ONIKA_HOST_DEVICE_FUNC inline Complexd operator - (const Complexd& b) const { return { r-b.r , i-b.i }; }
    ONIKA_HOST_DEVICE_FUNC inline Complexd operator + (const Complexd& b) const { return { r+b.r , i+b.i }; }

    ONIKA_HOST_DEVICE_FUNC inline Complexd operator + (double b) const { return { r+b , i }; }
    ONIKA_HOST_DEVICE_FUNC inline Complexd operator - (double b) const { return { r-b , i }; }
    ONIKA_HOST_DEVICE_FUNC inline Complexd operator * (double b) const { return { r*b , i*b }; }
    ONIKA_HOST_DEVICE_FUNC inline Complexd operator / (double b) const { return (*this) / Complexd{b,0.0}; }
  };

  ONIKA_HOST_DEVICE_FUNC inline double norm2(const Complexd& c) { return c.r*c.r + c.i*c.i ; }
  ONIKA_HOST_DEVICE_FUNC inline double norm(const Complexd& c) { return sqrt( norm2(c) ); }

  ONIKA_HOST_DEVICE_FUNC inline double real(const Complexd& c) { return c.r; }
  ONIKA_HOST_DEVICE_FUNC inline double imag(const Complexd& c) { return c.i; }

  ONIKA_HOST_DEVICE_FUNC inline Complexd operator + (double a, const Complexd &b ) { return { a+b.r , b.i }; }
  ONIKA_HOST_DEVICE_FUNC inline Complexd operator - (double a, const Complexd &b ) { return { a-b.r , -b.i }; }
  ONIKA_HOST_DEVICE_FUNC inline Complexd operator * (double a, const Complexd &b ) { return b * a; }
  ONIKA_HOST_DEVICE_FUNC inline Complexd operator / (double a, const Complexd &b ) { return Complexd{a,0.0} / b; }

  struct Double3d
  {
	  double x;
	  double y;
	  double z;
	  
    ONIKA_HOST_DEVICE_FUNC inline bool operator == (const Double3d& o) const { return x==o.x && y==o.y && z==o.z; }

#   define BINOP(op) \
      ONIKA_HOST_DEVICE_FUNC inline Double3d operator op ( const double& r ) const { return { x op r , y op r , z op r }; } \
      ONIKA_HOST_DEVICE_FUNC inline Double3d operator op ( const Double3d& r ) const { return { x op r.x , y op r.y , z op r.z }; }
    BINOP(+)
    BINOP(-)
    BINOP(*)
    BINOP(/)
#   undef BINOP
  };

  struct Complex3d
  {
	  Complexd x;
	  Complexd y;
	  Complexd z;

    ONIKA_HOST_DEVICE_FUNC inline bool operator == (const Complex3d& o) const { return x==o.x && y==o.y && z==o.z; }

    ONIKA_HOST_DEVICE_FUNC inline Complex3d operator - () const { return { -x , -y , -z }; }

#   define BINOP(op) \
      ONIKA_HOST_DEVICE_FUNC inline Complex3d operator op ( const double   & r ) const { return { x op r   , y op r   , z op r   }; } \
      ONIKA_HOST_DEVICE_FUNC inline Complex3d operator op ( const Complexd & r ) const { return { x op r   , y op r   , z op r   }; } \
      ONIKA_HOST_DEVICE_FUNC inline Complex3d operator op ( const Double3d & r ) const { return { x op r.x , y op r.y , z op r.z }; } \
      ONIKA_HOST_DEVICE_FUNC inline Complex3d operator op ( const Complex3d& r ) const { return { x op r.x , y op r.y , z op r.z }; }
    BINOP(+)
    BINOP(-)
    BINOP(*)
    BINOP(/)
#   undef BINOP
  };

  ONIKA_HOST_DEVICE_FUNC inline double norm(const Complex3d& c) { return sqrt( norm2(c.x) + norm2(c.y) + norm2(c.y) ); }

  ONIKA_HOST_DEVICE_FUNC inline Complexd conj( const Complexd &b ) { return { b.r , -b.i }; }
  ONIKA_HOST_DEVICE_FUNC inline Complex3d conj( const Complex3d &b ) { return { conj(b.x) , conj(b.y) , conj(b.z) }; }

# define EXTERNAL_BINOP(op) \
    ONIKA_HOST_DEVICE_FUNC inline Complex3d operator op ( const double  & l , const Complex3d& r ) { return { l   op r.x , l   op r.y , l   op r.z }; } \
    ONIKA_HOST_DEVICE_FUNC inline Complex3d operator op ( const Complexd& l , const Complex3d& r ) { return { l   op r.x , l   op r.y , l   op r.z }; } \
    ONIKA_HOST_DEVICE_FUNC inline Complex3d operator op ( const Double3d& l , const Complex3d& r ) { return { l.x op r.x , l.y op r.y , l.z op r.z }; } \
    ONIKA_HOST_DEVICE_FUNC inline Double3d  operator op ( const double  & l , const Double3d&  r ) { return { l   op r.x , l   op r.y , l   op r.z }; } \
    ONIKA_HOST_DEVICE_FUNC inline Complex3d operator op ( const Complexd& l , const Double3d&  r ) { return { l   op r.x , l   op r.y , l   op r.z }; } 
  EXTERNAL_BINOP(+)
  EXTERNAL_BINOP(-)
  EXTERNAL_BINOP(*)
  EXTERNAL_BINOP(/)  
# undef EXTERNAL_BINOP


  inline std::ostream& operator << ( std::ostream& out , const Complexd& o )
  {
    return out << "("<<o.r<<","<<o.i<<")";
  }

  inline std::ostream& operator << ( std::ostream& out , const Double3d& o )
  {
    return out << "{"<<o.x<<","<<o.y<<","<<o.z<<"}";
  }

  inline std::ostream& operator << ( std::ostream& out , const Complex3d& o )
  {
    return out << "{"<<o.x<<","<<o.y<<","<<o.z<<"}";
  }

}

