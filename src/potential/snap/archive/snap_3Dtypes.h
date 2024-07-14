#pragma once

#include <onika/cuda/cuda.h>

namespace SnapExt
{

struct Complexd
{
  double r;
  double i;

  Complexd(const Complexd&) = default;
  Complexd(Complexd&&) = default;
  ONIKA_HOST_DEVICE_FUNC inline constexpr Complexd(double _r=0.0, double _i=0.0) : r(_r), i(_i) {}

  Complexd& operator = (const Complexd& rhs) = default;
  Complexd& operator = (Complexd&& rhs) = default;

  ONIKA_HOST_DEVICE_FUNC inline double real() const { return r; }

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
    return Complexd( (a*c+b*d) / q , (b*c-a*d) / q );
  }

  ONIKA_HOST_DEVICE_FUNC inline Complexd operator - () const { return Complexd(-r,-i); }

  ONIKA_HOST_DEVICE_FUNC inline Complexd operator * (const Complexd& b) const { return Complexd( r*b.r - i*b.i , r*b.i + i*b.r ); }
  ONIKA_HOST_DEVICE_FUNC inline Complexd operator - (const Complexd& b) const { return Complexd( r-b.r , i-b.i ); }
  ONIKA_HOST_DEVICE_FUNC inline Complexd operator + (const Complexd& b) const { return Complexd( r+b.r , i+b.i ); }

  ONIKA_HOST_DEVICE_FUNC inline Complexd operator + (double b) const { return Complexd( r+b , i ); }
  ONIKA_HOST_DEVICE_FUNC inline Complexd operator - (double b) const { return Complexd( r-b , i ); }
  ONIKA_HOST_DEVICE_FUNC inline Complexd operator * (double b) const { return Complexd(r*b,i*b); }
  ONIKA_HOST_DEVICE_FUNC inline Complexd operator / (double b) const { return (*this) / Complexd(b); }
};

struct double3d
{
	double x;
	double y;
	double z;
  
	ONIKA_HOST_DEVICE_FUNC constexpr inline double3d() : x(0), y(0), z(0) {}
	double3d(double3d const &d) = default;
	ONIKA_HOST_DEVICE_FUNC inline double3d(double a) : x(a), y(a), z(a) {}
	ONIKA_HOST_DEVICE_FUNC inline double3d(double x, double y, double z) : x(x), y(y), z(z) {}
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator=(double d) {
		x = d;
		y = d;
		z = d;
		return *this;
	}
	double3d & operator=(double3d const &d) = default;
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator+=(double d) {
		x += d;
		y += d;
		z += d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator+=(double3d const &d) {
		x += d.x;
		y += d.y;
		z += d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator-=(double d) {
		x -= d;
		y -= d;
		z -= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator-=(double3d const &d) {
		x -= d.x;
		y -= d.y;
		z -= d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator*=(double d) {
		x *= d;
		y *= d;
		z *= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator*=(double3d const &d) {
		x *= d.x;
		y *= d.y;
		z *= d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator/=(double d) {
		x /= d;
		y /= d;
		z /= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline double3d & operator/=(double3d const &d) {
		x /= d.x;
		y /= d.y;
		z /= d.z;

		return *this;
	}
};

struct complex3d
{
	Complexd x;
	Complexd y;
	Complexd z;
	complex3d() = default;
	complex3d(const complex3d&) = default;
	ONIKA_HOST_DEVICE_FUNC constexpr inline complex3d(double d) : x(d), y(d), z(d) {}
	ONIKA_HOST_DEVICE_FUNC constexpr inline complex3d(double x, double y, double z) : x(x), y(y), z(z) {}
	ONIKA_HOST_DEVICE_FUNC constexpr inline complex3d(double3d const&a) : x(a.x), y(a.y), z(a.z) {}
	ONIKA_HOST_DEVICE_FUNC constexpr inline complex3d(Complexd const &c) : x(c), y(c), z(c) {}
	ONIKA_HOST_DEVICE_FUNC constexpr inline complex3d(Complexd const & x, Complexd const & y, Complexd const & z) : x(x), y(y), z(z) {}

	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator=(double d) {
		x = d;
		y = d;
		z = d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator=(double3d const &d) {
		x = d.x;
		y = d.y;
		z = d.z;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator=(Complexd   const &d) {
		x = d;
		y = d;
		z = d;
		return *this;
	}
	complex3d & operator=(complex3d const &d) = default;
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator+=(double d) {
		x += d;
		y += d;
		z += d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator+=(double3d const &d) {
		x += d.x;
		y += d.y;
		z += d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator+=(Complexd   const &d) {
		x += d;
		y += d;
		z += d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator+=(complex3d const &d) {
		x += d.x;
		y += d.y;
		z += d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator-=(double d) {
		x -= d;
		y -= d;
		z -= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator-=(double3d const &d) {
		x -= d.x;
		y -= d.y;
		z -= d.z;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator-=(Complexd   const &d) {
		x -= d;
		y -= d;
		z -= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator-=(complex3d const &d) {
		x -= d.x;
		y -= d.y;
		z -= d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator*=(double d) {
		x *= d;
		y *= d;
		z *= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator*=(double3d const &d) {
		x *= d.x;
		y *= d.y;
		z *= d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator*=(Complexd   const &d) {
		x *= d;
		y *= d;
		z *= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator*=(complex3d const &d) {
		x *= d.x;
		y *= d.y;
		z *= d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator/=(double d) {
		x /= d;
		y /= d;
		z /= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator/=(double3d const &d) {
		x /= d.x;
		y /= d.y;
		z /= d.z;

		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator/=(Complexd   const &d) {
		x /= d;
		y /= d;
		z /= d;
		return *this;
	}
	ONIKA_HOST_DEVICE_FUNC inline complex3d & operator/=(complex3d const &d) {
		x /= d.x;
		y /= d.y;
		z /= d.z;
		return *this;
	}
};

ONIKA_HOST_DEVICE_FUNC inline Complexd operator + (double a, const Complexd &b ) { return Complexd( a+b.r , b.i ); }
ONIKA_HOST_DEVICE_FUNC inline Complexd operator - (double a, const Complexd &b ) { return Complexd( a-b.r , -b.i ); }
ONIKA_HOST_DEVICE_FUNC inline Complexd operator * (double a, const Complexd &b ) { return b * a; }
ONIKA_HOST_DEVICE_FUNC inline Complexd operator / (double a, const Complexd &b ) { return Complexd(a) / b; }

ONIKA_HOST_DEVICE_FUNC inline Complexd conj( const Complexd &b ) { return Complexd(b.r,-b.i); }

// double3d binary operators
// addition
ONIKA_HOST_DEVICE_FUNC inline double3d operator+(double a, double3d const &b) {
	double3d r = b;
	r += a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  double3d operator+(double3d const &a, double b) {
	double3d r = a;
	r += b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  double3d operator+(double3d const &a, double3d const &b) {
	double3d r = a;
	r += b;
	return r;
}

// unary minus
ONIKA_HOST_DEVICE_FUNC inline  double3d operator-(double3d const &a) {
	double3d r;
	r.x = -a.x;
	r.y = -a.y;
	r.z = -a.z;
	return r;
}

// substraction
ONIKA_HOST_DEVICE_FUNC inline  double3d operator-(double a, double3d const &b) {
	double3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  double3d operator-(double3d const &a, double b) {
	double3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  double3d operator-(double3d const &a, double3d const &b) {
	double3d r = a;
	r -= b;
	return r;
}

// multiplication
ONIKA_HOST_DEVICE_FUNC inline  double3d operator*(double a, double3d const &b) {
	double3d r = b;
	r *= a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  double3d operator*(double3d const &a, double b) {
	double3d r = a;
	r *= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  double3d operator*(double3d const &a, double3d const &b) {
	double3d r = a;
	r *= b;
	return r;
}

// division
ONIKA_HOST_DEVICE_FUNC inline  double3d operator/(double a, double3d const &b) {
	double3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  double3d operator/(double3d const &a, double b) {
	double3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  double3d operator/(double3d const &a, double3d const &b) {
	double3d r = a;
	r /= b;
	return r;
}

// complex3d operators
// addition

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(double a, complex3d const &b) {
	complex3d r = b;
	r += a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(complex3d const &a, double b) {
	complex3d r = a;
	r += b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(Complexd  const &a, double3d const &b) {
	complex3d r = b;
	r += a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(double3d const &a, Complexd  const &b) {
	complex3d r = a;
	r += b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(Complexd  const &a, complex3d const &b) {
	complex3d r = b;
	r += a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(complex3d const &a, Complexd  const &b) {
	complex3d r = a;
	r += b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(double3d const & a, complex3d const &b) {
	complex3d r = b;
	r += a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(complex3d const &a, double3d const & b) {
	complex3d r = a;
	r += b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator+(complex3d const &a, complex3d const &b) {
	complex3d r = a;
	r += b;
	return r;
}

// unary minus
ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(complex3d const &a) { return complex3d(-a.x , -a.y , -a.z); }

// substraction
ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(double a, complex3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(complex3d const &a, double b) {
	complex3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(Complexd  const &a, double3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(double3d const &a, Complexd  const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(Complexd  const &a, complex3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(complex3d const &a, Complexd  const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(double3d const & a, complex3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(complex3d const &a, double3d const & b) {
	complex3d r = a;
	r -= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator-(complex3d const &a, complex3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

// multiplication
ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(double a, complex3d const &b) {
	complex3d r = b;
	r *= a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(complex3d const &a, double b) {
	complex3d r = a;
	r *= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(Complexd  const &a, double3d const &b) {
	complex3d r = b;
	r *= a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(double3d const &a, Complexd  const &b) {
	complex3d r = a;
	r *= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(Complexd  const &a, complex3d const &b) {
	complex3d r = b;
	r *= a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(complex3d const &a, Complexd  const &b) {
	complex3d r = a;
	r *= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(double3d const & a, complex3d const &b) {
	complex3d r = b;
	r *= a;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(complex3d const &a, double3d const & b) {
	complex3d r = a;
	r *= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator*(complex3d const &a, complex3d const &b) {
	complex3d r = a;
	r *= b;
	return r;
}

// division
ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(double a, complex3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(complex3d const &a, double b) {
	complex3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(Complexd  const &a, double3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(double3d const &a, Complexd  const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(Complexd  const &a, complex3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(complex3d const &a, Complexd  const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(double3d const & a, complex3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(complex3d const &a, double3d const & b) {
	complex3d r = a;
	r /= b;
	return r;
}

ONIKA_HOST_DEVICE_FUNC inline  complex3d operator/(complex3d const &a, complex3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

}

