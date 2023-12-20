#ifndef SNAP_3D_TYPES
#define SNAP_3D_TYPES
#include <complex>
using std::complex;

struct double3d
{
	double x;
	double y;
	double z;
	double3d() : x(0), y(0), z(0) {}
	double3d(double a) : x(a), y(a), z(a) {}
	double3d(double x, double y, double z) : x(x), y(y), z(z) {}
	double3d & operator=(double d) {
		this->x = d;
		this->y = d;
		this->z = d;
		return *this;
	}
	double3d & operator=(double3d const &d) = default;
	double3d & operator+=(double d) {
		this->x += d;
		this->y += d;
		this->z += d;
		return *this;
	}
	double3d & operator+=(double3d const &d) {
		this->x += d.x;
		this->y += d.y;
		this->z += d.z;

		return *this;
	}
	double3d & operator-=(double d) {
		this->x -= d;
		this->y -= d;
		this->z -= d;
		return *this;
	}
	double3d & operator-=(double3d const &d) {
		this->x -= d.x;
		this->y -= d.y;
		this->z -= d.z;

		return *this;
	}
	double3d & operator*=(double d) {
		this->x *= d;
		this->y *= d;
		this->z *= d;
		return *this;
	}
	double3d & operator*=(double3d const &d) {
		this->x *= d.x;
		this->y *= d.y;
		this->z *= d.z;

		return *this;
	}
	double3d & operator/=(double d) {
		this->x /= d;
		this->y /= d;
		this->z /= d;
		return *this;
	}
	double3d & operator/=(double3d const &d) {
		this->x /= d.x;
		this->y /= d.y;
		this->z /= d.z;

		return *this;
	}
};

struct complex3d
{
	complex<double> x;
	complex<double> y;
	complex<double> z;
	complex3d() = default;
	complex3d(double d) : x(d), y(d), z(d) {}
	complex3d(double x, double y, double z) : x(x), y(y), z(z) {}
	complex3d(double3d const&a) : x(a.x), y(a.y), z(a.z) {}
	complex3d(complex<double>const &c) : x(c), y(c), z(c) {}
	complex3d(complex<double>const & x, complex<double>const & y, complex<double>const & z) : x(x), y(y), z(z) {}
	complex3d & operator=(double d) {
		this->x = d;
		this->y = d;
		this->z = d;
		return *this;
	}
	complex3d & operator=(double3d const &d) {
		this->x = d.x;
		this->y = d.y;
		this->z = d.z;
		return *this;
	}
	complex3d & operator=(complex<double>  const &d) {
		this->x = d;
		this->y = d;
		this->z = d;
		return *this;
	}
	complex3d & operator=(complex3d const &d) = default;
	complex3d & operator+=(double d) {
		this->x += d;
		this->y += d;
		this->z += d;
		return *this;
	}
	complex3d & operator+=(double3d const &d) {
		this->x += d.x;
		this->y += d.y;
		this->z += d.z;

		return *this;
	}
	complex3d & operator+=(complex<double>  const &d) {
		this->x += d;
		this->y += d;
		this->z += d;
		return *this;
	}
	complex3d & operator+=(complex3d const &d) {
		this->x += d.x;
		this->y += d.y;
		this->z += d.z;

		return *this;
	}
	complex3d & operator-=(double d) {
		this->x -= d;
		this->y -= d;
		this->z -= d;
		return *this;
	}
	complex3d & operator-=(double3d const &d) {
		this->x -= d.x;
		this->y -= d.y;
		this->z -= d.z;

		return *this;
	}
	complex3d & operator-=(complex<double>  const &d) {
		this->x -= d;
		this->y -= d;
		this->z -= d;
		return *this;
	}
	complex3d & operator-=(complex3d const &d) {
		this->x -= d.x;
		this->y -= d.y;
		this->z -= d.z;

		return *this;
	}
	complex3d & operator*=(double d) {
		this->x *= d;
		this->y *= d;
		this->z *= d;
		return *this;
	}
	complex3d & operator*=(double3d const &d) {
		this->x *= d.x;
		this->y *= d.y;
		this->z *= d.z;

		return *this;
	}
	complex3d & operator*=(complex<double>  const &d) {
		this->x *= d;
		this->y *= d;
		this->z *= d;
		return *this;
	}
	complex3d & operator*=(complex3d const &d) {
		this->x *= d.x;
		this->y *= d.y;
		this->z *= d.z;

		return *this;
	}
	complex3d & operator/=(double d) {
		this->x /= d;
		this->y /= d;
		this->z /= d;
		return *this;
	}
	complex3d & operator/=(double3d const &d) {
		this->x /= d.x;
		this->y /= d.y;
		this->z /= d.z;

		return *this;
	}
	complex3d & operator/=(complex<double>  const &d) {
		this->x /= d;
		this->y /= d;
		this->z /= d;
		return *this;
	}
	complex3d & operator/=(complex3d const &d) {
		this->x /= d.x;
		this->y /= d.y;
		this->z /= d.z;

		return *this;
	}
};

// double3d binary operators
// addition
inline double3d operator+(double a, double3d const &b) {
	double3d r = b;
	r += a;
	return r;
}

inline double3d operator+(double3d const &a, double b) {
	double3d r = a;
	r += b;
	return r;
}

inline double3d operator+(double3d const &a, double3d const &b) {
	double3d r = a;
	r += b;
	return r;
}

// unary minus
inline double3d operator-(double3d const &a) {
	double3d r;
	r.x = -a.x;
	r.y = -a.y;
	r.z = -a.z;
	return r;
}

// substraction
inline double3d operator-(double a, double3d const &b) {
	double3d r = a;
	r -= b;
	return r;
}

inline double3d operator-(double3d const &a, double b) {
	double3d r = a;
	r -= b;
	return r;
}

inline double3d operator-(double3d const &a, double3d const &b) {
	double3d r = a;
	r -= b;
	return r;
}

// multiplication
inline double3d operator*(double a, double3d const &b) {
	double3d r = b;
	r *= a;
	return r;
}

inline double3d operator*(double3d const &a, double b) {
	double3d r = a;
	r *= b;
	return r;
}

inline double3d operator*(double3d const &a, double3d const &b) {
	double3d r = a;
	r *= b;
	return r;
}

// division
inline double3d operator/(double a, double3d const &b) {
	double3d r = a;
	r /= b;
	return r;
}

inline double3d operator/(double3d const &a, double b) {
	double3d r = a;
	r /= b;
	return r;
}

inline double3d operator/(double3d const &a, double3d const &b) {
	double3d r = a;
	r /= b;
	return r;
}

// complex3d operators
// addition

inline complex3d operator+(double a, complex3d const &b) {
	complex3d r = b;
	r += a;
	return r;
}

inline complex3d operator+(complex3d const &a, double b) {
	complex3d r = a;
	r += b;
	return r;
}

inline complex3d operator+(complex<double> const &a, double3d const &b) {
	complex3d r = b;
	r += a;
	return r;
}

inline complex3d operator+(double3d const &a, complex<double> const &b) {
	complex3d r = a;
	r += b;
	return r;
}

inline complex3d operator+(complex<double> const &a, complex3d const &b) {
	complex3d r = b;
	r += a;
	return r;
}

inline complex3d operator+(complex3d const &a, complex<double> const &b) {
	complex3d r = a;
	r += b;
	return r;
}

inline complex3d operator+(double3d const & a, complex3d const &b) {
	complex3d r = b;
	r += a;
	return r;
}

inline complex3d operator+(complex3d const &a, double3d const & b) {
	complex3d r = a;
	r += b;
	return r;
}

inline complex3d operator+(complex3d const &a, complex3d const &b) {
	complex3d r = a;
	r += b;
	return r;
}

// unary minus
inline complex3d operator-(complex3d const &a) {
	complex3d r;
	r.x = -a.x;
	r.y = -a.y;
	r.z = -a.z;
	return r;
}

// substraction
inline complex3d operator-(double a, complex3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

inline complex3d operator-(complex3d const &a, double b) {
	complex3d r = a;
	r -= b;
	return r;
}

inline complex3d operator-(complex<double> const &a, double3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

inline complex3d operator-(double3d const &a, complex<double> const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

inline complex3d operator-(complex<double> const &a, complex3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

inline complex3d operator-(complex3d const &a, complex<double> const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

inline complex3d operator-(double3d const & a, complex3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

inline complex3d operator-(complex3d const &a, double3d const & b) {
	complex3d r = a;
	r -= b;
	return r;
}

inline complex3d operator-(complex3d const &a, complex3d const &b) {
	complex3d r = a;
	r -= b;
	return r;
}

// multiplication
inline complex3d operator*(double a, complex3d const &b) {
	complex3d r = b;
	r *= a;
	return r;
}

inline complex3d operator*(complex3d const &a, double b) {
	complex3d r = a;
	r *= b;
	return r;
}

inline complex3d operator*(complex<double> const &a, double3d const &b) {
	complex3d r = b;
	r *= a;
	return r;
}

inline complex3d operator*(double3d const &a, complex<double> const &b) {
	complex3d r = a;
	r *= b;
	return r;
}

inline complex3d operator*(complex<double> const &a, complex3d const &b) {
	complex3d r = b;
	r *= a;
	return r;
}

inline complex3d operator*(complex3d const &a, complex<double> const &b) {
	complex3d r = a;
	r *= b;
	return r;
}

inline complex3d operator*(double3d const & a, complex3d const &b) {
	complex3d r = b;
	r *= a;
	return r;
}

inline complex3d operator*(complex3d const &a, double3d const & b) {
	complex3d r = a;
	r *= b;
	return r;
}

inline complex3d operator*(complex3d const &a, complex3d const &b) {
	complex3d r = a;
	r *= b;
	return r;
}

// division
inline complex3d operator/(double a, complex3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

inline complex3d operator/(complex3d const &a, double b) {
	complex3d r = a;
	r /= b;
	return r;
}

inline complex3d operator/(complex<double> const &a, double3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

inline complex3d operator/(double3d const &a, complex<double> const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

inline complex3d operator/(complex<double> const &a, complex3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

inline complex3d operator/(complex3d const &a, complex<double> const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

inline complex3d operator/(double3d const & a, complex3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

inline complex3d operator/(complex3d const &a, double3d const & b) {
	complex3d r = a;
	r /= b;
	return r;
}

inline complex3d operator/(complex3d const &a, complex3d const &b) {
	complex3d r = a;
	r /= b;
	return r;
}

#endif
