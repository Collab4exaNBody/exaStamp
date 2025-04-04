#include <exaStamp/ttm/source_term.h>
#include <onika/log.h>
#include <onika/math/basic_types_operators.h>

#include <iostream>
#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  // -----------------------------------------------
  // ------- Source term factory -------------------
  // -----------------------------------------------
  /*
   * uses 2D gaussian function as described below :
   * https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
   * 
   * X = distance from center 'c' , Y = time
   */
  class SphericalTemporalSourceTerm: public ScalarSourceTerm
  {
  public:
    inline SphericalTemporalSourceTerm(const Vec3d& c, double amplitude, double radius_mean, double radius_dev, double time_mean, double time_dev)
      : m_center(c)
      , m_amplitude(amplitude)
      , m_x0(radius_mean)
      , m_2_sigma_x_sqr( 2.0 * radius_dev * radius_dev )
      , m_y0(time_mean)
      , m_2_sigma_y_sqr( 2.0 * time_dev * time_dev )
      {}

    virtual inline double S( Vec3d r, double t ) const
    {
      double x = norm(r-m_center) - m_x0;
      double y = t - m_y0;
      return m_amplitude * std::exp( - ( (x*x)/m_2_sigma_x_sqr + (y*y)/m_2_sigma_y_sqr ) );
    }
    
  private:
    Vec3d m_center;
    double m_amplitude;
    double m_x0;
    double m_2_sigma_x_sqr;
    double m_y0;
    double m_2_sigma_y_sqr;
  };

  struct ConstantSourceTerm : public ScalarSourceTerm
  {
    inline ConstantSourceTerm(double s) : m_scalar(s) {}
    virtual inline double S( Vec3d r, double t ) const { return m_scalar; }
  private:
    double m_scalar = 0.0;
  };

  std::shared_ptr<ScalarSourceTerm> make_source_term( const std::string& stype, const std::vector<double>& p )
  {
//    ldbg << "source type = "<< stype << std::endl;
    if( stype == "null" )
    {
      return std::make_shared<ScalarSourceTerm>();
    }
    else if( stype == "sphere" )
    {
      if( p.size() != 8 )
      {
        lerr << "expected 8 parameters for sphere source type"<<std::endl;
        std::abort();
      }
      return std::make_shared<SphericalTemporalSourceTerm>( Vec3d{p[0],p[1],p[2]} , p[3] , p[4] , p[5] , p[6] , p[7] );
    }
    else if( stype == "constant" )
    {
      if( p.size() != 1 )
      {
        lerr << "expected 1 parameters for sphere source type"<<std::endl;
        std::abort();
      }
      return std::make_shared<ConstantSourceTerm>( p[0] );
    }
    else
    {
      lerr << "unrecognized source type '"<<stype<<"'"<<std::endl;
      std::abort();
    }
    return nullptr;
  }

}

