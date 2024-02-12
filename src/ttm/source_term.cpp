#include <exaStamp/ttm/source_term.h>
#include <exanb/core/log.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/core/quantity_yaml.h>
#include <exanb/core/basic_types_yaml.h>

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

    virtual inline double operator () ( Vec3d r, double t=0.0, int64_t id=-1 ) const override final
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

  struct WaveFrontSourceTerm : public ScalarSourceTerm
  {
    inline WaveFrontSourceTerm(const Plane3d& ref_plane, const Plane3d& wave_plane, double amplitude)
      : m_ref_plane(ref_plane)
      , m_wave_plane(wave_plane)
      , m_amplitude(amplitude)
    {}
    virtual inline double operator () ( Vec3d r, double t=0.0, int64_t id=-1 ) const override final
    {
      const double p = dot( r , m_ref_plane.N ) + m_ref_plane.D;
      const double w = dot( r , m_wave_plane.N ) + m_wave_plane.D;
      return p + std::sin( w ) * m_amplitude;
    }
  private:
    Plane3d m_ref_plane = { 1.0 , 0.0 , 0.0 , 0.0 };
    Plane3d m_wave_plane = { 0.0 , 1.0 , 0.0 , 0.0 };
    double m_amplitude = 1.0;
  };


  struct ConstantSourceTerm : public ScalarSourceTerm
  {
    inline ConstantSourceTerm(double s) : m_scalar(s) {}
    virtual inline double operator () ( Vec3d r, double t=0.0, int64_t id=-1 ) const override final
    {
      return m_scalar;
    }
  private:
    double m_scalar = 0.0;
  };


  ScalarSourceTermInstance make_source_term( const YAML::Node& node )
  {
    if( node.IsScalar() )
    {
      std::string type = node.as< std::string >();
      if( type == "null" ) return std::make_shared<ScalarSourceTerm>();
      else return nullptr;
    }

    if( ! node.IsMap() ) { return nullptr; }
    if( node.size() != 1 ) { return nullptr; }
          
    std::string type = node.begin()->first.as<std::string>();
    YAML::Node params = node.begin()->second;
    
    if( type == "sphere" )
    {
      return std::make_shared<SphericalTemporalSourceTerm>( params["center"].as<Vec3d>()
                                                          , params["amplitude"].as<Quantity>().convert()
                                                          , params["radius_mean"].as<Quantity>().convert()
                                                          , params["radius_dev"].as<Quantity>().convert()
                                                          , params["time_mean"].as<Quantity>().convert()
                                                          , params["time_dev"].as<Quantity>().convert()
                                                          );
    }
    else if( type == "wavefront" )
    {
      return std::make_shared<WaveFrontSourceTerm>( params["plane"].as<Plane3d>(), params["wave"].as<Plane3d>(), params["amplitude"].as<Quantity>().convert() );
    }
    else if( type == "constant" )
    {
      return std::make_shared<ConstantSourceTerm>( params.as<Quantity>().convert() );
    }    
    else
    {
      lerr << "unrecognized source type '"<<type<<"'"<<std::endl;
      std::abort();
    }
    return nullptr;
  }

}

