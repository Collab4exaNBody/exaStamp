#include <yaml-cpp/yaml.h>

#include <string>
#include <cstdlib>

struct ParticleReferenceValue
{
  uint64_t m_id = -1;
  double m_r[3] = { 0., 0., 0. };
  double m_a[3] = { 0., 0., 0. };
  double m_v[3] = { 0., 0., 0. };
  inline bool operator < (const ParticleReferenceValue& x) const { return m_id < x.m_id; }
};

namespace YAML
{
  template<> struct convert< ParticleReferenceValue >
  {
    static inline double as_hexfloat_double(const Node& node)
    {
      std::string s = node.as<std::string>();
      return std::strtod( s.c_str() , NULL ); // reading from std::istream seems to be buggy, using old style strtod reads correctly hexfloats
  /*      std::istringstream iss(s);
      double d = std::numeric_limits<double>::quiet_NaN();
      iss >> std::hexfloat >> d;
      return d;*/
    }
    static inline bool decode(const Node& node, ParticleReferenceValue& v)
    {
      if( ! node.IsSequence() ) { return false; }
      if( node.size() < 7 ) { return false; }
      v.m_id = node[0].as<uint64_t>();
      v.m_r[0] = as_hexfloat_double(node[1]);
      v.m_r[1] = as_hexfloat_double(node[2]);
      v.m_r[2] = as_hexfloat_double(node[3]);
      v.m_a[0] = as_hexfloat_double(node[4]);
      v.m_a[1] = as_hexfloat_double(node[5]);
      v.m_a[2] = as_hexfloat_double(node[6]);
      if( node.size() >= 9 )
      {
      v.m_v[0] = as_hexfloat_double(node[7]);
      v.m_v[1] = as_hexfloat_double(node[8]);
      v.m_v[2] = as_hexfloat_double(node[9]);
      }
      return true;
    }
  };
}

static inline std::ostream& print_particle_reference_value(std::ostream& out, const ParticleReferenceValue& p, std::ios_base&(*manip)(std::ios_base&) = std::defaultfloat )
{
  out << "[" << p.m_id
  << ", " << manip << p.m_r[0] 
  << ", " << manip << p.m_r[1] 
  << ", " << manip << p.m_r[2] 
  << ", " << manip << p.m_a[0] 
  << ", " << manip << p.m_a[1] 
  << ", " << manip << p.m_a[2]
  << ", " << manip << p.m_v[0] 
  << ", " << manip << p.m_v[1] 
  << ", " << manip << p.m_v[2]
  <<"]";
  return out;
}

static inline std::ostream& operator << ( std::ostream& out, const ParticleReferenceValue& p )
{
  return print_particle_reference_value(out,p,std::hexfloat);
}

struct ParticleError
{
  ParticleReferenceValue value;
  int rank=-1;
  double re=0.0;
  double ae=0.0;
  double ve=0.0;
};


