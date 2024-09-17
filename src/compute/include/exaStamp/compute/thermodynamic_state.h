#pragma once

#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <cstdlib>

#include <exanb/core/physics_constants.h> //NICO

namespace exaStamp
{
  using namespace exanb;

  class ThermodynamicState
  {
  public:
    inline void set_virial(const Mat3d& x) { m_virial = x; }
    inline void set_ke_tensor(const Mat3d& x) { m_ke_tensor = x; }
    inline const Mat3d& virial() const { return m_virial; }
    inline const Mat3d& ke_tensor() const { return m_ke_tensor; }    

    inline void set_pressure(const Vec3d& x) { m_pressure = x; }
    inline void set_pressurebis(const Vec3d& x) { m_pressurebis = x; }    
    inline void set_deviator(const Vec3d& x) { m_deviator = x; }
    inline const Vec3d& pressure() const { return m_pressure; }
    inline const Vec3d& pressurebis() const { return m_pressurebis; }    
    inline const Vec3d& deviator() const { return m_deviator; }
    inline double pressure_scal() const
    {
      double kinetic_pressure = temperature_scal() / volume();
      Vec3d virdiag = { m_virial.m11 , m_virial.m22, m_virial.m33 };
      double potential_pressure = ( virdiag.x + virdiag.y + virdiag.z ) / ( 3. * volume() );
//double conv_P=1.e4 * legacy_constant::atomicMass * 1e30; //NICO
//double conv_T=1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann / m_particle_count; //NICO
//printf("PRESS: %lf -> %e + %e = %e \n",temperature_scal()*conv_T,kinetic_pressure*conv_P,potential_pressure*conv_P,(kinetic_pressure+potential_pressure)*conv_P); //NICO
      return kinetic_pressure + potential_pressure;
    }

    inline void set_vonmises(const Vec3d& x) { m_vonmises = x; }
    inline const Vec3d& vonmises() const { return m_vonmises; }
    inline double vonmises_scal() const
    {
      Vec3d virdiag = { m_virial.m11 , m_virial.m22, m_virial.m33 };
      Vec3d virdeviatoric = { m_virial.m12 , m_virial.m13, m_virial.m23 };
      double vonmises = sqrt( 0.5 * (
                                     ( virdiag.x - virdiag.y ) * ( virdiag.x - virdiag.y ) 
                                   + ( virdiag.y - virdiag.z ) * ( virdiag.y - virdiag.z ) 
                                   + ( virdiag.z - virdiag.x ) * ( virdiag.z - virdiag.x ) 
                                   + 6.0 * ( virdeviatoric.x * virdeviatoric.x + virdeviatoric.y * virdeviatoric.y + virdeviatoric.z * virdeviatoric.z)
                                   ) ) 
                        / volume();
      return vonmises;
    }

    inline Mat3d stress_tensor() const
    {
      Mat3d S = (1./volume()) * m_virial;
      return S;
    }
    
    inline void set_kinetic_energy(const Vec3d& x) { m_kinetic_energy = x; }
    inline const Vec3d& kinetic_energy() const { return m_kinetic_energy; }
    inline double kinetic_energy_scal() const { return m_kinetic_energy.x + m_kinetic_energy.y + m_kinetic_energy.z; }

    inline void set_rotational_energy(const Vec3d& x) { m_rotational_energy = x; }
    inline const Vec3d& rotational_energy() const { return m_rotational_energy; }
    inline double rotational_energy_scal() const { return m_rotational_energy.x + m_rotational_energy.y + m_rotational_energy.z; }

    inline void set_ndof(const Vec3d& x) { m_ndof = x; }
    inline const Vec3d& ndof() const { return m_ndof; }
    inline double ndof_scal() const { return m_ndof.x + m_ndof.y + m_ndof.z; }
    
    inline void set_temperature(const Vec3d& x) { m_temperature = x; }
    inline const Vec3d& temperature() const { return m_temperature; }
    inline double temperature_scal() const { return ( m_temperature.x + m_temperature.y + m_temperature.z ) / 3. ; }
    inline double temperature_rigidmol_scal() const { return ( m_temperature.x + m_temperature.y + m_temperature.z ) / (3. + ndof_scal() / particle_count()) ; }

    inline void set_kinetic_temperature(const Vec3d& x) { m_kinetic_temperature = x; }
    inline const Vec3d& kinetic_temperature() const { return m_kinetic_temperature; }
    inline double kinetic_temperature_x() const { return m_kinetic_temperature.x ; }
    inline double kinetic_temperature_y() const { return m_kinetic_temperature.y ; }
    inline double kinetic_temperature_z() const { return m_kinetic_temperature.z ; }
    inline double kinetic_temperature_scal() const { return ( m_kinetic_temperature.x + m_kinetic_temperature.y + m_kinetic_temperature.z ) / (3. + ndof_scal() / particle_count()) ; }

    inline void set_rotational_temperature(const Vec3d& x) { m_rotational_temperature = x; }
    inline const Vec3d& rotational_temperature() const { return m_rotational_temperature; }
    inline double rotational_temperature_x() const { return m_rotational_temperature.x ; }
    inline double rotational_temperature_y() const { return m_rotational_temperature.y ; }
    inline double rotational_temperature_z() const { return m_rotational_temperature.z ; }
    inline double rotational_temperature_scal() const { return ( m_rotational_temperature.x + m_rotational_temperature.y + m_rotational_temperature.z ) / (3. + ndof_scal() / particle_count()) ; }

    inline void set_kinetic_momentum(const Vec3d& x) { m_kinetic_momentum = x; }
    inline const Vec3d& kinetic_momentum() const { return m_kinetic_momentum; }
    
    inline void set_potential_energy(double x) { m_potential_energy = x; }
    inline double potential_energy() const { return m_potential_energy; }

    inline void set_internal_energy(double x) { m_internal_energy = x; }
    inline double internal_energy() const { return m_internal_energy; }

    inline void set_chemical_energy(double x) { m_chemical_energy = x; }
    inline double chemical_energy() const { return m_chemical_energy; }
    
    inline double total_energy() const { return kinetic_energy_scal() + potential_energy() + internal_energy() + chemical_energy(); }
    inline double total_energy_rigidmol() const { return kinetic_energy_scal() + rotational_energy_scal() + potential_energy() + internal_energy() + chemical_energy(); }
    
    inline void set_mass(double x) { m_mass = x; }
    inline double mass() const { return m_mass; }

    inline void set_volume(double x) { m_volume = x; }
    inline double volume() const { return m_volume; }

    inline void set_particle_count(size_t x) { m_particle_count = x; }
    inline size_t particle_count() const { return m_particle_count; }
    
  private:
    Mat3d m_virial;
    Mat3d m_ke_tensor;
    Vec3d m_pressure;
    Vec3d m_pressurebis;    
    Vec3d m_deviator;
    Vec3d m_vonmises;
    Vec3d m_kinetic_energy;
    Vec3d m_rotational_energy;
    Vec3d m_temperature;
    Vec3d m_kinetic_temperature;
    Vec3d m_rotational_temperature;
    Vec3d m_kinetic_momentum;
    Vec3d m_ndof;
    double m_potential_energy = 0.;
    double m_internal_energy = 0.;
    double m_chemical_energy = 0.;
    double m_mass = 0.;
    double m_volume = 0.;
    size_t m_particle_count = 0;
  };
  
}
