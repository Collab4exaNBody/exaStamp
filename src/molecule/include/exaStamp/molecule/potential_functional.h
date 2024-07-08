#pragma once

//#include <utility>
#include <exaStamp/molecule/molecule_compute_param.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  /**
   * Signum function : give the sign of an expression (+1, -1 or 0)
   */
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value, int>::type
  ONIKA_HOST_DEVICE_FUNC
  inline constexpr signum(T x) noexcept {
    return T(0) < x;
  }

  template <typename T>
  typename std::enable_if<std::is_signed<T>::value, int>::type
  ONIKA_HOST_DEVICE_FUNC  
  inline constexpr signum(T x) noexcept {
    return (T(0) < x) - (x < T(0));
  }

  struct ScalarForceEnergy
  { // names are first and second for bacward compatibility with std::pair used until now
    double first; //force; 
    double second; //energy;
  };

  ONIKA_HOST_DEVICE_FUNC
  inline ScalarForceEnergy intramolecular_quar( const MoleculeGenericFuncParam& m , double t)
  {
    const auto k2 = m.p0;
    const auto k3 = m.p1;
    const auto k4 = m.p2;
    const auto t0 = m.x0;
    // const auto energy_scale = m.coeff;
    const double tmp = t - t0;
    return { 2 * k2 *      tmp
           + 3 * k3 *  pow(tmp,2)
           + 4 * k4 *  pow(tmp,3)
           ,
           ( k2 * pow(tmp,2)
           + k3 * pow(tmp,3)
           + k4 * pow(tmp,4) ) //* energy_scale
           };  
  }

  ONIKA_HOST_DEVICE_FUNC
  inline ScalarForceEnergy intramolecular_cos( const MoleculeGenericFuncParam& m , double _phi )
  {
    const auto k1 = m.p0;
    const auto k2 = m.p1;
    const auto k3 = m.p2;
    const auto phi0 = m.x0;
    const auto opls_fac = m.coeff;

    const double phi = _phi - phi0;
    const double sin_phi = sin(phi);
    const double sin_2phi = sin(2*phi);
    const double sin_3phi = sin(3*phi);
    const double cos_phi = cos(phi);
    const double cos_2phi = cos(2*phi);
    const double cos_3phi = cos(3*phi);
    
    return {
       opls_fac *     k1 * sin_phi
     +            2 * k2 * sin_2phi
     + opls_fac * 3 * k3 * sin_3phi
     ,
       k1 * ( 1 - opls_fac * cos_phi )
     + k2 * ( 1 -            cos_2phi )
     + k3 * ( 1 - opls_fac * cos_3phi )
     };
  }


  // base definition for an intra molecular ptoential. might also be used as a null potential
  class IntraMolecularPotentialFunctional
  {
  public:
    virtual inline ScalarForceEnergy force_energy(double) const
    {
      return { 0., 0. };
    }
    virtual inline MoleculeGenericFuncParam generic_parameters() const
    {
      return {0.,0.,0.,0.};
    }
    virtual ~IntraMolecularPotentialFunctional() = default;
  };

  // Harm potential
  class IntraMolecularHarmFunctional : public IntraMolecularPotentialFunctional
  {
  public:
    inline IntraMolecularHarmFunctional(double _k, double _t0) : k(_k), t0(_t0) {}
    inline ScalarForceEnergy force_energy(double t) const override final
    {
      return intramolecular_quar( generic_parameters() , t );
/*      return {       k *     (t - t0)
             , 0.5 * k * pow((t - t0),2) }; */
    }
    inline MoleculeGenericFuncParam generic_parameters() const override final
    {
      return {k,0.,0.,t0,1.0f};
    }
  private:
    double k = 0.0;
    float t0 = 0.0;
  };

  // OPLS Bond potential
  class IntraMolecularBondOPLSFunctional : public IntraMolecularPotentialFunctional
  {
  public:
    inline IntraMolecularBondOPLSFunctional(double _k, double _t0) : k(_k), t0(_t0) {}
    inline ScalarForceEnergy force_energy(double t) const override final
    {
      return intramolecular_quar( generic_parameters() , t );
/*      return { 2.0 * k *     (t - t0)
             ,       k * pow((t - t0),2) }; */
    }
    inline MoleculeGenericFuncParam generic_parameters() const override final
    {
      return {k,0.,0.,t0,1.0f};
    }
  private:
    double k = 0.0;
    float t0 = 0.0;
  };

  // Quar potential
  class IntraMolecularQuarFunctional : public IntraMolecularPotentialFunctional
  {
  public:
    inline IntraMolecularQuarFunctional(double _k2, double _k3, double _k4, double _t0) : k2(_k2), k3(_k3), k4(_k4), t0(_t0) {}
    inline ScalarForceEnergy force_energy(double t) const override final
    {
      return intramolecular_quar( generic_parameters() , t );
/*      double tmp = t-t0;
      return { 2 * k2 *      tmp
             + 3 * k3 *  pow(tmp,2)
             + 4 * k4 *  pow(tmp,3)
             ,
               k2 * pow(tmp,2)
             + k3 * pow(tmp,3)
             + k4 * pow(tmp,4)
             };*/
    }
    inline MoleculeGenericFuncParam generic_parameters() const override final
    {
      return {k2,k3,k4,t0,1.0f};
    }
  private:
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;
    float t0 = 0.0;
  };

  // Compass potential
  class IntraMolecularCompassFunctional : public IntraMolecularPotentialFunctional
  {
  public:
    inline IntraMolecularCompassFunctional(double _k1, double _k2, double _k3) : k1(_k1), k2(_k2), k3(_k3) {}
    inline ScalarForceEnergy force_energy(double phi) const override final
    {
      return intramolecular_cos( generic_parameters() , phi );
/*      return {
             k1 * sin(  phi)
       + 2 * k2 * sin(2*phi)
       + 3 * k3 * sin(3*phi)
       ,
         k1 * ( 1 - cos(  phi) )
       + k2 * ( 1 - cos(2*phi) )
       + k3 * ( 1 - cos(3*phi) )
       }; */
    }
    inline MoleculeGenericFuncParam generic_parameters() const override final
    {
      return {k1,k2,k3,0.0f,1.0f};
    }
  private:
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
  };

  // Half Compass potential
  class IntraMolecularHalfCompassFunctional : public IntraMolecularPotentialFunctional
  {
  public:
    inline IntraMolecularHalfCompassFunctional(double _k1, double _k2, double _k3) : k1(_k1), k2(_k2), k3(_k3) {}
    inline ScalarForceEnergy force_energy(double phi) const override final
    {
      return intramolecular_cos( generic_parameters() , phi );
/*      return {
           0.5 * k1 * sin(  phi)
         +       k2 * sin(2*phi)
         + 1.5 * k3 * sin(3*phi)
       ,
           k1/2 * ( 1 - cos(  phi) )
         + k2/2 * ( 1 - cos(2*phi) )
         + k3/2 * ( 1 - cos(3*phi) )
       }; */
    }
    inline MoleculeGenericFuncParam generic_parameters() const override final
    {
      return {k1/2,k2/2,k3/2,0.0f,1.0f};
    }
  private:
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
  };

  // OPLS potential
  class IntraMolecularCosOPLSFunctional : public IntraMolecularPotentialFunctional
  {
  public:
    inline IntraMolecularCosOPLSFunctional(double _k1, double _k2, double _k3) : k1(_k1), k2(_k2), k3(_k3) {}
    inline ScalarForceEnergy force_energy(double phi) const override final
    {
      return intramolecular_cos( generic_parameters() , phi );
/*      return {
         - 0.5 * k1 * sin(  phi)
         +       k2 * sin(2*phi)
         - 1.5 * k3 * sin(3*phi)
       ,
           0.5 * k1 * ( 1 + cos(  phi) )
         + 0.5 * k2 * ( 1 - cos(2*phi) )
         + 0.5 * k3 * ( 1 + cos(3*phi) )
       }; */
    }
    inline MoleculeGenericFuncParam generic_parameters() const override final
    {
      return {k1/2,k2/2,k3/2,0.0f,-1.0f};
    }
  private:
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
  };

  // CosTwo potential
  class IntraMolecularCosTwoFunctional : public IntraMolecularPotentialFunctional
  {
  public:
    inline IntraMolecularCosTwoFunctional(double _k, double _phi0) : k(_k), phi0(_phi0) {}
    inline ScalarForceEnergy force_energy(double phi) const override final
    {
      return intramolecular_cos( generic_parameters() , phi );
/*
      double tmp = phi - phi0;
      return {       k * sin(2*tmp)
             , 0.5 * k * ( 1 - cos(2*tmp) ) 
             };
*/
    }
    inline MoleculeGenericFuncParam generic_parameters() const override final
    {
      return {0.,k/2,0.,phi0,1.0f};
    }
  private:
    double k = 0.0;
    float phi0 = 0.0;
  };


}

