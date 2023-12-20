///@file
///@ brief Reactive equation of state

#ifndef __REACTIVE_EOS_HPP_INCLUDED
#define __REACTIVE_EOS_HPP_INCLUDED

#include "eos/eos.hpp"
#include "simd/reactiveEOS.hpp"


/// @class IReactiveEOS
/// @brief Interface for a reactive EOS
class IReactiveEOS : public IEOS {
public:
  virtual ~IReactiveEOS() {}

  /// @brief Get EOS type
  /// @return EOS type
  virtual Type getType() const {
    return Type::REACTIVE;
  }

  /// @brief Get EOS name
  /// @return EOS name
  virtual std::string getName() const {
    return "Reactive";
  }

  /// @brief Get first eos type
  /// @return Type of first equation of state
  virtual Type getType0() const = 0;

  /// @brief Get second eos type
  /// @return Type of second equation of state
  virtual Type getType1() const = 0;

  /// @brief Initialize energy
  /// @param [in] rho_ Density
  /// @param [in] T_ Temperature
  /// @return Internal energy
  virtual double initSample(double rho_, double T_) = 0;  

};

  
/// @class ReactiveEOS
/// @brief Class for Reactive equation of state
/// @tparam EOS0 Equation of state of the reactant
/// @tparam EOS1 Equation of state of the product
template <class EOS0, class EOS1>
class ReactiveEOS : public IReactiveEOS {
public:
  static constexpr bool has_simd_opt = EOS0::has_simd_opt && EOS1::has_simd_opt; ///< Indicates that there is a simd version of that potential
  /// @brief Shorcut for the simd version of the potential
  typedef simd::kernels::Reactive<double,EOS0,EOS1> simd_opt_t;

  /// @brief Constructor
  /// @param [in] size_ Mesoparticle size
  /// @param [in] mass_ Mesoparticle mass
  /// @param [in] eos0_ Equation of state of the reactant
  /// @param [in] eos1_ Equation of state of the product
  ReactiveEOS(const uint& size_, const double& mass_, EOS0* eos0_, EOS1* eos1_) :
    size(size_), mass(mass_), eos0(eos0_), eos1(eos1_) {
  }
  
  /// @brief Destructor
  virtual ~ReactiveEOS() {}

  /// @brief Get first eos type
  /// @return Type of first equation of state
  virtual Type getType0() const {
    return eos0->getType();
  }

  /// @brief Get second eos type
  /// @return Type of second equation of state
  virtual Type getType1() const {
    return eos1->getType();
  }

  /// @brief Get mesoparticle mass
  /// @return Mass
  inline double getMass() const {
    return mass;
  }

  /// @brief Get mesoparticle size
  /// @return Mesoparticle size
  inline uint getSize() const {
    return size;
  }

  virtual double initSample(double rho_, double T_);

  /// @brief Get energy from density and entropy
  /// @param [in] rho Density
  /// @param [in] S entropy
  /// @return Energy
  inline double getEnergy(const double& rho, const double& S) {
    double rho0, S0;
    newtonTP(rho,S,rho0,S0,&EOS0::getETPAndDerivatives,&EOS1::getETPAndDerivatives);

    return (1-lambda)*eos0->getEnergy(rho0,S0)+lambda*eos1->getEnergy(rho-(1-lambda)/lambda*rho0,S-(1-lambda)/lambda*S0);
  }

  /// @brief Get temperature from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Temperature
  inline double getTemperature(const double& rho, const double& e) {
    double rho0, e0;
    newtonTP(rho,e,rho0,e0,&EOS0::getSTPAndDerivatives,&EOS1::getSTPAndDerivatives);

    return eos0->getTemperature(rho0,e0);
  }

  /// @brief Get Inverse temperature from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Inverse temperature \f[ \beta \f]
  inline double getBeta(const double& rho, const double& e) {
    double rho0, e0;
    newtonTP(rho,e,rho0,e0,&EOS0::getSTPAndDerivatives,&EOS1::getSTPAndDerivatives);

    return eos0->getBeta(rho0,e0);
  }

  /// @brief Get entropy, temperature and pressure from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Entropy, temperature and pressure
  inline vec3<double> getSTP(const double& rho, const double& e) {

    double rho0, e0;
    newtonTP(rho,e,rho0,e0,&EOS0::getSTPAndDerivatives,&EOS1::getSTPAndDerivatives);

    vec3<double> stp = eos0->getSTP(rho0,e0);
    stp[0] = (1-lambda)*stp[0]+lambda*eos1->getEntropy(rho-(1-lambda)/lambda*rho0,e-(1-lambda)/lambda*e0);
    
    return stp;

  }

  /// @brief Get energy, temperature and pressure from density and entropy
  /// @param [in] rho Density
  /// @param [in] S Entropy
  /// @return Energy, temperature and pressure
  inline vec3<double> getETP(const double& rho, const double& S) {

    double rho0, S0;
    newtonTP(rho,S,rho0,S0,&EOS0::getETPAndDerivatives,&EOS1::getETPAndDerivatives);

    vec3<double> etp = eos0->getETP(rho0,S0);
    etp[0] = (1-lambda)*etp[0]+lambda*eos1->getEnergy(rho-(1-lambda)/lambda*rho0,S-(1-lambda)/lambda*S0);

    return etp;

  }

  /// @brief Set the parameters for the simd version of the equation of state
  /// @param [in,out] opt Simd version of the equation of state
  void setParameters(simd::kernels::Reactive<double,EOS0,EOS1>& opt) {
    opt.setParameters(eos0,eos1);
  }

private:

  typedef vec3<double> (EOS0::*FThermo0)(const double&, const double&, double&, double&, double&, double&);
  typedef vec3<double> (EOS1::*FThermo1)(const double&, const double&, double&, double&, double&, double&);
  void newtonTP(const double& rho, const double& e, double& rho0, double& e0, FThermo0 fct0, FThermo1 fct1);
  
  uint size; ///< Size of the mesoparticle
  double mass; ///< Mass of the mesoparticle
  double lambda; ///< Progress variable
  
  EOS0* eos0; ///< Reactant equation of state
  EOS1* eos1; ///< Products equation of state
  
};



/// @brief Initialize internal energy
/// @param [in] rho_ Density
/// @param [in] T_ Temperature
/// @return Energy
template <class EOS0, class EOS1>
inline double ReactiveEOS<EOS0,EOS1>::initSample(double rho_, double T_) {

  return eos0->initSample(rho_,T_);

}

/// @brief Find reactant state corresponding to mixed state
/// @param [in] rho Density
/// @param [in] var Internal variable (energy or entropy)
/// @param [out] rho0 Density for the reactant
/// @param [out] var0 Internel variable (energy or entropy) for the reactant
/// @param [in] fct0 Function to minimize for the first equation of state
/// @param [in] fct1 Function to minimize for the second equation of state
template <class EOS0, class EOS1>
void ReactiveEOS<EOS0,EOS1>::newtonTP(const double& rho, const double& var, double& rho0, double& var0, FThermo0 fct0, FThermo1 fct1) {
  
  rho0 = rho;
  var0 = var;

  if (lambda == 0)
    return;

  double ratio = (1-lambda)/lambda;

  uint iter = 0;
  while( iter < 10) { // TODO: Better stopping criterion

    vec3<double> stp0, stp1; // Entropy, Temperature, Pressure
    double drt0, det0, drp0, dep0; // dT0/drho, dT0/de, dP0/drho dP0/de
    double drt1, det1, drp1, dep1; // dT1/drho, dT1/de, dP1/drho dP1/de
    
    stp0 = (eos0->*fct0)(rho0,var0,drt0,det0,drp0,dep0);
    stp1 = (eos1->*fct1)(rho-ratio*rho0,var-ratio*var0,drt1,det1,drp1,dep1);
    
    // Error function
    // X=(rho0,e0) f(X)=(T0-T1,P0-P1)
    double f1 = stp0[1]-stp1[1];
    double f2 = stp0[2]-stp1[2];
    
    // Jacobian matrix
    // J=df/dX=(d(T0-T1)/drho0, d(T0-T1)/de0; d(P0-P1)/drho0, d(P0-P1)/de0)
    double j11 = drt0+ratio*drt1;
    double j12 = det0+ratio*det1;
    double j21 = drp0+ratio*drp1;
    double j22 = dep0+ratio*dep1;
    
    double Delta = j12*j21-j11*j22;
    
    rho0 -= ( j22*f1 - j12*f2)/Delta;
    var0   -= (-j21*f1 + j11*f2)/Delta;

  }

}


#endif // __REACTIVE_EOS_HPP_INCLUDED
