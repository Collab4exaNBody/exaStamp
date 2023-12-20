/// @file
/// @brief Gathering of all implementations of interface IEos

#ifndef __ALLEOS_HPP_INCLUDED
#define __ALLEOS_HPP_INCLUDED


#include "io/input.hpp"
#include "eos/eos.hpp"
#include "eos/dpdeEOS.hpp"
#include "eos/idealGasEOS.hpp"
#include "eos/mieGruneisenEOS.hpp"
#include "eos/hzEOS.hpp"
#include "eos/jwlEOS.hpp"
#include "eos/reactiveEOS.hpp"

/// @brief Structure gathering all data to initialize EOS
template<> struct Configuration<IEOS>  {
  /// @brief Default Constructor
  Configuration() {}
  
  /// @brief Destructor
  ~Configuration() {}
  
  /// @brief Constructor from an Input structure
  Configuration(Input & input);

  Array<std::string> DPDEtype; ///< eos for DPDE
  Array<double> DPDEcv; ///< cv for DPDE

  Array<std::string> IGtype; ///< ideal gas eos

  Array<std::string> MGtype; ///< Mie-Gruneisen eos
  Array<double> MGgamma0; ///< \f[ \Gamma_0 \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGgammaInf; ///< \f[ \Gamma_{\infty} \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGtheta0; ///< \f[ \theta_0 \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGq; ///< \f[ q \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGrho0; ///< \f[ \rho_0 \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGks; ///< \f[ K_S \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGns; ///< \f[ N_S \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGrhos; ///< \f[ \rho_S \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGur; ///< \f[ u_r \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGcvr; ///< \f[ Cv_r \f] parameter for the Mie-Gruneisen equation of state

  Array<std::string> HZtype; ///< HZ eos
  Array<double> HZgamma0; ///< Gruneisen parameter for the HZ equation of state
  Array<double> HZrho0; ///< Reference density for the HZ equation of state
  Array<double> HZc0; ///< \f[ c_0 \f] parameter for the HZ equation of state
  Array<double> HZcv; ///< Heat capacity for the HZ equation of state
  Array<double> HZs; ///< \f[ s \f] parameter for the HZ equation of state

  Array<std::string> JWLtype; ///< Jones-Wilkins-Lee eos
  Array<double> JWLgamma0; ///< Gruneisen parameter for the JWL equation of state
  Array<double> JWLrho0; ///< Reference density for the JWL equation of state
  Array<double> JWLe0; ///< Reference energy for the JWL equation of state
  Array<double> JWLdcj; ///< Detonation velocity
  Array<double> JWLpcj; ///< Pressure at the CJ point
  Array<double> JWLtcj; ///< Temperature at the CJ point
  Array<double> JWLcv; ///< Heat capacity for the JWL equation of state
  Array<double> JWLa; ///< \f[ a \f] parameter for the JWL equation of state
  Array<double> JWLb; ///< \f[ b \f] parameter for the JWL equation of state
  Array<double> JWLr1; ///< \f[ R_1 \f] parameter for the JWL equation of state
  Array<double> JWLr2; ///< \f[ R_2 \f] parameter for the JWL equation of state

  Array<std::string> Reactivetype; ///< Reactive eos
  Array<std::string> ReactiveEOS0; ///< First eos
  Array<std::string> ReactiveEOS1; ///< Second eos
  
};


// Constructor from an Input structure
inline Configuration<IEOS>::Configuration(Input& input)
  : DPDEtype(input.DPDEtype.size()),
    DPDEcv(input.DPDEcv.size()),
    IGtype(input.IGtype.size()),
    MGtype(input.MGtype.size()),
    MGgamma0(input.MGgamma0.size()),
    MGgammaInf(input.MGgammaInf.size()),
    MGtheta0(input.MGtheta0.size()),
    MGq(input.MGq.size()),
    MGrho0(input.MGrho0.size()),
    MGks(input.MGks.size()),
    MGns(input.MGns.size()),
    MGrhos(input.MGrhos.size()),
    MGur(input.MGur.size()),
    MGcvr(input.MGcvr.size()),
    HZtype(input.HZtype.size()),
    HZgamma0(input.HZgamma0.size()),
    HZrho0(input.HZrho0.size()),
    HZc0(input.HZc0.size()),
    HZcv(input.HZcv.size()),
    HZs(input.HZs.size()),
    JWLtype(input.JWLtype.size()),
    JWLgamma0(input.JWLgamma0.size()),
    JWLrho0(input.JWLrho0.size()),
    JWLe0(input.JWLe0.size()),
    JWLdcj(input.JWLdcj.size()),
    JWLpcj(input.JWLpcj.size()),
    JWLtcj(input.JWLtcj.size()),
    JWLcv(input.JWLcv.size()),
    JWLa(input.JWLa.size()),
    JWLb(input.JWLb.size()),
    JWLr1(input.JWLr1.size()),
    JWLr2(input.JWLr2.size()),
    Reactivetype(input.Reactivetype.size()),
    ReactiveEOS0(input.ReactiveEOS0.size()),
    ReactiveEOS1(input.ReactiveEOS1.size()) {

  // DPDE
  int numDPDE = DPDEtype.size();
  for (int i=0; i<numDPDE; ++i) {
    DPDEtype[i] = input.DPDEtype[i];
    DPDEcv[i] = input.DPDEcv[i];
  }

  // Ideal Gases
  int numIG = IGtype.size();
  for (int i=0; i<numIG; ++i)
    IGtype[i] = input.IGtype[i];

  // Mie-Grunesein EOS
  int numMG = MGtype.size();
  for (int i=0; i<numMG; ++i) {
    MGtype[i] = input.MGtype[i];
    MGgamma0[i] = input.MGgamma0[i];
    MGgammaInf[i] = input.MGgammaInf[i];
    MGtheta0[i] = input.MGtheta0[i];
    MGq[i] = input.MGq[i];
    MGrho0[i] = input.MGrho0[i];
    MGks[i] = input.MGks[i];
    MGns[i] = input.MGns[i];
    MGrhos[i] = input.MGrhos[i];
    MGur[i] = input.MGur[i];
    MGcvr[i] = input.MGcvr[i];
  }

  // HZ EOS
  int numHZ = HZtype.size();
  for (int i=0; i<numHZ; ++i) {
    HZtype[i] = input.HZtype[i];
    HZgamma0[i] = input.HZgamma0[i];
    HZrho0[i] = input.HZrho0[i];
    HZc0[i] = input.HZc0[i];
    HZcv[i] = input.HZcv[i];
    HZs[i] = input.HZs[i];
  }

  // JWL EOS
  int numJWL = JWLtype.size();
  for (int i=0; i<numJWL; ++i) {
    JWLtype[i] = input.JWLtype[i];
    JWLgamma0[i] = input.JWLgamma0[i];
    JWLrho0[i] = input.JWLrho0[i];
    JWLe0[i] = input.JWLe0[i];
    JWLdcj[i] = input.JWLdcj[i];
    JWLpcj[i] = input.JWLpcj[i];
    JWLtcj[i] = input.JWLtcj[i];
    JWLcv[i] = input.JWLcv[i];
    JWLa[i] = input.JWLa[i];
    JWLb[i] = input.JWLb[i];
    JWLr1[i] = input.JWLr1[i];
    JWLr2[i] = input.JWLr2[i];
  }

  // Reactive EOS
  int numReactive = Reactivetype.size();
  for (int i=0; i<numReactive; ++i) {
    Reactivetype[i] = input.Reactivetype[i];
    ReactiveEOS0[i] = input.ReactiveEOS0[i];
    ReactiveEOS1[i] = input.ReactiveEOS1[i];
  }
  
}

#endif // __ALLEOS_HPP_INCLUDED
