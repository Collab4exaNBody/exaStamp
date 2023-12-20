/// @file 
/// @brief Gathering of all implementations of interface IPotential

#ifndef __ALL_POTENTIALS_HPP_INCLUDED
#define __ALL_POTENTIALS_HPP_INCLUDED


#include "io/input.hpp"

#include "potential/eamVNIITFPotential.hpp"
#include "potential/meam.hpp"
#include "potential/exp6Potential.hpp"
#include "potential/idealGas.hpp"
#include "potential/lennardJonesPotential.hpp"
#include "potential/suttonChenPotential.hpp"
#include "potential/gaussianPotential.hpp"


/// @brief Temporary structure gathering all data to initialize the  potentials
template <> struct Configuration<IPotential> {

  /// @brief Default constructor
  Configuration() {}

  /// @brief Destructor (nothing to do)
  ~Configuration() {}

  Configuration(const Input& input);

  bool symmetrize; ///< Flag enabling symmetrization of the force computation

  // Ideal gas
  Array<std::string> IGtypeA; ///< First type of particle for each ideal gas potential
  Array<std::string> IGtypeB; ///< Second type of particle for each ideal gas potential

  // Lennard-Jones
  Array<std::string> LJtypeA; ///< First type of particle for each Lennard-Jones potential
  Array<std::string> LJtypeB; ///< Second type of particle for each Lennard-Jones potential

  Array<double> LJrcut;                       ///< Cutoff radius of particle for each Lennard-Jones potential
  Array<LennardJonesParameters> LJparameters; ///< Parameters for each Lennard-Jones potential

  // Exponential 6
  Array<std::string> Exp6typeA; ///< First type of particle for each Exp. 6 potential
  Array<std::string> Exp6typeB; ///< Second type of particle for each Exp. 6 potential

  Array<double> Exp6rcut;               ///< Cutoff radius for each Exp. 6 potential
  Array<Exp6Parameters> Exp6parameters; ///< Parameters for each Exp. 6 potential

  // Sutton-Chen
  Array<std::string> SCtypeA; ///< First type of particle for each Sutton-Chen potential
  Array<std::string> SCtypeB; ///< Second type of particle for each Sutton-Chen potential

  Array<double> SCrcut;                     ///< Cutoff radius for each Sutton-Chen potential
  Array<SuttonChenParameters> SCparameters; ///< Parameters for each Sutton-Chen potential

  // EAM VNIITF
  Array<std::string> EamVNIITFtypeA; ///< First type of particle for each EAM VNIITF potential
  Array<std::string> EamVNIITFtypeB; ///< Second type of particle for each EAM VNIITF potential

  Array<double> EamVNIITFrcut;                     ///< Cutoff radius for each EAM VNIITF potential
  Array<EamVniitfParameters> EamVNIITFparameters;  ///< Parameters for each EAM VNIITF potential

  // MEAM
  Array<std::string> MeamTypeA; ///< First type of particle for each MEAM potential
  Array<std::string> MeamTypeB; ///< Second type of particle for each MEAM potential

  Array<double> MeamRcut;                 ///< Cutoff radius for each MEAM potential
  Array<MeamParameters> MeamParam;			  ///< Parameters for each MEAM potential

  // Gaussian
  Array<std::string> GaussTypeA; ///< First type of particle for each Gaussian potential
  Array<std::string> GaussTypeB; ///< Second type of particle for each Gaussian potential

  Array<double> GaussRcut;                   ///< Cutoff radius of particle for each Gaussian potential
  Array<GaussianParameters> GaussParameters; ///< Parameters for each Gaussian potential


};


/// @brief Constructor from an Input structure
/// @param [in] input Input data
inline Configuration<IPotential>::Configuration(const Input& input)
  : symmetrize(input.symmetrization),
    IGtypeA(input.IGtypeA.size()), 
    IGtypeB(input.IGtypeB.size()),
    LJtypeA(input.LJtypeA.size()), 
    LJtypeB(input.LJtypeB.size()),
    LJrcut(input.LJrcut.size()), 
    LJparameters(input.LJepsilon.size()),
    Exp6typeA(input.Exp6typeA.size()), 
    Exp6typeB(input.Exp6typeB.size()),
    Exp6rcut(input.Exp6rcut.size()), 
    Exp6parameters(input.Exp6foo.size()),
    SCtypeA(input.SCtypeA.size()), 
    SCtypeB(input.SCtypeB.size()),
    SCrcut(input.SCrcut.size()), 
    SCparameters(input.SCc.size()),
    EamVNIITFtypeA(input.EamVNIITFtypeA.size()), 
    EamVNIITFtypeB(input.EamVNIITFtypeB.size()),
    EamVNIITFrcut(input.EamVNIITFrcut.size()), 
    EamVNIITFparameters(input.EamVNIITFrmax.size()),
    MeamTypeA(input.MeamTypeA.size()),
    MeamTypeB(input.MeamTypeB.size()),
    MeamRcut(input.MeamRcut.size()),
    MeamParam(input.MeamRmax.size())


 {

  // Ideal Gases
  int numIG = IGtypeA.size();
  for (int i=0; i<numIG; ++i) {
    IGtypeA[i] = input.IGtypeA[i];
    IGtypeB[i] = input.IGtypeB[i];
  }

  // Lennard-Jones potentials
  int numLJ = LJtypeA.size();
  for (int i=0; i<numLJ; ++i) {
    LJtypeA[i] = input.LJtypeA[i];
    LJtypeB[i] = input.LJtypeB[i];
    LJrcut[i]  = input.LJrcut[i];
    LJparameters[i] = LennardJonesParameters(input.LJepsilon[i], input.LJsigma[i]);
  }

  // Exp6 potentials
  int numExp6 = input.Exp6typeA.size();
  for (int i=0; i<numExp6; ++i) {
    Exp6typeA[i] = input.Exp6typeA[i];
    Exp6typeB[i] = input.Exp6typeB[i];
    Exp6rcut[i]  = input.Exp6rcut[i];
    // Warning : I'm quite confident that parameters of an exp. 6 potential should not systematically
    // be zero so this potential should not be used for now
    Exp6parameters[i] = Exp6Parameters(0.,1.,0.,0.);
  }

  // Sutton-Chen
  int numSC = SCtypeA.size();
  for (int i=0; i<numSC; ++i) {
    SCtypeA[i] = input.SCtypeA[i];
    SCtypeB[i] = input.SCtypeB[i];
    SCrcut[i]  = input.SCrcut[i];
    SCparameters[i] = SuttonChenParameters(input.SCc[i], input.SCepsilon[i], input.SCa0[i], input.SCn[i], input.SCm[i]);
  }

  // EAM VNIITF
  int numEamVNIITF = EamVNIITFtypeA.size();
  for (int i=0; i<numEamVNIITF; ++i) {
    EamVNIITFtypeA[i] = input.EamVNIITFtypeA[i];
    EamVNIITFtypeB[i] = input.EamVNIITFtypeB[i];
    EamVNIITFrcut[i]  = input.EamVNIITFrcut[i];
    EamVNIITFparameters[i] = EamVniitfParameters(input.EamVNIITFrmax[i], input.EamVNIITFrmin[i], input.EamVNIITFrt0[i], 
						 input.EamVNIITFEcoh[i], input.EamVNIITFE0[i], input.EamVNIITFbeta[i], 
						 input.EamVNIITFA[i], input.EamVNIITFZ[i], input.EamVNIITFn[i], 
						 input.EamVNIITFalpha[i], input.EamVNIITFD[i], input.EamVNIITFeta[i], 
						 input.EamVNIITFmu[i]);
  }
 // MEAM
  int numMeam = MeamTypeA.size();
  for (int i=0; i<numMeam; ++i) {
    MeamTypeA[i] = input.MeamTypeA[i];
    MeamTypeB[i] = input.MeamTypeB[i];
    MeamRcut[i]  = input.MeamRcut[i];
    MeamParam[i] = MeamParameters(input.MeamRmax[i], input.MeamRmin[i],
						 input.MeamEcoh[i] ,input.MeamE0[i]   , input.MeamA[i]   ,input.MeamR0[i], input.MeamAlpha[i], input.MeamDelta[i],
						 input.MeamBeta0[i],input.MeamBeta1[i],input.MeamBeta2[i],
						 input.MeamBeta3[i],input.MeamT0[i] ,input.MeamT1[i]   ,input.MeamT2[i]   ,input.MeamT3[i],
						 input.MeamS0[i],input.MeamS1[i]   ,input.MeamS2[i]   ,input.MeamS3[i],
						 input.MeamCmin[i], input.MeamCmax[i],
						 input.MeamZ[i],input.MeamRc[i],input.MeamRp[i]);
  	}

  // Gaussian potentials
  int numGauss = GaussTypeA.size();
  for (int i=0; i<numGauss; ++i) {
    GaussTypeA[i] = input.GaussTypeA[i];
    GaussTypeB[i] = input.GaussTypeB[i];
    GaussRcut[i]  = input.GaussRcut[i];
    GaussParameters[i] = GaussianParameters(input.GaussRatt[i], input.GaussRrep[i], input.GaussRatio[i], input.GaussEpsilon[i]);
  }

}

#endif // __ALL_POTENTIALS_HPP_INCLUDED
