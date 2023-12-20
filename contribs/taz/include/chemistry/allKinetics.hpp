/// @file
/// @brief Gathering of all implementations of interface IKinetics

#ifndef __ALLKINETICS_HPP_INCLUDED
#define __ALLKINETICS_HPP_INCLUDED


#include "io/input.hpp"
#include "chemistry/kinetics.hpp"
#include "chemistry/secondOrder.hpp"

/// @brief Structure gathering all data to initialize chemical reactions
template<> struct Configuration<IKinetics>  {
  /// @brief Default Constructor
  Configuration() {}
  
  /// @brief Destructor
  ~Configuration() {}
  
  /// @brief Constructor from an Input structure
  Configuration(Input & input);

  Array<std::string> SOtype; ///< Particle types reacting with a second order kinetics
  Array<double> SOzab; ///< Arrhenius prefactor for direct second order reaction (A->B)
  Array<double> SOzba; ///< Arrhenius prefactor for reverse second order reaction (A<-B)
  Array<double> SOeab; ///< Activation energy for direct second order reaction (A->B)
  Array<double> SOeba; ///< Activation energy for reverse second order reaction (A<-B)
    
};


// Constructor from an Input structure
inline Configuration<IKinetics>::Configuration(Input& input)
  : SOtype(input.SOtype.size()),
    SOzab(input.SOzab.size()),
    SOzba(input.SOzba.size()),
    SOeab(input.SOeab.size()),
    SOeba(input.SOeba.size()) {

  // Second order
  int numSO = SOtype.size();
  for (int i=0; i < numSO; ++i) {
    SOtype[i] = input.SOtype[i];
    SOzab[i] = input.SOzab[i];
    SOzba[i] = input.SOzba[i];
    SOeab[i] = input.SOeab[i];
    SOeba[i] = input.SOeba[i];
  }
  
}

#endif // __ALLKINETICS_HPP_INCLUDED
