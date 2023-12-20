/*
 * forceField.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: giarda
 */

/// @file
/// @brief Forcefield namespace : implementations


#include "forceField/forceField.hpp"


double ForceField::s_rcut=1.05;
const std::map<FF::Type,std::tuple<uint,uint,uint,uint>> ForceField::paramBuild={{FF::UFF,std::make_tuple(2,4,3,3)},{FF::COMPASS,std::make_tuple(0,0,0,0)},{FF::AMBER,std::make_tuple(0,0,0,0)},{FF::EXPERT,std::make_tuple(2,4,3,3)}};
const std::map<FF::FunctForm,int> ForceField::numParam={{FF::BHARM,2},{FF::ACOSTO3,4},{FF::DCOSN,3},{FF::ICOSTO2,3}};
const std::map<FF::Type,FF::ChargeMethod> ForceField::defaultCharge={{FF::UFF,FF::NONE},{FF::COMPASS,FF::FROMCOMPASS},{FF::AMBER,FF::OPENBABEL},{FF::EXPERT,FF::NONE}};
