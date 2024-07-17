#pragma once

#include <vector>
#include <array>
#include <cstdint>

#include <exanb/core/basic_types_def.h>
#include <exaStamp/atom_bond_connectivity.h>
#include <exaStamp/compute/force_energy.h>


namespace exaStamp
{
  using namespace exanb;

  using ChemicalBond      = std::array<uint64_t,2>;
  using ChemicalBonds     = std::vector< ChemicalBond >;
  using ChemicalAngle     = std::array<uint64_t,3>;
  using ChemicalAngles    = std::vector< ChemicalAngle >;
  using ChemicalTorsion   = std::array<uint64_t,4>;
  using ChemicalTorsions  = std::vector< ChemicalTorsion >;
  using ChemicalImproper  = std::array<uint64_t,4>;
  using ChemicalImpropers = std::vector< ChemicalImproper >;

  using MoleculeConnectivity = exaStamp::AtomBondConnectivity;
} // namespace exaStamp
