#pragma once

#include <vector>
#include <array>
#include <cstdint>

#include <exanb/core/basic_types_def.h>
#include <exaStamp/atom_bond_connectivity.h>

// #define XSTAMP_USE_BEND_FORCE_BUFFER 1

namespace exaStamp
{
  using namespace exanb;

  using ChemicalBonds     = std::vector< std::array<uint64_t,2> >;
  using ChemicalAngles    = std::vector< std::array<uint64_t,3> >;
  using ChemicalTorsions  = std::vector< std::array<uint64_t,4> >;
  using ChemicalImpropers = std::vector< std::array<uint64_t,4> >;

  using MoleculeConnectivity = exaStamp::AtomBondConnectivity;
  
  // Bonds acceleration structures
  
  struct ForceEnergy
  {
    Vec3d m_force;
    double m_energy;
  };
  
  struct BondComputeBuffer
  {
    std::vector< ForceEnergy > m_force_enegery;
    std::vector< Mat3d > m_virial;
  };

# ifdef XSTAMP_USE_BEND_FORCE_BUFFER
  struct Force2Energy
  {
    std::array<Vec3d,2> m_force;
    double m_energy;
  };
  struct BendComputeBuffer
  {
    std::vector< Force2Energy > m_force_enegery;
    std::vector< std::array<Mat3d,2> > m_virial;
  };
# endif

  // a chemical chain is either a bond, bending, or torsion/improper
  // encoding is : particle index (in cell) , number of attached bonds , bond_idx/pos_in_bond
  using GridCellAtomsToChemicalChain = std::vector< uint64_t >;
  using GridAtomsToChemicalChain = std::vector< GridCellAtomsToChemicalChain >;

} // namespace exaStamp
