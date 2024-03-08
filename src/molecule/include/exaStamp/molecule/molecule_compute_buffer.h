#pragma once

#include <exaStamp/molecule/molecule_species.h>
#include <exanb/core/basic_types_def.h>

namespace exaStamp
{
  struct MoleculeComputeBuffer
  {
    double rx[MAX_MOLECULE_ATOMS];
    double ry[MAX_MOLECULE_ATOMS];
    double rz[MAX_MOLECULE_ATOMS];
    double fx[MAX_MOLECULE_ATOMS];
    double fy[MAX_MOLECULE_ATOMS];
    double fz[MAX_MOLECULE_ATOMS];
    double ep[MAX_MOLECULE_ATOMS];
    Mat3d virial[0];
  };

  struct MoleculeVirialComputeBuffer
  {
    double rx[MAX_MOLECULE_ATOMS];
    double ry[MAX_MOLECULE_ATOMS];
    double rz[MAX_MOLECULE_ATOMS];
    double fx[MAX_MOLECULE_ATOMS];
    double fy[MAX_MOLECULE_ATOMS];
    double fz[MAX_MOLECULE_ATOMS];
    double ep[MAX_MOLECULE_ATOMS];
    Mat3d virial[MAX_MOLECULE_ATOMS];
  };

}

