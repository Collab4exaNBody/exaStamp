#pragma once

#include <vector>

namespace exaStamp
{

  struct CellParticleLocalStructuralMetrics
  {
    // The size of vectors below equals the number of particles in cell
    // All variables are per-atom quantities

    // Structural metrics used as descriptors in ML models (Bispectrum for now)
    std::vector< std::vector<double> > bispectrum; // Bispectrum
    std::vector< std::vector<double> > steinhart;  // Steinhart

    // Crystal structure assigned using logistic regression (ML model)
    std::vector< int > crystal_structure;
    
  };

  using GridParticleLocalStructuralMetrics = std::vector< CellParticleLocalStructuralMetrics >;

}


