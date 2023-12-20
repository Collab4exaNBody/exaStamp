#pragma once




/// @brief Namespace for SIMD implementations
namespace simdAMR {


  /// @brief Namespace for vectorized computations
  namespace kernels {


  /// @brief Correct distances with vectorisation
  /// @tparam T Type of the distance to correct
  /// @param [in,out] a_ Distances to correct
  /// @param [in] D_ Extension of the system
  /// @param [in] n Number of distances to correct
  inline void correct (double* a_, const double& D_, const size_t& n) {
      
      double tmp1, tmp2;
      
      double D (     D_);
      double pd( 0.5*D_);
      double md(-0.5*D_);
      
      for (size_t i=0; i<n; i++) {
	
        tmp1 = a_[i] > pd ? D : 0;
        tmp2 = a_[i] < md ? D : 0;
      	a_[i] += tmp2 - tmp1;
	
      }
      
    }

  /// @brief Correct positions with vectorisation
  /// @tparam T Type of the distance to correct
  /// @param [in] ref_ Positions to compare to
  /// @param [in,out] p_ Positions to correct
  /// @param [in] D_ Extension of the system
  /// @param [in] n Number of distances to correct
  inline void correctP (double* ref_, double* p_, const double D_, const size_t n) {

      double tmpd, tmp1, tmp2;

      const double D (     D_);
      const double pd( 0.5*D_);
      const double md(-0.5*D_);

      for (size_t i=0; i<n; i++) {
      	tmpd=p_[i]-ref_[i];
        tmp1 = tmpd > pd ? D : 0;
        tmp2 = tmpd < md ? D : 0;
      	p_[i] += tmp2 - tmp1;
      }

    }

  }


}
  

