/// @file
/// @brief Vectorization of neighbors computation

#ifndef __NEIGHBORS_AMR_HPP_INCLUDED
#define __NEIGHBORS_AMR_HPP_INCLUDED



/// @brief Namespace for SIMD implementations
namespace simdAMR {


  /// @brief Namespace for vectorized computations
  namespace kernels {


  	/// @brief Vectorized computation of the neighbors of a particle
  	/// @tparam T Type of the variables
  	/// @param [in] ri_ Position of the particle
  	/// @param [in] rx_ Y components of the positions of the potential neighbors
  	/// @param [in] ry_ Y components of the positions of the potential neighbors
  	/// @param [in] rz_ Z components of the positions of the potential neighbors
  	/// @param [in] rcut_ Cutoff radius for each neighbor
  	/// @param [out] out_ Boolean results (is or isn't a neighbor)
  	/// @param [in] n Number of potential neighbors
  	/// @warning The boolean results are stored into a array of T (double)
    inline void neighbors (const double *ri_, const double *rx_, const double *ry_, const double *rz_, const double *rcut_, double *out_, const uint n) {



      const double rxi =ri_[0];
      const double ryi =ri_[1];
      const double rzi =ri_[2];



      #pragma omp simd
      ////#pragma vector aligned
      for (int j=0; j<n; j++) {

       assert(rxi != rx_[j] || 
              ryi != ry_[j] || 
              rzi != rz_[j] );

      	double rxj = rxi - rx_[j];
      	double ryj = ryi - ry_[j];
      	double rzj = rzi - rz_[j];


        out_[j] = (rxj*rxj + ryj*ryj + rzj*rzj) < rcut_[j];

      }

    }


  }


}

#endif // __NEIGHBORS_AMR_HPP_INCLUDED
