/// @file 
/// @brief Definition and implementation of distance corrector (in case of periodic boundary conditions)

#ifndef __CORRECTER_HPP_INCLUDED
#define __CORRECTER_HPP_INCLUDED


#include "simdAMR/correct.hpp"

#include "utils/vec3/vec3.hpp"


/// @brief Distance corrector for boundary conditions case
class Correcter {

public:

  /// @brief Default constructor
  Correcter() {}

  /// @brief Constructor with data (not used)
  /// @param [in] setup Dimensions in which distance correction is necessary
  /// @param [in] extension_ System extension
  Correcter(const vec3<bool>& setup, const vec3<double>& extension_)
    : dim(setup), extension(extension_), extensionHalf(0.5*extension), minusExtensionHalf(-0.5*extension) {}

  /// @brief Destructor (nothing to do)
  ~Correcter() {}

  /// @brief Set Correcter from parameters
  /// @param [in] setup Dimensions in which distance correction is necessary
  /// @param [in] extension_ System extension
  inline void set(const vec3<bool>& setup, const vec3<double>& extension_) {
    dim                = setup;
    extension          = extension_;
    extensionHalf      = 0.5*extension;
    minusExtensionHalf = -1.*extensionHalf;
  }

  /// @brief Correct a single distance
  /// @param [in,out] distance Distance to correct
  inline void correct(vec3<double>& distance) const {
    if (dim.x) correct<0>(distance.x);
    if (dim.y) correct<1>(distance.y);
    if (dim.z) correct<2>(distance.z);
  }

  /// @brief Correct multiple distances, call simd version
  /// @param [in] size Number of distances to correct
  /// @param [in,out] rx Pointer to x components
  /// @param [in,out] ry Pointer to y components
  /// @param [in,out] rz Pointer to z components
  inline void correct(const uint size, double* rx, double* ry, double* rz) const {
    if (dim.x) simdAMR::kernels::correct(rx, extension[0], size); 
    if (dim.y) simdAMR::kernels::correct(ry, extension[1], size); 
    if (dim.z) simdAMR::kernels::correct(rz, extension[2], size); 
  }

  /// @brief Correct the position of a particle to put it near another
  /// @param [in] position1 Position to compare to
  /// @param [in,out] position2 Position to correct
  inline void correctP(const vec3<double>& position1, vec3<double> position2) const {
    if (dim.x) correctP<0>(position1.x, position2.x);
    if (dim.y) correctP<1>(position1.y, position2.y);
    if (dim.z) correctP<2>(position1.z, position2.z);
  }

  /// @brief Correct multiple positions according to multiple others, call simd version
  /// @param [in] size Number of positions to correct
  /// @param [in] x1 Pointer to x components of the position to compare to
  /// @param [in] y1 Pointer to y components of the position to compare to
  /// @param [in] z1 Pointer to z components of the position to compare to
  /// @param [in,out] x2 Pointer to x components of the position to correct
  /// @param [in,out] y2 Pointer to e components of the position to correct
  /// @param [in,out] z2 Pointer to z components of the position to correct
  inline void correctP(const uint size, double* x1, double* y1, double* z1, double* x2, double* y2, double* z2) const {
    if (dim.x) simdAMR::kernels::correctP(x1, x2, extension[0], size);
    if (dim.y) simdAMR::kernels::correctP(y1, y2, extension[1], size);
    if (dim.z) simdAMR::kernels::correctP(z1, z2, extension[2], size);
  }
  
  
  /// @brief Correct multiple positions according to multiple others, call simd version
  /// @param [in] size Number of positions to correct
  /// @param [in] x1 Pointer to x components of the position to compare to
  /// @param [in] y1 Pointer to y components of the position to compare to
  /// @param [in] z1 Pointer to z components of the position to compare to
  /// @param [in,out] x2 Pointer to x components of the position to correct
  /// @param [in,out] y2 Pointer to e components of the position to correct
  /// @param [in,out] z2 Pointer to z components of the position to correct
  inline void correctP(const uint size, vec3<double>& position1, double* x2, double* y2, double* z2) const {
 
    double* ptr;
    double tmp, tmp1, tmp2;

    for(int d_ = 0 ; d_ < 3 ; ++d_)
      if(dim[d_])
      {
        const double D (extension[d_]);
        const double pd( 0.5*D);
        const double md(-0.5*D);   
        
        if(d_==0) ptr=x2;
        else if(d_==1) ptr=y2;
        else ptr=z2;
        
        #pragma omp simd
        //#pragma vector aligned
        for(uint i = 0 ; i < size ; ++i)
        {
          tmp = ptr[i]-position1[d_];
          tmp1 = tmp > pd ? D : 0;
          tmp2 = tmp < md ? D : 0;
      	  ptr[i] += tmp2 - tmp1;
        }
      }

  }


private:

  /// @brief Correct a component to a distance
  /// @tparam D Identifier of the dimension (1 for x, 2 for y and 3 for z)
  /// @param [in,out] distance Component to correct
  template <const uint8_t D> inline void correct(double& distance) const {
    if      (distance<minusExtensionHalf[D]) distance+=extension[D];
    else if (distance>     extensionHalf[D]) distance-=extension[D];
  }

  /// @brief Correct a component to a position
  /// @tparam D Identifier of the dimension (1 for x, 2 for y and 3 for z)
  /// @param [in] position1 Component to compare to
  /// @param [in,out] position2 Component to correct
  template <const uint8_t D> inline void correctP(const double& position1, double& position2) const {
    if      (position2-position1<minusExtensionHalf[D]) position2+=extension[D];
    else if (position2-position1>     extensionHalf[D]) position2-=extension[D];
  }

  vec3<bool>   dim; ///< Identify if the distances must be corrected in each dimension
  vec3<double> extension; ///< System extension
  vec3<double> extensionHalf; ///< Half of the system extension
  vec3<double> minusExtensionHalf; ///< Minus half of the system extension

};

#endif // __CORRECTER_HPP_INCLUDED
