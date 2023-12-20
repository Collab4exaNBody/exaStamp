/*
 * dcosn.hpp
 *
 *  Created on: Oct 9, 2017
 *      Author: giarda
 */

/// @file
/// @brief Dihedral angle force computer where E=C0+C1*cos(n*phi), n=1,2,3 or 6

#ifndef DCOSN_HPP_
#define DCOSN_HPP_


#include "simd/forceField/dihedralDerivation.hpp"
#include "simd/polynoms.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
	namespace kernels {

		/// @brief Dihedral angle force computer where E=C0+C1*cos(n*phi), n=1,2,3 or 6
		template <class T> class dcosn {

		public :

			void operator() (T* fx1_, T* fy1_, T* fz1_, T* fx2_, T* fy2_, T* fz2_, T* fx4_, T* fy4_, T* fz4_, T* e_,
					const T* rx21_, const T* ry21_, const T* rz21_, const T* rx23_, const T* ry23_, const T* rz23_, const T* rx34_, const T* ry34_, const T* rz34_,
					const T* c0_, const T* c1_, const T* n_, const uint size) const;

		private :

  		static const vector_t<T> v0; ///< Vector of 0
  		static const vector_t<T> v1; ///< Vector of 1
  		static const vector_t<T> v2; ///< Vector of 2
    	static const vector_t<T> v3; ///< Vector of 3
    	static const vector_t<T> v4; ///< Vector of 4
    	static const vector_t<T> v6; ///< Vector of 6
    	static const vector_t<T> v8; ///< Vector of 8
    	static const vector_t<T> vm1; ///< Vector of -1
    	static const vector_t<T> vhalf; ///< Vector of 1/2
    	static const vector_t<T> v3rd; ///< Vector of 1/3
    	static const vector_t<T> v4th; ///< Vector of 1/4
    	static const vector_t<T> v5th; ///< Vector of 1/5
    	static const vector_t<T> v6th; ///< Vector of 1/6
    	static const vector_t<T> v10th; ///< Vector of 1/10

		};

		template <class T> const vector_t<T> dcosn<T>::v0  	(  0.);
		template <class T> const vector_t<T> dcosn<T>::v1  	(  1.);
		template <class T> const vector_t<T> dcosn<T>::v2  	(  2.);
		template <class T> const vector_t<T> dcosn<T>::v3  	(  3.);
		template <class T> const vector_t<T> dcosn<T>::v4  	(  4.);
		template <class T> const vector_t<T> dcosn<T>::v6  	(  6.);
		template <class T> const vector_t<T> dcosn<T>::v8  	(  8.);
		template <class T> const vector_t<T> dcosn<T>::vm1  (	-1.);
		template <class T> const vector_t<T> dcosn<T>::vhalf(1/2.);
		template <class T> const vector_t<T> dcosn<T>::v3rd (1/3.);
		template <class T> const vector_t<T> dcosn<T>::v4th (1/4.);
		template <class T> const vector_t<T> dcosn<T>::v5th (1/5.);
		template <class T> const vector_t<T> dcosn<T>::v6th (1/6.);
		template <class T> const vector_t<T> dcosn<T>::v10th (1/10.);


		/// @brief Function call operator : compute the dihedral angle forces and energies from the distances
		///
		/// Interaction between four particles 1, 2, 3 and 4, where 2 and 3 are central
		/// @param [out] fx1_ X components of the forces on 1
		/// @param [out] fy1_ Y components of the forces on 1
		/// @param [out] fz1_ Z components of the forces on 1
		/// @param [out] fx2_ X components of the forces on 2
		/// @param [out] fy2_ Y components of the forces on 2
		/// @param [out] fz2_ Z components of the forces on 2
		/// @param [out] fx4_ X components of the forces on 4
		/// @param [out] fy4_ Y components of the forces on 4
		/// @param [out] fz4_ Z components of the forces on 4
		/// @param [out] e_ Energies
		/// @param [in] rx21_ X components of the distance vector between 1 and 2
		/// @param [in] ry21_ Y components of the distance vector between 1 and 2
		/// @param [in] rz21_ Z components of the distance vector between 1 and 2
		/// @param [in] rx23_ X components of the distance vector between 3 and 2
		/// @param [in] ry23_ Y components of the distance vector between 3 and 2
		/// @param [in] rz23_ Z components of the distance vector between 3 and 2
		/// @param [in] rx34_ X components of the distance vector between 4 and 3
		/// @param [in] ry34_ Y components of the distance vector between 4 and 3
		/// @param [in] rz34_ Z components of the distance vector between 4 and 3
		/// @param [in] c0_ Constants 0
		/// @param [in] c1_ Constants 1
		/// @param [in] n_ Polynom degree
		/// @param [in] size Number of elements
		template <class T> void dcosn<T>::operator() (T* fx1_, T* fy1_, T* fz1_, T* fx2_, T* fy2_, T* fz2_, T* fx4_, T* fy4_, T* fz4_, T* e_,
				const T* rx21_, const T* ry21_, const T* rz21_, const T* rx23_, const T* ry23_, const T* rz23_, const T* rx34_, const T* ry34_, const T* rz34_,
				const T* c0_, const T* c1_, const T* n_, const uint size) const {
    	vector_t<T> rx21, ry21, rz21, rx23, ry23, rz23, rx34, ry34, rz34, n;
    	vector_t<T> tmp0, fx, fy, fz;
    	vector_t<T> tmp1, tmp2, tmp3, tmp4, tmp6, cos;
    	vector_t<T> xA, yA, zA, xB, yB, zB, dxA, dyA, dzA, dxB, dyB, dzB;
    	for (uint i=0; i<size; i+=vector_t<T>::chunk_size) {
    		// Load data
    		rx21.load(rx21_+i);
    		ry21.load(ry21_+i);
    		rz21.load(rz21_+i);
    		rx23.load(rx23_+i);
    		ry23.load(ry23_+i);
    		rz23.load(rz23_+i);
    		rx34.load(rx34_+i);
    		ry34.load(ry34_+i);
    		rz34.load(rz34_+i);
    		tmp0.load(c0_+i);
    		tmp1.load(c1_+i);
    		n.load(n_+i);
    		// A = 21*23
    		cross(rx21,ry21,rz21,rx23,ry23,rz23,xA,yA,zA);
    		// B= 32*34 = 34*23
    		cross(rx34,ry34,rz34,rx23,ry23,rz23,xB,yB,zB);
    		// Cosine of the dihedral angle
    		cos=cosine(xA,yA,zA,xB,yB,zB);
    		// Calculate E
    		//	 E=c0+c1*cosn with cosn a polynom of cos
    		// Coefficients stored in tmp0, tmp1, tmp2, tmp3, tmp4 and tmp6
    		tmp0=tmp0+tmp1*((v1-n)*(v3-n)*(v6-n)*v4th+(v1-n)*(v2-n)*(v3-n)*v6th*v10th);
    		tmp2=tmp1*((n-v1)*(n-v3)*(n-v6)*vhalf+(n-v1)*(n-v2)*(n-v3)*v3*v10th);
    		tmp3=tmp1*(v1-n)*(v2-n)*(v6-n)*v2*v3rd;
    		tmp4=tmp1*(v1-n)*(v2-n)*(v3-n)*v4*v5th;
    		tmp6=tmp1*(n-v1)*(n-v2)*(n-v3)*v8*v3rd*v5th;
    		tmp1=tmp1*((v2-n)*(v3-n)*(v6-n)*v10th+(n-v1)*(n-v2)*(n-v6)*vhalf);
    		// Energy stored in tmp0 and energy derivative stored in tmp1
    		tmp0=Polynoms<T>::UpTo6(cos, tmp0, tmp1, tmp2, tmp3, tmp4, v0, tmp6);
    		tmp1=vm1*Polynoms<T>::Derivative6(cos, tmp1, tmp2, tmp3, tmp4, v0, tmp6);
    		// Store result
    		tmp0.store(e_+i);
    		// Calculate f1
    		// 	f1x = dE/dx1 = dcos/dx1 * dE/dcos
    		fx=derivativeCosine1Variable(xA,yA,zA,xB,yB,zB,v0,v0-rz23,ry23,cos)*tmp1;
    		fy=derivativeCosine1Variable(xA,yA,zA,xB,yB,zB,rz23,v0,v0-rx23,cos)*tmp1;
    		fz=derivativeCosine1Variable(xA,yA,zA,xB,yB,zB,v0-ry23,rx23,v0,cos)*tmp1;
    		// Store result
    		fx.store(fx1_+i);
    		fy.store(fy1_+i);
    		fz.store(fz1_+i);
    		// Calculate f4
    		fx=derivativeCosine1Variable(xB,yB,zB,xA,yA,zA,v0,v0-rz23,ry23,cos)*tmp1;
    		fy=derivativeCosine1Variable(xB,yB,zB,xA,yA,zA,rz23,v0,v0-rx23,cos)*tmp1;
    		fz=derivativeCosine1Variable(xB,yB,zB,xA,yA,zA,v0-ry23,rx23,v0,cos)*tmp1;
    		// Store result
    		fx.store(fx4_+i);
    		fy.store(fy4_+i);
    		fz.store(fz4_+i);
    		// Calculate f2
    		fx=derivativeCosine2Variables(xA,yA,zA,xB,yB,zB,v0,rz23-rz21,ry21-ry23,v0,v0-rz34,ry34,cos)*tmp1;
    		fy=derivativeCosine2Variables(xA,yA,zA,xB,yB,zB,rz21-rz23,v0,rx23-rx21,rz34,v0,v0-rx34,cos)*tmp1;
    		fz=derivativeCosine2Variables(xA,yA,zA,xB,yB,zB,ry23-ry21,rx21-rx23,v0,v0-ry34,rx34,v0,cos)*tmp1;
    		// Store result
    		fx.store(fx2_+i);
    		fy.store(fy2_+i);
    		fz.store(fz2_+i);
    	}
		}


	}  // namespace kernels


}  // namespace simd

#endif /* DCOSN_HPP_ */
