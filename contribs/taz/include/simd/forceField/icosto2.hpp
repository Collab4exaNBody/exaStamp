/*
 * icosto2.hpp
 *
 *  Created on: Oct 9, 2017
 *      Author: giarda
 */

/// @file
/// @brief Improper torsion force computer where E=C0+C1*cospsi+C2*cos2psi

#ifndef ICOSTO2_HPP_
#define ICOSTO2_HPP_


#include "simd/forceField/improperDerivation.hpp"
#include "simd/polynoms.hpp"


/// @brief Namespace for SIMD implementations
namespace simd {


  /// @brief Namespace for vectorized computations
	namespace kernels {

		/// @brief Improper torsion force computer where E=C0+C1*cospsi+C2*cos2psi
		template <class T> class icosto2 {

		public :

			void operator() (T* fx2_, T* fy2_, T* fz2_, T* fx3_, T* fy3_, T* fz3_, T* fx4_, T* fy4_, T* fz4_, T* e_,
							const T* rx21_, const T* ry21_, const T* rz21_, const T* rx31_, const T* ry31_, const T* rz31_, const T* rx41_, const T* ry41_, const T* rz41_,
							const T* c0_, const T* c1_, const T* c2_, const uint size) const;


		private :

  		static const vector_t<T> v0; ///< Vector of 0
  		static const vector_t<T> v1; ///< Vector of 1
  		static const vector_t<T> v2; ///< Vector of 2
    	static const vector_t<T> vm1; ///< Vector of -1
    	static const vector_t<T> v3rd; ///< Vector of 1/3

		};

		template <class T> const vector_t<T> icosto2<T>::v0  	(  0.);
		template <class T> const vector_t<T> icosto2<T>::v1  	(  1.);
		template <class T> const vector_t<T> icosto2<T>::v2  	(  2.);
		template <class T> const vector_t<T> icosto2<T>::vm1 	( -1.);
		template <class T> const vector_t<T> icosto2<T>::v3rd (1/3.);


		/// @brief Function call operator : compute the improper torsion forces and energies from the distances
		///
		/// Interaction between three particles 1, 2, 3 and 4, where 1 is the central one
		/// @param [out] fx2_ X components of the forces on 2
		/// @param [out] fy2_ Y components of the forces on 2
		/// @param [out] fz2_ Z components of the forces on 2
		/// @param [out] fx3_ X components of the forces on 3
		/// @param [out] fy3_ Y components of the forces on 3
		/// @param [out] fz3_ Z components of the forces on 3
		/// @param [out] fx4_ X components of the forces on 4
		/// @param [out] fy4_ Y components of the forces on 4
		/// @param [out] fz4_ Z components of the forces on 4
		/// @param [out] e_ Energies
		/// @param [in] rx21_ X components of the distance vector between 1 and 2
		/// @param [in] ry21_ Y components of the distance vector between 1 and 2
		/// @param [in] rz21_ Z components of the distance vector between 1 and 2
		/// @param [in] rx31_ X components of the distance vector between 1 and 3
		/// @param [in] ry31_ Y components of the distance vector between 1 and 3
		/// @param [in] rz31_ Z components of the distance vector between 1 and 3
		/// @param [in] rx41_ X components of the distance vector between 1 and 4
		/// @param [in] ry41_ Y components of the distance vector between 1 and 4
		/// @param [in] rz41_ Z components of the distance vector between 1 and 4
		/// @param [in] c0_ Constants 0
		/// @param [in] c1_ Constants 1
		/// @param [in] c2_ Constants 2
		/// @param [in] size Number of elements
		template <class T> void icosto2<T>::operator() (T* fx2_, T* fy2_, T* fz2_, T* fx3_, T* fy3_, T* fz3_, T* fx4_, T* fy4_, T* fz4_, T* e_,
						const T* rx21_, const T* ry21_, const T* rz21_, const T* rx31_, const T* ry31_, const T* rz31_, const T* rx41_, const T* ry41_, const T* rz41_,
						const T* c0_, const T* c1_, const T* c2_, const uint size) const {
		    	vector_t<T> rx21, ry21, rz21, rx31, ry31, rz31, rx41, ry41, rz41;
		    	vector_t<T> c0, c1, c2, e, de, cos;
		    	vector_t<T> fx2, fy2, fz2, fx3, fy3, fz3, fx4, fy4, fz4;
		    	vector_t<T> xn, yn, zn, xp, yp, zp, dxp, dyp, dzp;
		    	vector_t<T> dnxdx2, dnxdy2, dnxdz2, dnxdx3, dnxdy3, dnxdz3, dnxdx4, dnxdy4, dnxdz4;
		    	vector_t<T> dnydx2, dnydy2, dnydz2, dnydx3, dnydy3, dnydz3, dnydx4, dnydy4, dnydz4;
		    	vector_t<T> dnzdx2, dnzdy2, dnzdz2, dnzdx3, dnzdy3, dnzdz3, dnzdx4, dnzdy4, dnzdz4;
		    	for (uint i=0; i<size; i+=vector_t<T>::chunk_size) {
		    		// Load data
		    		rx21.load(rx21_+i);
		    		ry21.load(ry21_+i);
		    		rz21.load(rz21_+i);
		    		rx31.load(rx31_+i);
		    		ry31.load(ry31_+i);
		    		rz31.load(rz31_+i);
		    		rx41.load(rx41_+i);
		    		ry41.load(ry41_+i);
		    		rz41.load(rz41_+i);
		    		c0.load(c0_+i);
		    		c1.load(c1_+i);
		    		c2.load(c2_+i);
		    		// Coefficients for the energy calculation stored in c0, c1, c2
		    		c0=c0-c2;
		    		c2=v2*c2;

		    		// N234 = 23*34 = (21-31)*(31-41)
		    		// n234 = N234 / ||N234||
		    		// Compute N234
		    		cross(rx21-rx31,ry21-ry31,rz21-rz31,rx31-rx41,ry31-ry41,rz31-rz41,xn,yn,zn);

		    		// Compute the derivative of n234 relative to x2
		    		//	 Derivative of N234 = 23*34
		    		dnxdx2=v0;
		    		dnydx2=rz31-rz41;
		    		dnzdx2=ry41-ry31;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdx2,dnydx2,dnzdx2);
		    		// Compute the derivative of n234 relative to y2
		    		//	 Derivative of N234 = 23*34
		    		dnxdy2=rz41-rz31;
		    		dnydy2=v0;
		    		dnzdy2=rx31-rx41;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdy2,dnydy2,dnzdy2);
		    		// Compute the derivative of n234 relative to z2
		    		//	 Derivative of N234 = 23*34
		    		dnxdz2=ry31-ry41;
		    		dnydz2=rx41-rx31;
		    		dnzdz2=v0;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdz2,dnydz2,dnzdz2);
		    		// Compute the derivative of n234 relative to x3
		    		//	 Derivative of N234 = 34*42
		    		dnxdx3=v0;
		    		dnydx3=rz41-rz21;
		    		dnzdx3=ry21-ry41;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdx3,dnydx3,dnzdx3);
		    		// Compute the derivative of n234 relative to y3
		    		//	 Derivative of N234 = 34*42
		    		dnxdy3=rz21-rz41;
		    		dnydy3=v0;
		    		dnzdy3=rx41-rx21;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdy3,dnydy3,dnzdy3);
		    		// Compute the derivative of n234 relative to z3
		    		//	 Derivative of N234 = 34*42
		    		dnxdz3=ry41-ry21;
		    		dnydz3=rx21-rx41;
		    		dnzdz3=v0;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdz3,dnydz3,dnzdz3);
		    		// Compute the derivative of n234 relative to x4
		    		//	 Derivative of N234 = 42*23
		    		dnxdx4=v0;
		    		dnydx4=rz21-rz31;
		    		dnzdx4=ry31-ry21;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdx4,dnydx4,dnzdx4);
		    		// Compute the derivative of n234 relative to y4
		    		//	 Derivative of N234 = 42*23
		    		dnxdy4=rz31-rz21;
		    		dnydy4=v0;
		    		dnzdy4=rx21-rx31;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdy4,dnydy4,dnzdy4);
		    		// Compute the derivative of n234 relative to z2
		    		//	 Derivative of N234 = 42*23
		    		dnxdz4=ry21-ry31;
		    		dnydz4=rx31-rx21;
		    		dnzdz4=v0;
		    		//	 Derivative of the unit vector
		    		derivativeUnit(xn,yn,zn,dnxdz4,dnydz4,dnzdz4);

		    		// Then compute n234
		    		unit(xn,yn,zn);

		    		// Energy and forces for the angle between 41 and 234
		    		//	Project 41 on n234
		    		proj(rx41,ry41,rz41,xn,yn,zn,xp,yp,zp);
		    		//	omega is the angle between 41 and 41+proj
		    		cos=cosine(rx41,ry41,rz41,rx41-xp,ry41-yp,rz41-zp);
		    		//	Calculate E
		    		//		E=c0+c1*cos+c2*cos2 with cos2 is a polynom of cos
		    		//	Energy stored in e and energy derivative stored in tmp1
		    		e=Polynoms<T>::UpTo3(cos, c0, c1, c2, v0);
		    		de=vm1*Polynoms<T>::Derivative3(cos, c1, c2, v0);
		    		//	Compute forces
		    		//		Compute projected vector derivatives to x2
		    		dxp=vm1*scalar(rx41,ry41,rz41,dnxdx2,dnydx2,dnzdx2)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdx2;
		    		dyp=vm1*scalar(rx41,ry41,rz41,dnxdx2,dnydx2,dnzdx2)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydx2;
		    		dzp=vm1*scalar(rx41,ry41,rz41,dnxdx2,dnydx2,dnzdx2)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdx2;
		    		//		Compute fx2
		    		fx2=derivativeCosine1Variable(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to y2
		    		dxp=vm1*scalar(rx41,ry41,rz41,dnxdy2,dnydy2,dnzdy2)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdy2;
		    		dyp=vm1*scalar(rx41,ry41,rz41,dnxdy2,dnydy2,dnzdy2)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydy2;
		    		dzp=vm1*scalar(rx41,ry41,rz41,dnxdy2,dnydy2,dnzdy2)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdy2;
		    		//		Compute fy2
		    		fy2=derivativeCosine1Variable(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to z2
		    		dxp=vm1*scalar(rx41,ry41,rz41,dnxdz2,dnydz2,dnzdz2)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdz2;
		    		dyp=vm1*scalar(rx41,ry41,rz41,dnxdz2,dnydz2,dnzdz2)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydz2;
		    		dzp=vm1*scalar(rx41,ry41,rz41,dnxdz2,dnydz2,dnzdz2)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdz2;
		    		//		Compute fz2
		    		fz2=derivativeCosine1Variable(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to x3
		    		dxp=vm1*scalar(rx41,ry41,rz41,dnxdx3,dnydx3,dnzdx3)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdx3;
		    		dyp=vm1*scalar(rx41,ry41,rz41,dnxdx3,dnydx3,dnzdx3)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydx3;
		    		dzp=vm1*scalar(rx41,ry41,rz41,dnxdx3,dnydx3,dnzdx3)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdx3;
		    		//		Compute fx3
		    		fx3=derivativeCosine1Variable(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to y3
		    		dxp=vm1*scalar(rx41,ry41,rz41,dnxdy3,dnydy3,dnzdy3)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdy3;
		    		dyp=vm1*scalar(rx41,ry41,rz41,dnxdy3,dnydy3,dnzdy3)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydy3;
		    		dzp=vm1*scalar(rx41,ry41,rz41,dnxdy3,dnydy3,dnzdy3)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdy3;
		    		//		Compute fy3
		    		fy3=derivativeCosine1Variable(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to z3
		    		dxp=vm1*scalar(rx41,ry41,rz41,dnxdz3,dnydz3,dnzdz3)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdz3;
		    		dyp=vm1*scalar(rx41,ry41,rz41,dnxdz3,dnydz3,dnzdz3)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydz3;
		    		dzp=vm1*scalar(rx41,ry41,rz41,dnxdz3,dnydz3,dnzdz3)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdz3;
		    		//		Compute fz3
		    		fz3=derivativeCosine1Variable(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to x4
		    		dxp=xn*xn+vm1-scalar(rx41,ry41,rz41,dnxdx4,dnydx4,dnzdx4)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdx4;
		    		dyp=xn*yn-scalar(rx41,ry41,rz41,dnxdx4,dnydx4,dnzdx4)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydx4;
		    		dzp=xn*zn-scalar(rx41,ry41,rz41,dnxdx4,dnydx4,dnzdx4)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdx4;
		    		//		Compute fx4
		    		fx4=derivativeCosine2Variables(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,vm1,v0,v0,cos)*de;
		    		//		Compute projected vector derivatives to y4
		    		dxp=yn*xn-scalar(rx41,ry41,rz41,dnxdy4,dnydy4,dnzdy4)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdy4;
		    		dyp=yn*yn+vm1-scalar(rx41,ry41,rz41,dnxdy4,dnydy4,dnzdy4)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydy4;
		    		dzp=yn*zn-scalar(rx41,ry41,rz41,dnxdy4,dnydy4,dnzdy4)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdy4;
		    		//		Compute fy4
		    		fy4=derivativeCosine2Variables(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,v0,vm1,v0,cos)*de;
		    		//		Compute projected vector derivatives to z4
		    		dxp=zn*xn-scalar(rx41,ry41,rz41,dnxdz4,dnydz4,dnzdz4)*xn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnxdz4;
		    		dyp=zn*yn-scalar(rx41,ry41,rz41,dnxdz4,dnydz4,dnzdz4)*yn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnydz4;
		    		dzp=zn*zn+vm1-scalar(rx41,ry41,rz41,dnxdz4,dnydz4,dnzdz4)*zn-scalar(rx41,ry41,rz41,xn,yn,zn)*dnzdz4;
		    		//		Compute fz4
		    		fz4=derivativeCosine2Variables(rx41-xp,ry41-yp,rz41-zp,rx41,ry41,rz41,dxp,dyp,dzp,v0,v0,vm1,cos)*de;

		    		// Energy and forces for the angle between 31 and 234
		    		//	Project 31 on n234
		    		proj(rx31,ry31,rz31,xn,yn,zn,xp,yp,zp);
		    		//	omega is the angle between 31 and 31+proj
		    		cos=cosine(rx31,ry31,rz31,rx31-xp,ry31-yp,rz31-zp);
		    		//	Calculate E
		    		//		E=c0+c1*cos+c2*cos2 with cos2 is a polynom of cos
		    		//	Energy stored in e and energy derivative stored in tmp1
		    		e=e+Polynoms<T>::UpTo3(cos, c0, c1, c2, v0);
		    		de=vm1*Polynoms<T>::Derivative3(cos, c1, c2, v0);
		    		//	Compute forces
		    		//		Compute projected vector derivatives to x4
		    		dxp=vm1*scalar(rx31,ry31,rz31,dnxdx4,dnydx4,dnzdx4)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdx4;
		    		dyp=vm1*scalar(rx31,ry31,rz31,dnxdx4,dnydx4,dnzdx4)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydx4;
		    		dzp=vm1*scalar(rx31,ry31,rz31,dnxdx4,dnydx4,dnzdx4)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdx4;
		    		//		Compute fx4
		    		fx4=fx4+derivativeCosine1Variable(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to y4
		    		dxp=vm1*scalar(rx31,ry31,rz31,dnxdy4,dnydy4,dnzdy4)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdy4;
		    		dyp=vm1*scalar(rx31,ry31,rz31,dnxdy4,dnydy4,dnzdy4)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydy4;
		    		dzp=vm1*scalar(rx31,ry31,rz31,dnxdy4,dnydy4,dnzdy4)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdy4;
		    		//		Compute fy4
		    		fy4=fy4+derivativeCosine1Variable(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to z4
		    		dxp=vm1*scalar(rx31,ry31,rz31,dnxdz4,dnydz4,dnzdz4)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdz4;
		    		dyp=vm1*scalar(rx31,ry31,rz31,dnxdz4,dnydz4,dnzdz4)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydz4;
		    		dzp=vm1*scalar(rx31,ry31,rz31,dnxdz4,dnydz4,dnzdz4)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdz4;
		    		//		Compute fz4
		    		fz4=fz4+derivativeCosine1Variable(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to x2
		    		dxp=vm1*scalar(rx31,ry31,rz31,dnxdx2,dnydx2,dnzdx2)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdx2;
		    		dyp=vm1*scalar(rx31,ry31,rz31,dnxdx2,dnydx2,dnzdx2)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydx2;
		    		dzp=vm1*scalar(rx31,ry31,rz31,dnxdx2,dnydx2,dnzdx2)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdx2;
		    		//		Compute fx2
		    		fx2=fx2+derivativeCosine1Variable(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to y2
		    		dxp=vm1*scalar(rx31,ry31,rz31,dnxdy2,dnydy2,dnzdy2)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdy2;
		    		dyp=vm1*scalar(rx31,ry31,rz31,dnxdy2,dnydy2,dnzdy2)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydy2;
		    		dzp=vm1*scalar(rx31,ry31,rz31,dnxdy2,dnydy2,dnzdy2)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdy2;
		    		//		Compute fy2
		    		fy2=fy2+derivativeCosine1Variable(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to z2
		    		dxp=vm1*scalar(rx31,ry31,rz31,dnxdz2,dnydz2,dnzdz2)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdz2;
		    		dyp=vm1*scalar(rx31,ry31,rz31,dnxdz2,dnydz2,dnzdz2)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydz2;
		    		dzp=vm1*scalar(rx31,ry31,rz31,dnxdz2,dnydz2,dnzdz2)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdz2;
		    		//		Compute fz2
		    		fz2=fz2+derivativeCosine1Variable(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to x3
		    		dxp=xn*xn+vm1-scalar(rx31,ry31,rz31,dnxdx3,dnydx3,dnzdx3)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdx3;
		    		dyp=xn*yn-scalar(rx31,ry31,rz31,dnxdx3,dnydx3,dnzdx3)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydx3;
		    		dzp=xn*zn-scalar(rx31,ry31,rz31,dnxdx3,dnydx3,dnzdx3)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdx3;
		    		//		Compute fx3
		    		fx3=fx3+derivativeCosine2Variables(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,vm1,v0,v0,cos)*de;
		    		//		Compute projected vector derivatives to y3
		    		dxp=yn*xn-scalar(rx31,ry31,rz31,dnxdy3,dnydy3,dnzdy3)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdy3;
		    		dyp=yn*yn+vm1-scalar(rx31,ry31,rz31,dnxdy3,dnydy3,dnzdy3)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydy3;
		    		dzp=yn*zn-scalar(rx31,ry31,rz31,dnxdy3,dnydy3,dnzdy3)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdy3;
		    		//		Compute fy3
		    		fy3=fy3+derivativeCosine2Variables(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,v0,vm1,v0,cos)*de;
		    		//		Compute projected vector derivatives to z3
		    		dxp=zn*xn-scalar(rx31,ry31,rz31,dnxdz3,dnydz3,dnzdz3)*xn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnxdz3;
		    		dyp=zn*yn-scalar(rx31,ry31,rz31,dnxdz3,dnydz3,dnzdz3)*yn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnydz3;
		    		dzp=zn*zn+vm1-scalar(rx31,ry31,rz31,dnxdz3,dnydz3,dnzdz3)*zn-scalar(rx31,ry31,rz31,xn,yn,zn)*dnzdz3;
		    		//		Compute fz3
		    		fz3=fz3+derivativeCosine2Variables(rx31-xp,ry31-yp,rz31-zp,rx31,ry31,rz31,dxp,dyp,dzp,v0,v0,vm1,cos)*de;

		    		// Energy and forces for the angle between 21 and 234
		    		//	Project 21 on n234
		    		proj(rx21,ry21,rz21,xn,yn,zn,xp,yp,zp);
		    		//	omega is the angle between 21 and 21+proj
		    		cos=cosine(rx21,ry21,rz21,rx21-xp,ry21-yp,rz21-zp);
		    		//	Calculate E
		    		//		E=c0+c1*cos+c2*cos2 with cos2 is a polynom of cos
		    		//	Energy stored in e and energy derivative stored in tmp1
		    		e=e+Polynoms<T>::UpTo3(cos, c0, c1, c2, v0);
		    		de=vm1*Polynoms<T>::Derivative3(cos, c1, c2, v0);
		    		//	Compute forces
		    		//		Compute projected vector derivatives to x3
		    		dxp=vm1*scalar(rx21,ry21,rz21,dnxdx3,dnydx3,dnzdx3)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdx3;
		    		dyp=vm1*scalar(rx21,ry21,rz21,dnxdx3,dnydx3,dnzdx3)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydx3;
		    		dzp=vm1*scalar(rx21,ry21,rz21,dnxdx3,dnydx3,dnzdx3)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdx3;
		    		//		Compute fx3
		    		fx3=fx3+derivativeCosine1Variable(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to y3
		    		dxp=vm1*scalar(rx21,ry21,rz21,dnxdy3,dnydy3,dnzdy3)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdy3;
		    		dyp=vm1*scalar(rx21,ry21,rz21,dnxdy3,dnydy3,dnzdy3)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydy3;
		    		dzp=vm1*scalar(rx21,ry21,rz21,dnxdy3,dnydy3,dnzdy3)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdy3;
		    		//		Compute fy3
		    		fy3=fy3+derivativeCosine1Variable(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to z3
		    		dxp=vm1*scalar(rx21,ry21,rz21,dnxdz3,dnydz3,dnzdz3)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdz3;
		    		dyp=vm1*scalar(rx21,ry21,rz21,dnxdz3,dnydz3,dnzdz3)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydz3;
		    		dzp=vm1*scalar(rx21,ry21,rz21,dnxdz3,dnydz3,dnzdz3)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdz3;
		    		//		Compute fz3
		    		fz3=fz3+derivativeCosine1Variable(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to x4
		    		dxp=vm1*scalar(rx21,ry21,rz21,dnxdx4,dnydx4,dnzdx4)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdx4;
		    		dyp=vm1*scalar(rx21,ry21,rz21,dnxdx4,dnydx4,dnzdx4)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydx4;
		    		dzp=vm1*scalar(rx21,ry21,rz21,dnxdx4,dnydx4,dnzdx4)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdx4;
		    		//		Compute fx4
		    		fx4=fx4+derivativeCosine1Variable(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to y4
		    		dxp=vm1*scalar(rx21,ry21,rz21,dnxdy4,dnydy4,dnzdy4)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdy4;
		    		dyp=vm1*scalar(rx21,ry21,rz21,dnxdy4,dnydy4,dnzdy4)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydy4;
		    		dzp=vm1*scalar(rx21,ry21,rz21,dnxdy4,dnydy4,dnzdy4)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdy4;
		    		//		Compute fy4
		    		fy4=fy4+derivativeCosine1Variable(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to z4
		    		dxp=vm1*scalar(rx21,ry21,rz21,dnxdz4,dnydz4,dnzdz4)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdz4;
		    		dyp=vm1*scalar(rx21,ry21,rz21,dnxdz4,dnydz4,dnzdz4)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydz4;
		    		dzp=vm1*scalar(rx21,ry21,rz21,dnxdz4,dnydz4,dnzdz4)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdz4;
		    		//		Compute fz4
		    		fz4=fz4+derivativeCosine1Variable(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,cos)*de;
		    		//		Compute projected vector derivatives to x2
		    		dxp=xn*xn+vm1-scalar(rx21,ry21,rz21,dnxdx2,dnydx2,dnzdx2)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdx2;
		    		dyp=xn*yn-scalar(rx21,ry21,rz21,dnxdx2,dnydx2,dnzdx2)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydx2;
		    		dzp=xn*zn-scalar(rx21,ry21,rz21,dnxdx2,dnydx2,dnzdx2)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdx2;
		    		//		Compute fx2
		    		fx2=fx2+derivativeCosine2Variables(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,vm1,v0,v0,cos)*de;
		    		//		Compute projected vector derivatives to y2
		    		dxp=yn*xn-scalar(rx21,ry21,rz21,dnxdy2,dnydy2,dnzdy2)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdy2;
		    		dyp=yn*yn+vm1-scalar(rx21,ry21,rz21,dnxdy2,dnydy2,dnzdy2)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydy2;
		    		dzp=yn*zn-scalar(rx21,ry21,rz21,dnxdy2,dnydy2,dnzdy2)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdy2;
		    		//		Compute fy2
		    		fy2=fy2+derivativeCosine2Variables(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,v0,vm1,v0,cos)*de;
		    		//		Compute projected vector derivatives to z2
		    		dxp=zn*xn-scalar(rx21,ry21,rz21,dnxdz2,dnydz2,dnzdz2)*xn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnxdz2;
		    		dyp=zn*yn-scalar(rx21,ry21,rz21,dnxdz2,dnydz2,dnzdz2)*yn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnydz2;
		    		dzp=zn*zn+vm1-scalar(rx21,ry21,rz21,dnxdz2,dnydz2,dnzdz2)*zn-scalar(rx21,ry21,rz21,xn,yn,zn)*dnzdz2;
		    		//		Compute fz2
		    		fz2=fz2+derivativeCosine2Variables(rx21-xp,ry21-yp,rz21-zp,rx21,ry21,rz21,dxp,dyp,dzp,v0,v0,vm1,cos)*de;

		    		// Average of the 3
		    		e=e*v3rd;
		    		fx2=fx2*v3rd;
		    		fy2=fy2*v3rd;
		    		fz2=fz2*v3rd;
		    		fx3=fx3*v3rd;
		    		fy3=fy3*v3rd;
		    		fz3=fz3*v3rd;
		    		fx4=fx4*v3rd;
		    		fy4=fy4*v3rd;
		    		fz4=fz4*v3rd;
		    		// Store results
		    		e.store(e_+i);
		    		fx2.store(fx2_+i);
		    		fy2.store(fy2_+i);
		    		fz2.store(fz2_+i);
		    		fx3.store(fx3_+i);
		    		fy3.store(fy3_+i);
		    		fz3.store(fz3_+i);
		    		fx4.store(fx4_+i);
		    		fy4.store(fy4_+i);
		    		fz4.store(fz4_+i);
		    	}
				}


	}  // namespace kernels


}  // namespace simd

#endif /* ICOSTO2_HPP_ */
