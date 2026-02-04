/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <onika/math/basic_types_def.h>
#include <onika/math/basic_types_constructors.h>
#include <onika/math/basic_types_operators.h>

#include <algorithm>
#include <cmath>
#include <array>
#include <functional>
#include <onika/cuda/cuda.h>

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

namespace exaStamp
{

  using namespace exanb;

  ONIKA_HOST_DEVICE_FUNC inline int Jacobi ( const Mat3d FF, Mat3d& vv, Vec3d& dd)
  {

    int j,iq,ip,i, nb_trials = 0;
    double tresh,theta,tau,t,sm,s,h,g,c,b[3],z[3], F[3][3], v[3][3], d[3];

    F[0][0] = FF.m11;
    F[0][1] = FF.m12;
    F[0][2] = FF.m13;

    F[1][0] = FF.m21;
    F[1][1] = FF.m22;
    F[1][2] = FF.m23;

    F[2][0] = FF.m31;
    F[2][1] = FF.m32;
    F[2][2] = FF.m33;
    
    v[0][0] = v[1][1] = v[2][2] = 1.;
    v[0][1] = v[0][2] = v[1][0] = v[1][2] = v[2][0] = v[2][1] = 0.;

    /* b = d = diag of this. */
    for (ip=0;ip<3;ip++)
      {
	b[ip]=d[ip]=F[ip][ip];
	z[ip]=0.0;
      }

    for (i=0;i<50;i++)
      {
	sm=0.0; /* Sum of the off-diag terms. */
	for (ip=0;ip<2;ip++)
	  {
	    for (iq=ip+1;iq<3;iq++)
	      sm += fabs(F[ip][iq]);
	  }

	if (fabs(sm) <= 1e-12) goto sortie;	// Convergence

	if (i < 4)
	  tresh=0.2*sm/9.0;
	else
	  tresh=0.0;
	for (ip=0;ip<2;ip++)
	  {
	    for (iq=ip+1;iq<3;iq++)
	      {
		g=100.0*fabs(F[ip][iq]);
		if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
		    && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
		  F[ip][iq]=0.0;
		else if (fabs(F[ip][iq]) > tresh)
		  {
		    h=d[iq]-d[ip];
		    if ((float)(fabs(h)+g) == (float)fabs(h))
		      t=(F[ip][iq])/h;
		    else
		      {
			theta=0.5*h/(F[ip][iq]);
			t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
			if (theta < 0.0) t = -t;
		      }

		    /* Rotation Tensor_o2. */
		    c=1.0/sqrt(1.0+t*t);
		    s=t*c;
		    tau=s/(1.0+c);
		    h=t*F[ip][iq];
		    z[ip] -= h;
		    z[iq] += h;
		    d[ip] -= h;
		    d[iq] += h;
		    F[ip][iq]=0.0;
		    for (j=0;j<=ip-1;j++)     {ROTATE(F,  j, ip, j, iq)}
		    for (j=ip+1;j<=iq-1;j++)  {ROTATE(F, ip, j,  j, iq)}
		    for (j=iq+1;j<3;j++)      {ROTATE(F, ip, j,  iq, j)}
		    for (j=0;j<3;j++)         {ROTATE(v,     j, ip, j, iq)}
		    nb_trials += 1;
		  }
	      }
	  }
	for (ip=0;ip<3;ip++)
	  {
	    b[ip] += z[ip];
	    d[ip]=b[ip];
	    z[ip]=0.0;
	  }
      }/* for (i=0;i<50;i++) */    
    {
      
      lerr << "Diagonalization of a Tensor_o2 is not possible in less than " << nb_trials << " steps" << std::endl;
	
    }
  sortie:;
    
    /* Give back the off-diag components of a. */
    F[0][1] = F[1][0];
    F[0][2] = F[2][0];
    F[1][2] = F[2][1];

    vv.m11 = v[0][0];
    vv.m12 = v[0][1];
    vv.m13 = v[0][2];    

    vv.m21 = v[1][0];
    vv.m22 = v[1][1];
    vv.m23 = v[1][2];    

    vv.m31 = v[2][0];
    vv.m32 = v[2][1];
    vv.m33 = v[2][2];    

    vv = transpose(vv);

    dd.x = d[0];
    dd.y = d[1];
    dd.z = d[2];
    
    return nb_trials;
    
  }
  // @brief component polar RU decomposition of F, if an estimation of R is provided (R can be I)
  ONIKA_HOST_DEVICE_FUNC inline int RU_decomposition ( const Mat3d F, Mat3d& R, Mat3d& U)
  {
    int i,j,k,count;
    Vec3d lambda;
    Mat3d UU, vp;
    UU = transpose(F) * F;

    count = Jacobi(UU, vp, lambda);
    
    if (lambda.x < 1e-9 or lambda.y < 1e-9 or lambda.z < 1e-9)
      {
    	std::cout << "Negative eigen value for the square of a stretch strain!" << std::endl;
      }
    lambda.x = std::max(lambda.x, 0.0);
    lambda.y = std::max(lambda.y, 0.0);
    lambda.z = std::max(lambda.z, 0.0);

    double Umat[3][3]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double lambdavec[3]={lambda.x,lambda.y,lambda.z};
    double vpmat[3][3]={vp.m11,vp.m12,vp.m13,vp.m21,vp.m22,vp.m23,vp.m31,vp.m32,vp.m33};

    for (k=0; k<3; k++)
      for (i=0; i<3; i++)
	for (j=0; j<3; j++)
	  Umat[i][j] += sqrt(lambdavec[k]) * vpmat[k][i] * vpmat[k][j];

    U.m11 = Umat[0][0];
    U.m12 = Umat[0][1];
    U.m13 = Umat[0][2];    

    U.m21 = Umat[1][0];
    U.m22 = Umat[1][1];
    U.m23 = Umat[1][2];    

    U.m31 = Umat[2][0];
    U.m32 = Umat[2][1];
    U.m33 = Umat[2][2];    

    R = F * inverse(U);

    return count;
  }
  
}
