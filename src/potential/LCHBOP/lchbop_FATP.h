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

#include <cstdint>
#include <vector>

#include "lchbop_getEF.h"
//#include "lchbop_utils.h"
#include <onika/math/basic_types_def.h>
 

using std::vector; 

struct FATP
{
  double oneby3=1.0/3.0,sqrt3=1.732050807569;
  uint8_t inbcnb[2][2];
//  static size_t nbrm=400;
  size_t ijp[2];
  int imat[4][40],nfrc3[2];
//  uint16_t imat[4][40],nfrc3[2];
  double cf1,p2sm1,p2sm2;
  double fcnj[4][4][2];
  double xmelmin1[3],dxmel1[3],xmelmin2[3][3],dxmel2[3][3];
  double aab1,aab2[3][3];
  double At10,At11,At12;
  double At20,At21,At22;
  double Bt1,Bt2,Bt3,Bt4;
  double FcnjH2,pHH0;
  double pcc[4][4],pch[4][4],phc0;

  uint8_t nijbr[400],lstcnbbr[400][3],lstcnb[2][400][3];
  double wijbr[400],wijstor[2][400];
  double Nelijbr[400],Nelijstor[2][400];
  double dNelijbrdki[400][2],dNelijdkistor[2][400][2];

  void subFcnj_CH(int,int,int,int,double,Vec3d,double,
  vector<uint8_t>&,LNBStr*,double&,double*,FATPLoc&,NBStr&);
  void F1new(int,int,double,double,double&,double*);
  void frcsFcnj(int,int,int,int,double,double,double,vector<uint8_t>&,
  LNBStr*,FATPLoc&,NBStr&);
  void frcsFcnjFull(size_t ijp[2],int,int,double,double,double,double,
  double,double,double,LNBStr*,FATPLoc&,vector<Vec3d>&,vector<Mat3d>&);
  void torsion1(size_t,size_t,double,double,double,double,Vec3d,
  double&,double&,double&,LNBStr*,vector<Vec3d>&);
  void torsion2(size_t ijp[2],double,double,double,double,Vec3d,
  double&,double&,double&,LNBStr*,vector<Vec3d>&,vector<Mat3d>&);
  void TdTyz1(double,double,double&,double&,double&);
  void TdTyz2(double,double,double&,double&,double&);
  void AdAz(int,int,double,double,double&,double&,double&);
  void subPCH(int,int,int,int,int,double,double&,LNBStr*,
  FATPLoc&,NBStr&);
  void subPHC(int,int,double,double&,LNBStr*,NBStr&);
  void subPHH(int,double,double&,LNBStr*,NBStr&);
  void subFcnj_HH(int,double,double&,LNBStr*,FATPLoc&,NBStr&);
};

//-----------------------------------------------------------------------//
