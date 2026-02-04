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

#include <iostream>
#include <sstream> // stringstream.
#include <fstream> // ifstream.
#include <cmath> // fabs.

#include <exaStamp/potential/snaplegacy/SnapLegacyCG.h>

using namespace std;

int main (int argc , char ** argv) {
  if (argc != 3) {cerr << "Error: need jmax as argument" << endl; return 1;}
  stringstream ss1(argv[1]);
  double jmax = 0.; ss1 >> jmax;
  stringstream ss2(argv[2]);
  int nt = 0; ss2 >> nt;
  if (!ss1 || jmax <= 0) {cerr << "Error: bad jmax" << endl; return 1;}

  cout << "snap-compute-cg-cpp: jmax = " << jmax << endl;
  cout << "snap-compute-cg-cpp: nt = " << nt << endl;

  SnapLegacyCG cg(jmax,nt);
  int rc= cg.compute();
    cout << "snap-compute-cg-cpp: calcul cg fait. Comparaison" << endl;
  if (rc != 0) {cerr << "Error: compute_cg KO" << endl; return 1;}

  double valcg;
  valcg=0.;
  int nt_jmax=floor(nt*jmax);
  for (int j = 0; j <= nt_jmax; ++j) {
    for (int j1 = 0; j1 <= nt_jmax; ++j1) {
      for (int j2 = 0; j2 <= nt_jmax; ++j2) {
        for (int m = -j; m <= j; m=m+nt) {
          for (int m1 = -j1; m1 <= j1; m1=m1+nt) {
            for (int m2 = -j2; m2 <= j2; m2=m2+nt) {
              if (j1+j2-j  < 0)  continue; // undefined
              if (j-j1+j2  < 0)  continue; // undefined
              if (j-j2+j1  < 0)  continue; // undefined
	      if ((j+j1+j2)%nt==1) continue;
	      if (m-m1-m2!=0) continue; //nul values
	      valcg= cg.val(j,j1,j2,m,m1,m2);
	      cout << "j,j1,j2,m,m1,m2,cg" << "\t" << j << "\t " << j1 <<  "\t "<< j2 << "\t "<< m << "\t "<< m1 << "\t "<< m2 << "\t" << valcg << endl;
  	    }
	  }
	}
      }
    }
  }
  int cgsize=cg.size();
  cout << "snap-compute-cg-cpp: read cg.ref" << endl;

//  vector<double> cg_ref;
//  cg_ref.reserve(cgsize);
  ifstream ifs("cg.ref");
  while (ifs) {
    int j1 = 0, j2 = 0, j = 0, m1 = 0, m2 = 0, m = 0;
    double cgr = 0.;
    ifs >> j >> j1 >> j2 >> m >> m1 >> m2 >> cgr;
//    cout << "j,j1,j2,m,m1,m2,cgr:" << "\t" << j << "\t " << j1 <<  "\t "<< j2 << "\t "<< m << "\t "<< m1 << "\t "<< m2 << "\t" << cgr << endl;
if (ifs) cout << "j,j1,j2,m,m1,m2,cgr,cg_val:" << "\t" << j << "\t " << j1 <<  "\t "<< j2 << "\t "<< m << "\t "<< m1 << "\t "<< m2 << "\t" << cgr  << "\t " << cg.val(j,j1,j2,m,m1,m2) << endl;

//    if (ifs) cg_ref.push_back(cgr);
    if (ifs) { if(fabs(cgr-cg.val(j,j1,j2,m,m1,m2))>1.e-5) {cerr << "snap-compute-cg: cmp diff KO" << endl; return 1;}}
  }

  cout << "snap-compute-cg-cpp: compare size" << endl;

//  cout << "f90 " << cg_ref.size() << ", cpp " << cgsize << endl;
//  if (cg_ref.size() != cgsize) {cerr << "snap-compute-cg: cmp size KO" << endl; return 1;}

  cout << "snap-compute-cg-cpp: compare value" << endl;

//  for (size_t i = 0; i < cg_ref.size(); ++i) {
//   cout << "f90 " << cg_ref[i] << ", cpp " << cg[i] << endl;
//   if (fabs(cg_ref[i] - cg.cg_val(j,j1,j2,m,m1,m2) > 1.e-6) {cerr << "snap-compute-cg: cmp diff KO" << endl; return 1;}
//  }

  cout << "snap-compute-cg-cpp: OK" << endl;

  return 0;
}
