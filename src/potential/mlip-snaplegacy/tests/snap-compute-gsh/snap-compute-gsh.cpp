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
#include <iomanip>
#include <sstream> // stringstream.
#include <fstream> // ifstream.
#include <cmath> // fabs.
#include <complex>

#include <exaStamp/potential/snaplegacy/SnapLegacyGSH.h>

using namespace std;

int main (int argc , char ** argv) {
  cout.precision(12);
  if (argc != 6) {cerr << "Error: need x, y, z, rcut, jmax as argument" << endl; return 1;}
  double3d ran;
  stringstream ss1(argv[1]);
  ss1 >> ran.x;
  if (!ss1) {cerr << "Error: bad x" << endl; return 1;}
  stringstream ss2(argv[2]);
  ss2 >> ran.y;
  if (!ss2) {cerr << "Error: bad y" << endl; return 1;}
  stringstream ss3(argv[3]);
  ss3 >> ran.z;
  if (!ss3) {cerr << "Error: bad z" << endl; return 1;}
  double rcut = 0.;
  stringstream ss4(argv[4]);
  ss4 >> rcut;
  if (!ss4) {cerr << "Error: bad rcut" << endl; return 1;}
  stringstream ss5(argv[5]);
  double jmax = 0; ss5 >> jmax;
  if (!ss5 || jmax <= 0) {cerr << "Error: bad jmax" << endl; return 1;}

  cout << "snap-compute-gsh-cpp: ran = " << ran.x << ", " << ran.y << ", " << ran.z << endl;
  cout << "snap-compute-gsh-cpp: rcut = " << rcut << endl;
  cout << "snap-compute-gsh-cpp: jmax = " << jmax << endl;

  SnapLegacyGSH gsh(jmax);
  cout << "Objet gsh cree" << endl;
  int rc = gsh.compute_gsh(ran, rcut);
  if (rc != 0) {cerr << "Error: compute_gsh KO" << endl; return 1;}
  int two_jmax=floor(2*jmax);
  cout << "snap-compute-gsh-cpp: read gsh.ref" << endl;

  vector<complex<double>> gsh_ref;
  vector<complex3d> dgsh_ref;
  gsh_ref.resize(gsh.size());
  dgsh_ref.resize(3*gsh.size());
  ifstream ifs("gsh.ref");
  while (ifs) {
    int m1 = 0, m2 = 0, j = 0;
    double gshr = 0., gshi = 0.;
    ifs >> j >> m1 >> m2 >> gshr >> gshi;
    if (ifs) gsh_ref[gsh.idx(j, m1, m2)] = complex<double>(gshr, gshi);
  }

  ifstream ifs2("dgsh.ref");
  while (ifs2) {
    int m1 = 0, m2 = 0, j = 0;
    double gshr = 0., gshi = 0.;
    ifs2 >> j >> m1 >> m2 >> std::setprecision(10) >> gshr >> gshi;
    if (ifs2) dgsh_ref[gsh.idx(j, m1, m2)].x = complex<double>(gshr, gshi);
    ifs2 >> j >> m1 >> m2 >> gshr >> gshi;
    if (ifs2) dgsh_ref[gsh.idx(j, m1, m2)].y = complex<double>(gshr, gshi);
    ifs2 >> j >> m1 >> m2 >> gshr >> gshi;
    if (ifs2) dgsh_ref[gsh.idx(j, m1, m2)].z = complex<double>(gshr, gshi);
  }
  cout << "snap-compute-gsh-cpp: compare value" << endl;
  double crit=1.0e-6;

  const char* labels [6] = { "realx", "realy", "realz", "iamgx", "imagy", "imagz" };

  for (int j=0; j<=two_jmax;j++) {
    for (int m2=-j; m2<=j;m2 += 2) {
      for (int m1=-j; m1<=j; m1 += 2) {
      
      double ref[6] = {
        real(dgsh_ref[gsh.idx(j, m1, m2)].x) ,
        real(dgsh_ref[gsh.idx(j, m1, m2)].y) ,
        real(dgsh_ref[gsh.idx(j, m1, m2)].z) ,
        imag(dgsh_ref[gsh.idx(j, m1, m2)].x) ,
        imag(dgsh_ref[gsh.idx(j, m1, m2)].y) ,
        imag(dgsh_ref[gsh.idx(j, m1, m2)].z)
      };
      
      double val[6] = {
        real(gsh.dgsh_val(j, m1, m2).x) ,
        real(gsh.dgsh_val(j, m1, m2).y) , 
        real(gsh.dgsh_val(j, m1, m2).z) ,
        imag(gsh.dgsh_val(j, m1, m2).x) ,
        imag(gsh.dgsh_val(j, m1, m2).y) , 
        imag(gsh.dgsh_val(j, m1, m2).z) 
        };
 
      for(int k=0;k<6;k++)
      {
        double diff = abs( (ref[k]-val[k]) / ref[k] );
        if( diff > crit )
        {
          std::cerr << "snap-compute-gsh: cmp KO - diff_"<<labels[k]<< "\t " <<j << "\t " << m1 << "\t " << m2 << "\t " <<diff << "\t" << ref[k] << "\t" << val[k] << std::endl;
          std::abort();
        }
      }

      double diff_gsh = abs((gsh_ref[gsh.idx(j, m1, m2)]-gsh.gsh_val(j, m1, m2))/gsh_ref[gsh.idx(j, m1, m2)]);
      if (diff_gsh>crit)   {cerr << "snap-compute-gsh: cmp KO - diff_gsh   j,m1,m2" << "\t " <<j << "\t " << m1 << "\t " << m2 << "\t " <<diff_gsh   << endl; std::abort();}
      }
    }
  }

  cout << "snap-compute-gsh-cpp: OK" << endl;

  return 0;
}
