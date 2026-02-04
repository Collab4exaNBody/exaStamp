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

//#include "SnapLegacyCG.h"
//#include "SnapLegacyGSH.h"
#include <exaStamp/potential/snaplegacy/SnapLegacyBS.h>

using namespace std;

int main (int argc , char ** argv) {
  int rc=0;
  int ntype=2;
  double jmax=3.5;
  double rcut=5.;
  SnapLegacyCG cg(jmax,ntype);
  cg.compute();
  cout << "Reading input datas..."<< endl;
  double rx_tmp[1000];
  double ry_tmp[1000];
  double rz_tmp[1000];
  double factor_tmp[1000];
  ifstream ifs("atom_positions.ref");
  int i=0;
  while(ifs) {
  cout << i << endl;
  ifs >> rx_tmp[i] >> ry_tmp[i] >> rz_tmp[i] >> factor_tmp[i];
  cout << "rx,ry,rz:" <<  rx_tmp[i] << "\t" << ry_tmp[i] << "\t" << rz_tmp[i] << endl;
  if (ifs) {i+=1;}
//  i+=1;
  }
  size_t N_atom=i;
  double rx[N_atom];
  double ry[N_atom];
  double rz[N_atom];
  for (size_t j=0; j<N_atom; j+=1) {
    rx[j]=rx_tmp[j];
    ry[j]=ry_tmp[j];
    rz[j]=rz_tmp[j];
    cout << "rx,ry,rz:" <<  rx[j] << "\t" << ry[j] << "\t" << rz[j] << endl;
  }
  int myspecy=0;
  int n_species=1;
  int bs_size=26;
  int species[N_atom];
  double coefs[bs_size];
  double factor[n_species];
  factor[0]=1.;
  for (int i=0; i<n_species; i++) {coefs[i]=1.;}
  for (int i=0; i<N_atom; i++) {species[i]=0;}
//  SnapLegacyBS bs(jmax,rcut,myspecy);
  SnapLegacyBS bs(jmax /*,rcut*/ /*,myspecy*/ ,&coefs[0],&factor[0]);
  cout << "Refresh neighbours"<< endl;
  rc=bs.set_neighbours(&rx[0],&ry[0],&rz[0],&species[0],rcut,N_atom);
  
  
//  rc=bs.set_ml_coefs(&coefs[0],&factor[0],n_species);
  if (rc != 0) {cerr << "Error: refresh neighbours KO" << endl; return 1;}
  cout << "Compute CMM for N_neigh=" << bs.get_n_atom() << endl;
  rc=bs.compute_cmm(rcut);
//  rc = sp.compute_gsh(ran, rcut, gsh, dgsh);
  if (rc != 0) {cerr << "Error: compute cmm KO" << endl; return 1;}
  cout << "Compute BS"<< endl;
  rc=bs.compute_bs(myspecy,rcut,cg);
  if (rc != 0) {cerr << "Error: compute bs KO" << endl; return 1;}
  cout << "bs_size: " << bs_size << endl;
  cout << "Bispectrum coefficients: " << endl;
  for (i=0; i<bs_size; i++) {
   cout << bs.bs_val(i) << endl;
  }
  i=0;
  vector<double> bs_ref_real;
  bs_ref_real.resize(bs_size);
  ifstream ifs2("bs.ref2"); //Ajout de division par N_neigh dans SnapLegacyBS
  while (ifs2) {
  ifs2 >> bs_ref_real[i];
  if (ifs2) {i+=1;}
  }
  vector<double3d> dbs_ref_real;
  dbs_ref_real.resize(bs_size*N_atom);
  ifstream ifs3("dbs.ref");
  i=0;
  while (ifs3) {
  ifs3 >> dbs_ref_real[i].x >> dbs_ref_real[i].y >> dbs_ref_real[i].z ;
  if (ifs3) {i+=1;}
  }
   
  double diff_bs_real;
  double crit=1.e-7;
  for (i=0; i<bs_size; i++) {
    diff_bs_real=abs((bs_ref_real[i]-real(bs.bs_val(i)))/bs_ref_real[i]);
    if (diff_bs_real > crit) {
	cout << "PB: diff>crit for " << "\t" << i << "\t" << "th element" << endl; rc=1;
        cout << "bs_ref_real[i], real(bs.bs_val(i))" << "\t" << bs_ref_real[i] << "\t" << real(bs.bs_val(i)) << endl;
        cerr << "Error: compute_bs KO" << endl; return 1;
    }
  }
 cout << "Compute BS OK" << endl;

 double3d diff_dbs_real;
 for (int j=0; j<N_atom; j++) {
   for (i=0; i<bs_size; i++) {
     int k=i+j*bs_size;
     diff_dbs_real.x=dbs_ref_real[k].x-real(bs.dbs_val(k).x);
     diff_dbs_real.y=dbs_ref_real[k].y-real(bs.dbs_val(k).y);
     diff_dbs_real.z=dbs_ref_real[k].z-real(bs.dbs_val(k).z);
     if (diff_dbs_real.x > crit) { cout << "PB: diff>crit for dbs.x with k=" << "\t" << k << " ; bs.dbs_val(k).x=" << "\t" << bs.dbs_val(k).x << "\t" << "dbs_ref_real[k].x=" << "\t" << dbs_ref_real[k].x << endl; cerr << "Error: compute_dbs KO" << endl; return 1;}
     if (diff_dbs_real.y > crit) { cout << "PB: diff>crit for dbs.y with k=" << "\t" << k << " ; bs.dbs_val(k).y=" << "\t" << bs.dbs_val(k).y << "\t" << "dbs_ref_real[k].y=" << "\t" << dbs_ref_real[k].y << endl; cerr << "Error: compute_dbs KO" << endl; return 1;}
     if (diff_dbs_real.z > crit) { cout << "PB: diff>crit for dbs.z with k=" << "\t" << k << " ; bs.dbs_val(k).z=" << "\t" << bs.dbs_val(k).z << "\t" << "dbs_ref_real[k].z=" << "\t" << dbs_ref_real[k].z << endl; cerr << "Error: compute_dbs KO" << endl; return 1;}
   }
 }
cout << "Compute dBS OK" << endl;

return 0;

}
