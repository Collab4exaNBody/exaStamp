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

#include <exaStamp/potential/snaplegacy/SnapLegacyCG.h>
#include <iostream>
#include <iomanip> // setw.
#include <cmath> // sqrt.
            
using namespace std;                
                                              
//Constructor
SnapLegacyCG::SnapLegacyCG(double const &jmax, int const &nt) :  m_cg_tab(0.), m_tab_idxj(0), m_nt(nt), m_nt_jmax(floor(nt*jmax)), m_jmax(jmax) {}; 

//Member functions
//
//Public member functions 
//
//Returns cg coefficient for (j,j1,j2,m,m1,m2)
//Computes all the cg coefficients
int SnapLegacyCG::compute() {
  int cnt_cg=0;
# ifndef NDEBUG
  cout << "compute_cg: calcul des cg, m_nt_jmax=" << m_nt_jmax << endl;
  cout << "m_nt=" << m_nt << endl;
  cout << "m_jmax=" << m_jmax << endl;
  cout << "floor(m_nt*_jmax)" << floor(m_nt*m_jmax) << endl;
# endif
  m_nt_jmax=floor(m_nt*m_jmax);
  m_tab_idxj.resize(pow(m_nt_jmax+1,3));
  for (int j = 0; j <= m_nt_jmax; ++j) {
    for (int j1 = 0; j1 <= m_nt_jmax; ++j1) {
      for (int j2 = 0; j2 <= m_nt_jmax; ++j2) {
	      m_tab_idxj[idx_tabj(j,j1,j2)]=cnt_cg;
        if (j1+j2-j  < 0) {m_tab_idxj[idx_tabj(j,j1,j2)]=-99; continue;} // undefined
        if (j-j1+j2  < 0) {m_tab_idxj[idx_tabj(j,j1,j2)]=-99; continue;} // undefined
        if (j-j2+j1  < 0) {m_tab_idxj[idx_tabj(j,j1,j2)]=-99; continue;} // undefined
        if ((j+j1+j2)%m_nt==1) {m_tab_idxj[idx_tabj(j,j1,j2)]=-99; continue;}  //undefined: j+j1+j2 integer
        for (int m = -j; m <= j; m+=m_nt) {
          for (int m1 = -j1; m1 <= j1; m1+=m_nt) {
            for (int m2 = -j2; m2 <= j2; m2+=m_nt) {
              double cgval = 0.;
              if (m-m1-m2 != 0) {m_cg_tab.push_back(cgval); cnt_cg+=1; continue;} // push_back all defined values
              int rc = compute_cg(j1, m1, j2, m2, j, m, cgval);
              if (rc != 0) {cerr << "Error: compute_cg KO" << endl; return 1;}
              m_cg_tab.push_back(cgval);
	            cnt_cg+=1;
  	    }
	  }
	}
      }
    }
  }
  return 0;
}

//Private member functions
//
//Returns the index in the array of cg coefficients
//Computes one cg coefficient
int SnapLegacyCG::compute_cg(int const j1, int const m1,
                       int const j2, int const m2,
                       int const j,  int const m,
                       double & onecg) {
  onecg = 0.; // Reset out parameter.
  // http://mathworld.wolfram.com/Wigner3j-Symbol.html (19)
  double w3j = 0.;
  int rc = wigner3j(j1, m1, j2, m2, j, -m, w3j);
  if (rc != 0) {cerr << "Error: wigner3j KO" << endl; return 1;}
  onecg = pow(-1., (m+j1-j2)/m_nt)*sqrt((2.*j)/m_nt+1.)*w3j;
  return 0;
}

//Computes Wigner3j symbol for (j1,m1,j2,m2,j,m)
int SnapLegacyCG::wigner3j(int const j1, int const m1,
                     int const j2, int const m2,
                     int const j,  int const m,
                     double & w3j) {
  w3j = 0.; // Reset out parameter.

  // http://mathworld.wolfram.com/TriangleCoefficient.html

  int rc = 0;
  double tcn = 1.; // Triangle coefficient: numerator.
  rc = factorial( (j1+j2-j)/m_nt,    tcn); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  rc = factorial( (j1-j2+j)/m_nt,    tcn); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  rc = factorial((-j1+j2+j)/m_nt,    tcn); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  double tcd = 1.; // Triangle coefficient: denominator.
  rc = factorial( (j1+j2+j)/m_nt +1, tcd); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  double tc = tcn/tcd; // Triangle coefficient.

  // http://mathworld.wolfram.com/Wigner3j-Symbol.html (7)

  double coef = 1.;
  rc = factorial((j1+m1)/m_nt, coef); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  rc = factorial((j1-m1)/m_nt, coef); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  rc = factorial((j2+m2)/m_nt, coef); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  rc = factorial((j2-m2)/m_nt, coef); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  rc = factorial((j +m )/m_nt, coef); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
  rc = factorial((j -m )/m_nt, coef); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}

  double sum = 0.;
  int tmin = max((j2-j -m1)/m_nt, max((j1+m2-j)/m_nt,     0      ));
  int tmax = min((j1+j2- j)/m_nt, min((j1-m1)/m_nt  , (j2+m2)/m_nt));
  for (int t = tmin; t <= tmax; ++t) {
    double x = 1.;
    rc = factorial(               +t, x); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
    rc = factorial((-j2+j+m1)/m_nt +t, x); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
    rc = factorial((-j1+j-m2)/m_nt +t, x); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
    rc = factorial(( j1+j2-j)/m_nt -t, x); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
    rc = factorial(( j1-m1  )/m_nt -t, x); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
    rc = factorial(( j2+m2  )/m_nt -t, x); if (rc != 0) {cerr << "Error: bad factorial" << endl; return 1;}
//    cout << "w3j : x:" <<  "\t" << x << endl;
    sum += pow(-1., t)/x;
  }

  w3j = pow(-1., (j1-j2-m)/m_nt)*sqrt(tc*coef)*sum;

  return 0;
}

//Returns f*n!
int SnapLegacyCG::factorial(int const n, double & f) {
  if (n > 168) {cerr << "Error: factorial > 168 KO" << endl; return 1;}
  f *= m_factorial[n]; // Tabulated factorial (for performance).
  return 0;
}

//Tabulated n!
double const SnapLegacyCG::m_factorial[168] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  1.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,
  8.15915283247898e+47,
  3.34525266131638e+49,
  1.40500611775288e+51,
  6.04152630633738e+52,
  2.65827157478845e+54,
  1.1962222086548e+56,
  5.50262215981209e+57,
  2.58623241511168e+59,
  1.24139155925361e+61,
  6.08281864034268e+62,
  3.04140932017134e+64,
  1.55111875328738e+66,
  8.06581751709439e+67,
  4.27488328406003e+69,
  2.30843697339241e+71,
  1.26964033536583e+73,
  7.10998587804863e+74,
  4.05269195048772e+76,
  2.35056133128288e+78,
  1.3868311854569e+80,
  8.32098711274139e+81,
  5.07580213877225e+83,
  3.14699732603879e+85,
  1.98260831540444e+87,
  1.26886932185884e+89,
  8.24765059208247e+90,
  5.44344939077443e+92,
  3.64711109181887e+94,
  2.48003554243683e+96,
  1.71122452428141e+98,
  1.19785716699699e+100,
  8.50478588567862e+101,
  6.12344583768861e+103,
  4.47011546151268e+105,
  3.30788544151939e+107,
  2.48091408113954e+109,
  1.88549470166605e+111,
  1.45183092028286e+113,
  1.13242811782063e+115,
  8.94618213078297e+116,
  7.15694570462638e+118,
  5.79712602074737e+120,
  4.75364333701284e+122,
  3.94552396972066e+124,
  3.31424013456535e+126,
  2.81710411438055e+128,
  2.42270953836727e+130,
  2.10775729837953e+132,
  1.85482642257398e+134,
  1.65079551609085e+136,
  1.48571596448176e+138,
  1.3520015276784e+140,
  1.24384140546413e+142,
  1.15677250708164e+144,
  1.08736615665674e+146,
  1.03299784882391e+148,
  9.91677934870949e+149,
  9.61927596824821e+151,
  9.42689044888324e+153,
  9.33262154439441e+155,
  9.33262154439441e+157,
  9.42594775983835e+159,
  9.61446671503512e+161,
  9.90290071648618e+163,
  1.02990167451456e+166,
  1.08139675824029e+168,
  1.14628056373471e+170,
  1.22652020319614e+172,
  1.32464181945183e+174,
  1.44385958320249e+176,
  1.58824554152274e+178,
  1.76295255109024e+180,
  1.97450685722107e+182,
  2.23119274865981e+184,
  2.54355973347219e+186,
  2.92509369349301e+188,
  3.3931086844519e+190,
  3.96993716080872e+192,
  4.68452584975429e+194,
  5.5745857612076e+196,
  6.68950291344912e+198,
  8.09429852527344e+200,
  9.8750442008336e+202,
  1.21463043670253e+205,
  1.50614174151114e+207,
  1.88267717688893e+209,
  2.37217324288005e+211,
  3.01266001845766e+213,
  3.8562048236258e+215,
  4.97450422247729e+217,
  6.46685548922047e+219,
  8.47158069087882e+221,
  1.118248651196e+224,
  1.48727070609069e+226,
  1.99294274616152e+228,
  2.69047270731805e+230,
  3.65904288195255e+232,
  5.01288874827499e+234,
  6.91778647261949e+236,
  9.61572319694109e+238,
  1.34620124757175e+241,
  1.89814375907617e+243,
  2.69536413788816e+245,
  3.85437071718007e+247,
  5.5502938327393e+249,
  8.04792605747199e+251,
  1.17499720439091e+254,
  1.72724589045464e+256,
  2.55632391787286e+258,
  3.80892263763057e+260,
  5.71338395644585e+262,
  8.62720977423323e+264,
  1.31133588568345e+267,
  2.00634390509568e+269,
  3.08976961384735e+271,
  4.78914290146339e+273,
  7.47106292628289e+275,
  1.17295687942641e+278,
  1.85327186949373e+280,
  2.94670227249504e+282,
  4.71472363599206e+284,
  7.59070505394721e+286,
  1.22969421873945e+289,
  2.0044015765453e+291,
  3.28721858553429e+293,
  5.42391066613159e+295,
  9.00369170577843e+297,
  1.503616514865e+300, // nmax factorial = 168
};
