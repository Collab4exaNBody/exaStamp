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

# ifdef DEBUG_PROLOG_ADDITIONAL_CODE
  DEBUG_PROLOG_ADDITIONAL_CODE
# endif

  double _ep = 0.;
  double _fx = 0.;
  double _fy = 0.;
  double _fz = 0.;

  Mat3d _vir; // default constructor defines all elements to 0
  assert( is_zero(_vir) );

# ifndef DEBUG_PER_NBH_ADDITIONAL_CODE
# pragma omp simd reduction(+:_ep,_fx,_fy,_fz,_vir)
# endif
  for(size_t i=0;i<n;i++)
  {
    const double r = std::sqrt(tab.d2[i]);
    // optional weighting
    const auto weight = tab.nbh_data.get(i);
    double e=0.0, de=0.0;

#   if USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
    USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de , weight );
#   else
    USTAMP_POTENTIAL_COMPUTE( p, pair_params, r, e, de );
    e *= weight;
    de *= weight;
#   endif
    e -= ecut * weight;
    de /= r;
    
    // force = dE * dr->
    const double drx = tab.drx[i];
    const double dry = tab.dry[i];
    const double drz = tab.drz[i];
    const double fe_x = de * drx;
    const double fe_y = de * dry;
    const double fe_z = de * drz;

    _fx += fe_x;
    _fy += fe_y;
    _fz += fe_z;
    _ep += .5 * e;

    _vir += tensor( Vec3d{fe_x,fe_y,fe_z}, Vec3d{drx,dry,drz} ) * -0.5;

#   ifdef DEBUG_PER_NBH_ADDITIONAL_CODE
    DEBUG_PER_NBH_ADDITIONAL_CODE
#   endif
  }

# ifdef DEBUG_ADDITIONAL_CODE
  DEBUG_ADDITIONAL_CODE
# endif

  ep += _ep;
  ax += _fx ; 
  ay += _fy ; 
  az += _fz ; 
  virial += _vir;

