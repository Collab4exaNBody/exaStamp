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

#include <exaStamp/potential/coulombic/pppm_fft.h>
#include <complex>
#include <cstddef>

// pocketfft's own thread pool would compete with OpenMP
#define POCKETFFT_NO_MULTITHREADING
#include <exaStamp/potential/coulombic/pocketfft/pocketfft_hdronly.h>

namespace exaStamp
{
inline namespace coulombic_ewald
{
  static_assert( sizeof(Complexd) == sizeof(std::complex<double>) );

  void PPPMFFT::resize( int nx, int ny, int nz )
  {
    m_nx = nx; m_ny = ny; m_nz = nz;
  }

  static inline void pppm_c2c( int nx, int ny, int nz, Complexd* data, bool forward )
  {
    const pocketfft::shape_t shape = { size_t(nz) , size_t(ny) , size_t(nx) };
    const ptrdiff_t cs = sizeof(std::complex<double>);
    const pocketfft::stride_t stride = { ptrdiff_t(ny)*nx*cs , ptrdiff_t(nx)*cs , cs };
    auto* c = reinterpret_cast<std::complex<double>*>( data );
    pocketfft::c2c( shape , stride , stride , { 0 , 1 , 2 } , forward , c , c , 1.0 , 1 );
  }

  void PPPMFFT::forward( Complexd* data ) const { pppm_c2c( m_nx , m_ny , m_nz , data , true ); }
  void PPPMFFT::backward( Complexd* data ) const { pppm_c2c( m_nx , m_ny , m_nz , data , false ); }
}
}
