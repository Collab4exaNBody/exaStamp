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

// 3D complex-to-complex FFT on a nx*ny*nz mesh stored with x fastest (index (iz*ny+iy)*nx+ix).
// Unnormalized, forward = exp(-i k.x), backward = exp(+i k.x), as LAMMPS FFT3d used by PPPM.
// GPU : cuFFT (Cuda builds, EXASTAMP_PPPM_CUFFT), in place on unified memory, asynchronous on the given stream :
// call sync() before the results are used from another stream or from the host.
// CPU : pocketfft (vendored in pocketfft/). Both are only included by pppm_fft.cpp.

#include <onika/math/basic_types.h>

namespace exaStamp
{
inline namespace coulombic_ewald
{
  using exanb::Complexd;

  class PPPMFFT
  {
  public:
    PPPMFFT() = default;
    PPPMFFT( const PPPMFFT& ) = delete;
    PPPMFFT& operator = ( const PPPMFFT& ) = delete;
    ~PPPMFFT();

    // use_gpu : run with cuFFT (only possible when compiled with EXASTAMP_PPPM_CUFFT), stream = cudaStream_t to run on
    void resize( int nx, int ny, int nz, bool use_gpu = false, void* stream = nullptr );
    void forward( Complexd* data ) const;
    void backward( Complexd* data ) const;
    void sync() const; // wait for transforms launched on the GPU stream (no-op on CPU)
    inline bool on_gpu() const { return m_gpu; }

    static bool gpu_support();

  private:
    void exec( Complexd* data, bool forward ) const;
    void release();

    int m_nx = 0, m_ny = 0, m_nz = 0;
    bool m_gpu = false;
    void* m_stream = nullptr;
    int m_plan = 0;          // cufftHandle
    bool m_has_plan = false;
  };

}
}
