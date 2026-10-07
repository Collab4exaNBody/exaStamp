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
#include <onika/log.h>
#include <complex>
#include <cstddef>

// pocketfft's own thread pool would compete with OpenMP
#define POCKETFFT_NO_MULTITHREADING
#include <exaStamp/potential/coulombic/pocketfft/pocketfft_hdronly.h>

#ifdef EXASTAMP_PPPM_CUFFT
#include <cuda_runtime.h>
#include <cufft.h>
#endif

namespace exaStamp
{
inline namespace coulombic_ewald
{
  static_assert( sizeof(Complexd) == sizeof(std::complex<double>) );

#ifdef EXASTAMP_PPPM_CUFFT
  static_assert( sizeof(Complexd) == sizeof(cufftDoubleComplex) );
  static_assert( sizeof(cufftHandle) == sizeof(int) );

  static inline void pppm_cufft_check( cufftResult r, const char* what )
  {
    if( r != CUFFT_SUCCESS ) ::onika::fatal_error() << "PPPM cuFFT error "<<int(r)<<" in "<<what << std::endl;
  }
  static inline void pppm_cuda_check( cudaError_t r, const char* what )
  {
    if( r != cudaSuccess ) ::onika::fatal_error() << "PPPM cuda error '"<<cudaGetErrorString(r)<<"' in "<<what << std::endl;
  }
#endif

  bool PPPMFFT::gpu_support()
  {
#ifdef EXASTAMP_PPPM_CUFFT
    return true;
#else
    return false;
#endif
  }

  void PPPMFFT::release()
  {
#ifdef EXASTAMP_PPPM_CUFFT
    if( m_has_plan ) cufftDestroy( static_cast<cufftHandle>(m_plan) );
#endif
    m_has_plan = false;
  }

  PPPMFFT::~PPPMFFT()
  {
    release();
  }

  void PPPMFFT::resize( int nx, int ny, int nz, bool use_gpu, void* stream )
  {
    use_gpu = use_gpu && gpu_support();
    if( nx == m_nx && ny == m_ny && nz == m_nz && use_gpu == m_gpu && stream == m_stream ) return;
    release();
    m_nx = nx; m_ny = ny; m_nz = nz;
    m_gpu = use_gpu;
    m_stream = stream;
#ifdef EXASTAMP_PPPM_CUFFT
    if( m_gpu )
    {
      cufftHandle plan;
      // row major nz x ny x nx : x is the fastest index, as the mesh storage
      pppm_cufft_check( cufftPlan3d( &plan , nz , ny , nx , CUFFT_Z2Z ) , "cufftPlan3d" );
      pppm_cufft_check( cufftSetStream( plan , static_cast<cudaStream_t>(stream) ) , "cufftSetStream" );
      m_plan = plan;
      m_has_plan = true;
    }
#endif
  }

  void PPPMFFT::exec( Complexd* data, bool forward ) const
  {
#ifdef EXASTAMP_PPPM_CUFFT
    if( m_gpu )
    {
      auto* c = reinterpret_cast<cufftDoubleComplex*>( data );
      pppm_cufft_check( cufftExecZ2Z( static_cast<cufftHandle>(m_plan) , c , c , forward ? CUFFT_FORWARD : CUFFT_INVERSE ) , "cufftExecZ2Z" );
      pppm_cuda_check( cudaStreamSynchronize( static_cast<cudaStream_t>(m_stream) ) , "cudaStreamSynchronize" );
      return;
    }
#endif
    const pocketfft::shape_t shape = { size_t(m_nz) , size_t(m_ny) , size_t(m_nx) };
    const ptrdiff_t cs = sizeof(std::complex<double>);
    const pocketfft::stride_t stride = { ptrdiff_t(m_ny)*m_nx*cs , ptrdiff_t(m_nx)*cs , cs };
    auto* c = reinterpret_cast<std::complex<double>*>( data );
    pocketfft::c2c( shape , stride , stride , { 0 , 1 , 2 } , forward , c , c , 1.0 , 1 );
  }

  void PPPMFFT::forward( Complexd* data ) const { exec( data , true ); }
  void PPPMFFT::backward( Complexd* data ) const { exec( data , false ); }
}
}
