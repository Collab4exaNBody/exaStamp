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
#include <omp.h>
#include <algorithm>
#include <cstdint>
#include <type_traits>

// pocketfft's own thread pool would compete with OpenMP
#define POCKETFFT_NO_MULTITHREADING
#include <exaStamp/potential/coulombic/pocketfft/pocketfft_hdronly.h>

#if defined(EXASTAMP_PPPM_CUFFT)
#include <cuda_runtime.h>
#include <cufft.h>
#elif defined(EXASTAMP_PPPM_HIPFFT)
#include <hip/hip_runtime_api.h>
#if __has_include(<hipfft/hipfft.h>)
#include <hipfft/hipfft.h>
#else
#include <hipfft.h>
#endif
#endif

#if defined(EXASTAMP_PPPM_CUFFT) || defined(EXASTAMP_PPPM_HIPFFT)
#define EXASTAMP_PPPM_GPUFFT 1
#endif

namespace exaStamp
{
inline namespace coulombic_ewald
{
  static_assert( sizeof(Complexd) == sizeof(std::complex<double>) );

#ifdef EXASTAMP_PPPM_GPUFFT
  // GPU FFT backend : cuFFT and hipFFT have the same API up to the prefix. Plans are stored as an opaque integer
  // (cufftHandle is an int, hipfftHandle a pointer).
  namespace gpufft
  {
#   if defined(EXASTAMP_PPPM_CUFFT)
    using Handle = cufftHandle;
    using DoubleComplex = cufftDoubleComplex;
    using Stream = cudaStream_t;
    static constexpr const char* name = "cuFFT";
    static inline bool ok( cufftResult r ) { return r == CUFFT_SUCCESS; }
    static inline auto plan_3d( Handle* p, int nz, int ny, int nx ) { return cufftPlan3d( p , nz , ny , nx , CUFFT_Z2Z ); }
    static inline auto plan_many( Handle* p, int rank, int* n, int dist, int batch ) { return cufftPlanMany( p , rank , n , nullptr , 1 , dist , nullptr , 1 , dist , CUFFT_Z2Z , batch ); }
    static inline auto set_stream( Handle p, Stream s ) { return cufftSetStream( p , s ); }
    static inline auto exec( Handle p, DoubleComplex* c, bool forward ) { return cufftExecZ2Z( p , c , c , forward ? CUFFT_FORWARD : CUFFT_INVERSE ); }
    static inline void destroy( Handle p ) { cufftDestroy( p ); }
    static inline void stream_sync( Stream s )
    {
      const cudaError_t r = cudaStreamSynchronize( s );
      if( r != cudaSuccess ) ::onika::fatal_error() << "PPPM cuda error '"<<cudaGetErrorString(r)<<"' in cudaStreamSynchronize" << std::endl;
    }
#   else
    using Handle = hipfftHandle;
    using DoubleComplex = hipfftDoubleComplex;
    using Stream = hipStream_t;
    static constexpr const char* name = "hipFFT";
    static inline bool ok( hipfftResult r ) { return r == HIPFFT_SUCCESS; }
    static inline auto plan_3d( Handle* p, int nz, int ny, int nx ) { return hipfftPlan3d( p , nz , ny , nx , HIPFFT_Z2Z ); }
    static inline auto plan_many( Handle* p, int rank, int* n, int dist, int batch ) { return hipfftPlanMany( p , rank , n , nullptr , 1 , dist , nullptr , 1 , dist , HIPFFT_Z2Z , batch ); }
    static inline auto set_stream( Handle p, Stream s ) { return hipfftSetStream( p , s ); }
    static inline auto exec( Handle p, DoubleComplex* c, bool forward ) { return hipfftExecZ2Z( p , c , c , forward ? HIPFFT_FORWARD : HIPFFT_BACKWARD ); }
    static inline void destroy( Handle p ) { hipfftDestroy( p ); }
    static inline void stream_sync( Stream s )
    {
      const hipError_t r = hipStreamSynchronize( s );
      if( r != hipSuccess ) ::onika::fatal_error() << "PPPM hip error '"<<hipGetErrorString(r)<<"' in hipStreamSynchronize" << std::endl;
    }
#   endif

    static_assert( sizeof(Complexd) == sizeof(DoubleComplex) );
    static_assert( sizeof(Handle) <= sizeof(std::intptr_t) );

    template<class R> static inline void check( R r, const char* what )
    {
      if( ! ok(r) ) ::onika::fatal_error() << "PPPM "<<name<<" error "<<int(r)<<" in "<<what << std::endl;
    }
    // templates, so that only the cast matching the handle type is instantiated
    template<class H = Handle> static inline std::intptr_t store( H p )
    {
      if constexpr ( std::is_pointer_v<H> ) return reinterpret_cast<std::intptr_t>( p );
      else return static_cast<std::intptr_t>( p );
    }
    template<class H = Handle> static inline H load( std::intptr_t p )
    {
      if constexpr ( std::is_pointer_v<H> ) return reinterpret_cast<H>( p );
      else return static_cast<H>( p );
    }
    static inline std::intptr_t make_plan_3d( int nz, int ny, int nx, void* stream )
    {
      Handle p;
      check( plan_3d( &p , nz , ny , nx ) , "plan 3d" );
      check( set_stream( p , static_cast<Stream>(stream) ) , "set stream" );
      return store( p );
    }
    static inline std::intptr_t make_plan_many( int rank, int* n, int dist, int batch, void* stream, const char* what )
    {
      Handle p;
      check( plan_many( &p , rank , n , dist , batch ) , what );
      check( set_stream( p , static_cast<Stream>(stream) ) , "set stream" );
      return store( p );
    }
    static inline void run( std::intptr_t plan, Complexd* data, bool forward, const char* what )
    {
      check( exec( load(plan) , reinterpret_cast<DoubleComplex*>(data) , forward ) , what );
    }
  }
#endif

  bool PPPMFFT::gpu_support()
  {
#ifdef EXASTAMP_PPPM_GPUFFT
    return true;
#else
    return false;
#endif
  }

  void PPPMFFT::release()
  {
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_has_plan ) gpufft::destroy( gpufft::load(m_plan) );
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
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_gpu )
    {
      // row major nz x ny x nx : x is the fastest index, as the mesh storage
      m_plan = gpufft::make_plan_3d( nz , ny , nx , stream );
      m_has_plan = true;
    }
#endif
  }

  void PPPMFFT::exec( Complexd* data, bool forward ) const
  {
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_gpu )
    {
      gpufft::run( m_plan , data , forward , "exec 3d" );
      return;
    }
#endif
    const ptrdiff_t cs = sizeof(std::complex<double>);
    auto* c = reinterpret_cast<std::complex<double>*>( data );
    if( omp_get_max_threads() <= 1 )
    {
      const pocketfft::shape_t shape = { size_t(m_nz) , size_t(m_ny) , size_t(m_nx) };
      const pocketfft::stride_t stride = { ptrdiff_t(m_ny)*m_nx*cs , ptrdiff_t(m_nx)*cs , cs };
      pocketfft::c2c( shape , stride , stride , { 0 , 1 , 2 } , forward , c , c , 1.0 , 1 );
      return;
    }
    // OpenMP : 2D transforms of the xy planes, then 1D transforms along z, one xz plane per iteration
    const size_t nxy = size_t(m_ny) * m_nx;
    const pocketfft::shape_t plane_shape = { size_t(m_ny) , size_t(m_nx) };
    const pocketfft::stride_t plane_stride = { ptrdiff_t(m_nx)*cs , cs };
#   pragma omp parallel for schedule(static)
    for( int iz = 0 ; iz < m_nz ; iz++ )
    {
      pocketfft::c2c( plane_shape , plane_stride , plane_stride , { 0 , 1 } , forward , c + iz*nxy , c + iz*nxy , 1.0 , 1 );
    }
    const pocketfft::shape_t zx_shape = { size_t(m_nz) , size_t(m_nx) };
    const pocketfft::stride_t zx_stride = { ptrdiff_t(nxy)*cs , cs };
#   pragma omp parallel for schedule(static)
    for( int iy = 0 ; iy < m_ny ; iy++ )
    {
      pocketfft::c2c( zx_shape , zx_stride , zx_stride , { 0 } , forward , c + size_t(iy)*m_nx , c + size_t(iy)*m_nx , 1.0 , 1 );
    }
  }

  // ---------------------------------- distributed FFT building blocks ----------------------------------

  void PPPMDistFFT::release()
  {
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_has_plan_planes ) gpufft::destroy( gpufft::load(m_plan_planes) );
    if( m_has_plan_columns ) gpufft::destroy( gpufft::load(m_plan_columns) );
#endif
    m_has_plan_planes = m_has_plan_columns = false;
  }

  PPPMDistFFT::~PPPMDistFFT()
  {
    release();
  }

  void PPPMDistFFT::resize( int nx, int ny, int nz, int nzl, int ncol, bool use_gpu, void* stream )
  {
    use_gpu = use_gpu && PPPMFFT::gpu_support();
    if( nx == m_nx && ny == m_ny && nz == m_nz && nzl == m_nzl && ncol == m_ncol && use_gpu == m_gpu && stream == m_stream ) return;
    release();
    m_nx = nx; m_ny = ny; m_nz = nz; m_nzl = nzl; m_ncol = ncol;
    m_gpu = use_gpu;
    m_stream = stream;
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_gpu )
    {
      if( nzl > 0 )
      {
        int n[2] = { ny , nx };
        m_plan_planes = gpufft::make_plan_many( 2 , n , nx*ny , nzl , stream , "plan planes" );
        m_has_plan_planes = true;
      }
      if( ncol > 0 )
      {
        int n[1] = { nz };
        m_plan_columns = gpufft::make_plan_many( 1 , n , nz , ncol , stream , "plan columns" );
        m_has_plan_columns = true;
      }
    }
#endif
  }

  void PPPMDistFFT::planes( Complexd* data, bool forward ) const
  {
    if( m_nzl == 0 ) return;
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_gpu )
    {
      gpufft::run( m_plan_planes , data , forward , "exec planes" );
      return;
    }
#endif
    const ptrdiff_t cs = sizeof(std::complex<double>);
    auto* c = reinterpret_cast<std::complex<double>*>( data );
    const size_t nxy = size_t(m_ny) * m_nx;
    const pocketfft::shape_t plane_shape = { size_t(m_ny) , size_t(m_nx) };
    const pocketfft::stride_t plane_stride = { ptrdiff_t(m_nx)*cs , cs };
#   pragma omp parallel for schedule(static)
    for( int iz = 0 ; iz < m_nzl ; iz++ )
    {
      pocketfft::c2c( plane_shape , plane_stride , plane_stride , { 0 , 1 } , forward , c + iz*nxy , c + iz*nxy , 1.0 , 1 );
    }
  }

  void PPPMDistFFT::columns( Complexd* data, bool forward ) const
  {
    if( m_ncol == 0 ) return;
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_gpu )
    {
      gpufft::run( m_plan_columns , data , forward , "exec columns" );
      return;
    }
#endif
    const ptrdiff_t cs = sizeof(std::complex<double>);
    auto* c = reinterpret_cast<std::complex<double>*>( data );
    // chunks of columns, one pocketfft call per chunk (vectorized over columns inside pocketfft)
    const int nthreads = omp_get_max_threads();
    const int nchunks = std::min( m_ncol , std::max( 1 , 4*nthreads ) );
#   pragma omp parallel for schedule(static)
    for( int ch = 0 ; ch < nchunks ; ch++ )
    {
      const int c0 = int( ( int64_t(ch) * m_ncol ) / nchunks ), c1 = int( ( int64_t(ch+1) * m_ncol ) / nchunks );
      if( c1 <= c0 ) continue;
      const pocketfft::shape_t shape = { size_t(c1-c0) , size_t(m_nz) };
      const pocketfft::stride_t stride = { ptrdiff_t(m_nz)*cs , cs };
      pocketfft::c2c( shape , stride , stride , { 1 } , forward , c + size_t(c0)*m_nz , c + size_t(c0)*m_nz , 1.0 , 1 );
    }
  }

  void PPPMDistFFT::sync() const
  {
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_gpu ) gpufft::stream_sync( static_cast<gpufft::Stream>(m_stream) );
#endif
  }

  // ---------------------------------- replicated mesh FFT ----------------------------------

  void PPPMFFT::forward( Complexd* data ) const { exec( data , true ); }
  void PPPMFFT::backward( Complexd* data ) const { exec( data , false ); }

  void PPPMFFT::sync() const
  {
#ifdef EXASTAMP_PPPM_GPUFFT
    if( m_gpu ) gpufft::stream_sync( static_cast<gpufft::Stream>(m_stream) );
#endif
  }
}
}
