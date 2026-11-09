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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>

#include <exaStamp/potential/coulombic/ewald.h>
#include <exaStamp/unit_system.h>

#include <exanb/core/config.h>

#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_context.h>
#include <exanb/core/xform.h>
#include <onika/flat_tuple.h>
#include <onika/parallel/parallel_for.h>
#include <mpi.h>

namespace exaStamp
{
inline namespace coulombic_ewald
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  template<bool PerAtomCharge, class ChargeOrTypeT>
  ONIKA_HOST_DEVICE_FUNC static inline double ewald_particle_charge( const ParticleSpecie* __restrict__ species, ChargeOrTypeT ct )
  {
    if constexpr ( PerAtomCharge ) return ct;
    else return species[ct].m_charge;
  }

  ONIKA_HOST_DEVICE_FUNC static inline Complexd ewald_cmul( const Complexd& a , const Complexd& b )
  {
    return { a.r * b.r - a.i * b.i , a.r * b.i + a.i * b.r };
  }

  // Reciprocal space evaluation follows LAMMPS Ewald : exp(i G.r) = exp(2.i.pi n.u) with u the fractional coordinates,
  // is the product of per direction factors exp(2.i.pi nx.ux) exp(2.i.pi ny.uy) exp(2.i.pi nz.uz), tabulated once per
  // particle for 0 <= n <= kdmax (negative n are complex conjugates).
  // Local (non ghost) particles get a compact index a = m_offset[cell] + index in cell. Tables are stored entry major :
  // entry j of particle a is at j * n + a, so that consecutive particles read consecutive addresses.
  struct EwaldTable
  {
    Complexd * __restrict__ t = nullptr;
    size_t n = 0;     // number of local particles (leading dimension)
    int oy = 0;       // entry offset of the y factors
    int oz = 0;       // entry offset of the z factors
    int kmax[3] = { 0, 0, 0 };

    static inline int entries( const ReadOnlyEwaldParameters& p ) { return p.kxmax + p.kymax + p.kzmax + 3; }

    ONIKA_HOST_DEVICE_FUNC inline Complexd& at( int j , size_t a ) const { return t[ j * n + a ]; }

    // x and y factors of exp(i G.r_a) (nx >= 0)
    ONIKA_HOST_DEVICE_FUNC inline Complexd eikr_xy( size_t a , int nx , int ny ) const
    {
      const Complexd cx = at( nx , a );
      Complexd cy = at( oy + ( ny >= 0 ? ny : -ny ) , a );
      if( ny < 0 ) cy.i = - cy.i;
      return ewald_cmul( cx , cy );
    }
    ONIKA_HOST_DEVICE_FUNC inline Complexd eikr_z( size_t a , int nz ) const
    {
      Complexd cz = at( oz + ( nz >= 0 ? nz : -nz ) , a );
      if( nz < 0 ) cz.i = - cz.i;
      return cz;
    }
    // exp(i G.r_a) for the k vector g (nx >= 0)
    ONIKA_HOST_DEVICE_FUNC inline Complexd eikr( size_t a , const EwaldCoeffs& g ) const
    {
      return ewald_cmul( eikr_xy( a , g.nx , g.ny ) , eikr_z( a , g.nz ) );
    }
  };

  // 1. per particle : charge and table of exp(2.i.pi n.u_d), u = (r - bmin) / bounds_size (grid coordinates).
  // real positions are xform.r, the cell is H = xform.diag(bounds_size), so that G.(xform.r) = 2.pi n.u + constant phase.
  template<bool PerAtomCharge>
  struct EwaldTableFunc
  {
    const ParticleSpecie * __restrict__ m_species = nullptr;
    const size_t * __restrict__ m_offset = nullptr;
    EwaldTable m_table;
    Vec3d m_bmin;
    Vec3d m_inv_size;
    double * __restrict__ m_q = nullptr;

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( size_t cell, unsigned int p, double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      const size_t a = m_offset[cell] + p;
      m_q[a] = ewald_particle_charge<PerAtomCharge>( m_species , ct );
      const double u[3] = { ( rx - m_bmin.x ) * m_inv_size.x , ( ry - m_bmin.y ) * m_inv_size.y , ( rz - m_bmin.z ) * m_inv_size.z };
      const int off[3] = { 0 , m_table.oy , m_table.oz };
      for(int d=0;d<3;d++)
      {
        double s,c;
        sincos( 2.0 * M_PI * u[d] , &s , &c );
        const Complexd e1 = { c , s };
        Complexd e = { 1.0 , 0.0 };
        m_table.at( off[d] , a ) = e;
        for(int n=1;n<=m_table.kmax[d];n++)
        {
          e = ewald_cmul( e , e1 );
          m_table.at( off[d] + n , a ) = e;
        }
      }
    }
  };

  // 2. structure factor S(k) = sum_a q_a exp(i G.r_a). Work item = (particle tile, group of k vectors sharing (nx,ny)) :
  // the x.y factor is computed once per particle and group, tile sums are accumulated locally, then added to S(k) with
  // one atomic per k vector and item (consecutive items share the tile).
  static constexpr size_t EWALD_RHO_CPU_TILE = 128;
  static constexpr size_t EWALD_RHO_GPU_TILE = 32;
  struct EwaldRhoFunc
  {
    const EwaldKGroup * __restrict__ m_groups = nullptr;
    const EwaldCoeffs * __restrict__ m_gdata = nullptr;
    EwaldTable m_table;
    size_t m_ngroups = 0;
    size_t m_tile = 0;
    const double * __restrict__ m_q = nullptr;
    Complexd * __restrict__ m_rho = nullptr;

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t i ) const
    {
      const size_t tile = i / m_ngroups;
      const EwaldKGroup grp = m_groups[ i % m_ngroups ];
      const EwaldCoeffs& g0 = m_gdata[ grp.k0 ];
      const int nx = g0.nx, ny = g0.ny, nz0 = g0.nz;
      const unsigned int cnt = grp.count;
      const size_t a0 = tile * m_tile;
      const size_t a1 = ( a0 + m_tile < m_table.n ) ? a0 + m_tile : m_table.n;
      Complexd acc[EWALD_MAX_KGROUP];
#     ifdef ONIKA_GPU_DEVICE_COMPILE
      // GPU : one particle at a time, x.y factor kept in registers
      for(unsigned int j=0;j<cnt;j++) acc[j] = Complexd{ 0.0 , 0.0 };
      for(size_t a=a0;a<a1;a++)
      {
        Complexd exy = m_table.eikr_xy( a , nx , ny );
        exy.r *= m_q[a];
        exy.i *= m_q[a];
        for(unsigned int j=0;j<cnt;j++)
        {
          const Complexd e = ewald_cmul( exy , m_table.eikr_z( a , nz0 + int(j) ) );
          acc[j].r += e.r;
          acc[j].i += e.i;
        }
      }
#     else
      // CPU : q.x.y factors of the tile first, then contiguous (vectorizable) loops over particles for each k
      static_assert( EWALD_RHO_CPU_TILE <= 128 );
      double xr[128], xi[128];
      const size_t na = a1 - a0;
      for(size_t b=0;b<na;b++)
      {
        const Complexd exy = m_table.eikr_xy( a0+b , nx , ny );
        xr[b] = exy.r * m_q[a0+b];
        xi[b] = exy.i * m_q[a0+b];
      }
      for(unsigned int j=0;j<cnt;j++)
      {
        const int nz = nz0 + int(j);
        const Complexd * __restrict__ row = m_table.t + ( m_table.oz + ( nz >= 0 ? nz : -nz ) ) * m_table.n + a0;
        const double sg = nz >= 0 ? 1.0 : -1.0;
        double sr = 0.0, si = 0.0;
#       pragma omp simd reduction(+:sr,si)
        for(size_t b=0;b<na;b++)
        {
          const double cr = row[b].r , ci = sg * row[b].i;
          sr += xr[b] * cr - xi[b] * ci;
          si += xr[b] * ci + xi[b] * cr;
        }
        acc[j] = Complexd{ sr , si };
      }
#     endif
      for(unsigned int j=0;j<cnt;j++)
      {
        ONIKA_CU_ATOMIC_ADD( m_rho[grp.k0+j].r , acc[j].r );
        ONIKA_CU_ATOMIC_ADD( m_rho[grp.k0+j].i , acc[j].i );
      }
    }
  };

  template<class T>
  struct EwaldZeroFunc
  {
    T * __restrict__ m_data = nullptr;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t i ) const { m_data[i] = T{}; }
  };

  // per particle results in compact arrays, added to particle fields by EwaldScatterFunc
  struct EwaldParticleResults
  {
    double * __restrict__ f = nullptr;   // 3 per particle
    double * __restrict__ e = nullptr;   // 1 per particle
    double * __restrict__ vir = nullptr; // 6 per particle : xx yy zz xy xz yz
  };

  // per k vector data for the force pass, built once S(k) is known : A + iB = Gc S(k)
  struct EwaldKForce
  {
    double A, B;
    double Gx, Gy, Gz;
    int nx, ny, nz;
  };
  // Gv G (x) G : xx yy zz xy xz yz
  struct EwaldKVirial
  {
    double xx, yy, zz, xy, xz, yz;
  };

  struct EwaldKDataFunc
  {
    const EwaldCoeffs * __restrict__ m_gdata = nullptr;
    const Complexd * __restrict__ m_rho = nullptr;
    EwaldKForce * __restrict__ m_kf = nullptr;
    EwaldKVirial * __restrict__ m_kv = nullptr; // optional
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t k ) const
    {
      const EwaldCoeffs& g = m_gdata[k];
      m_kf[k] = EwaldKForce{ g.Gc * m_rho[k].r , g.Gc * m_rho[k].i , g.Gx , g.Gy , g.Gz , g.nx , g.ny , g.nz };
      if( m_kv != nullptr )
      {
        m_kv[k] = EwaldKVirial{ g.Gv * g.Gx * g.Gx , g.Gv * g.Gy * g.Gy , g.Gv * g.Gz * g.Gz ,
                                g.Gv * g.Gx * g.Gy , g.Gv * g.Gx * g.Gz , g.Gv * g.Gy * g.Gz };
      }
    }
  };

  // 3. reciprocal forces, and optionaly per particle energy (reciprocal + self + background) and reciprocal virial.
  // with exp(i G.r_a) = c + i.s and A + iB = Gc S(k) :
  // force                : F_a = 2 q_a sum_k ( A s - B c ) G
  // per particle energy  : e_a = q_a sum_k ( A c + B s ) , sums up to sum_k Gc |S(k)|^2
  // per particle virial  : W_a = q_a sum_k ( A c + B s ) ( I - Gv G (x) G ) (same as LAMMPS Ewald per atom virial)
  // Two work decompositions :
  //  - m_block == 0 (GPU) : work item = (k chunk, particle). With several chunks (to expose enough parallelism),
  //    partial sums are added atomically to zeroed result arrays ; with one chunk results are stored directly.
  //    k vectors are ordered with nz fastest : the x.y factor is reused while (nx,ny) does not change.
  //  - m_block > 0 (CPU) : work item = block of m_block consecutive particles, all k vectors, by groups sharing (nx,ny),
  //    with contiguous (vectorizable) loops over the particles of the block.
  static constexpr size_t EWALD_FORCE_CPU_BLOCK = 64;

  template<bool ComputeEnergy, bool ComputeVirial>
  struct EwaldLongRangeForceComputeFunc
  {
    ReadOnlyEwaldParameters p;
    EwaldTable m_table;
    size_t m_kchunk = 0;   // number of k vectors per chunk (m_block == 0)
    bool m_atomic = false; // more than one chunk (m_block == 0)
    size_t m_block = 0;    // particles per work item (CPU), 0 for one particle per work item (GPU)
    const double * __restrict__ m_q = nullptr;
    const EwaldKForce * __restrict__ m_kf = nullptr;
    const EwaldKVirial * __restrict__ m_kv = nullptr;
    EwaldParticleResults m_out;

    ONIKA_HOST_DEVICE_FUNC inline void store( double& dst , double v ) const
    {
      if( m_atomic ) { ONIKA_CU_ATOMIC_ADD( dst , v ); }
      else dst = v;
    }

    // q-less sums of particle a over the k range [k0,k1), and their conversion to results
    struct Sums
    {
      double fx = 0.0, fy = 0.0, fz = 0.0, ep = 0.0;
      double wxx = 0.0, wyy = 0.0, wzz = 0.0, wxy = 0.0, wxz = 0.0, wyz = 0.0;
    };

    ONIKA_HOST_DEVICE_FUNC inline void finalize( size_t a , const Sums& r , bool add_self ) const
    {
      const double q = m_q[a];
      const double q2 = 2.0 * q;
      store( m_out.f[3*a  ] , q2 * r.fx );
      store( m_out.f[3*a+1] , q2 * r.fy );
      store( m_out.f[3*a+2] , q2 * r.fz );
      if constexpr ( ComputeEnergy )
      {
        double e = q * r.ep;
        if( add_self ) e += ewald_self_energy( p , q );
        store( m_out.e[a] , e );
        if constexpr ( ComputeVirial )
        {
          store( m_out.vir[6*a  ] , q * ( r.ep - r.wxx ) );
          store( m_out.vir[6*a+1] , q * ( r.ep - r.wyy ) );
          store( m_out.vir[6*a+2] , q * ( r.ep - r.wzz ) );
          store( m_out.vir[6*a+3] , - q * r.wxy );
          store( m_out.vir[6*a+4] , - q * r.wxz );
          store( m_out.vir[6*a+5] , - q * r.wyz );
        }
      }
    }

    ONIKA_HOST_DEVICE_FUNC inline void one_particle( size_t i ) const
    {
      const size_t n = m_table.n;
      const size_t chunk = i / n;
      const size_t a = i % n;
      const size_t k0 = chunk * m_kchunk;
      const size_t nk = p.nknz;
      const size_t k1 = ( k0 + m_kchunk < nk ) ? k0 + m_kchunk : nk;

      Sums r;
      int px = -1, py = 0;
      Complexd exy = { 0.0 , 0.0 };
      for(size_t k=k0;k<k1;k++)
      {
        const EwaldKForce& kf = m_kf[k];
        if( kf.nx != px || kf.ny != py )
        {
          px = kf.nx; py = kf.ny;
          exy = m_table.eikr_xy( a , px , py );
        }
        const Complexd e = ewald_cmul( exy , m_table.eikr_z( a , kf.nz ) );
        const double t = kf.A * e.i - kf.B * e.r;
        r.fx += t * kf.Gx;
        r.fy += t * kf.Gy;
        r.fz += t * kf.Gz;
        if constexpr ( ComputeEnergy )
        {
          const double ek = kf.A * e.r + kf.B * e.i;
          r.ep += ek;
          if constexpr ( ComputeVirial )
          {
            const EwaldKVirial& kv = m_kv[k];
            r.wxx += ek * kv.xx;
            r.wyy += ek * kv.yy;
            r.wzz += ek * kv.zz;
            r.wxy += ek * kv.xy;
            r.wxz += ek * kv.xz;
            r.wyz += ek * kv.yz;
          }
        }
      }
      finalize( a , r , chunk == 0 );
    }

#   ifndef ONIKA_GPU_DEVICE_COMPILE
    inline void particle_block( size_t ib ) const
    {
      static constexpr size_t B = EWALD_FORCE_CPU_BLOCK;
      const size_t a0 = ib * m_block;
      const size_t na = std::min( m_block , m_table.n - a0 );
      double xr[B], xi[B];
      double fx[B], fy[B], fz[B], ep[B], wxx[B], wyy[B], wzz[B], wxy[B], wxz[B], wyz[B];
      for(size_t b=0;b<na;b++)
      {
        fx[b] = fy[b] = fz[b] = ep[b] = 0.0;
        wxx[b] = wyy[b] = wzz[b] = wxy[b] = wxz[b] = wyz[b] = 0.0;
      }
      for(size_t ig=0;ig<p.nkgroups;ig++)
      {
        const EwaldKGroup grp = p.kgroups[ig];
        const EwaldKForce& kf0 = m_kf[grp.k0];
        for(size_t b=0;b<na;b++)
        {
          const Complexd exy = m_table.eikr_xy( a0+b , kf0.nx , kf0.ny );
          xr[b] = exy.r; xi[b] = exy.i;
        }
        for(unsigned int j=0;j<grp.count;j++)
        {
          const EwaldKForce kf = m_kf[grp.k0+j];
          const int nz = kf.nz;
          const Complexd * __restrict__ row = m_table.t + ( m_table.oz + ( nz >= 0 ? nz : -nz ) ) * m_table.n + a0;
          const double sg = nz >= 0 ? 1.0 : -1.0;
          EwaldKVirial kv = {};
          if constexpr ( ComputeVirial ) kv = m_kv[grp.k0+j];
#         pragma omp simd
          for(size_t b=0;b<na;b++)
          {
            const double cr = row[b].r , ci = sg * row[b].i;
            const double er = xr[b] * cr - xi[b] * ci;
            const double ei = xr[b] * ci + xi[b] * cr;
            const double t = kf.A * ei - kf.B * er;
            fx[b] += t * kf.Gx;
            fy[b] += t * kf.Gy;
            fz[b] += t * kf.Gz;
            if constexpr ( ComputeEnergy )
            {
              const double ek = kf.A * er + kf.B * ei;
              ep[b] += ek;
              if constexpr ( ComputeVirial )
              {
                wxx[b] += ek * kv.xx;
                wyy[b] += ek * kv.yy;
                wzz[b] += ek * kv.zz;
                wxy[b] += ek * kv.xy;
                wxz[b] += ek * kv.xz;
                wyz[b] += ek * kv.yz;
              }
            }
          }
        }
      }
      for(size_t b=0;b<na;b++)
      {
        const Sums r = { fx[b], fy[b], fz[b], ep[b], wxx[b], wyy[b], wzz[b], wxy[b], wxz[b], wyz[b] };
        finalize( a0+b , r , true );
      }
    }
#   endif

    ONIKA_HOST_DEVICE_FUNC inline void operator () ( size_t i ) const
    {
#     ifdef ONIKA_GPU_DEVICE_COMPILE
      one_particle( i );
#     else
      if( m_block == 0 ) one_particle( i );
      else particle_block( i );
#     endif
    }
  };

  // 4. add compact per particle results to particle fields
  struct EwaldScatterFunc
  {
    const size_t * __restrict__ m_offset = nullptr;
    EwaldParticleResults m_res;

    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( size_t cell, unsigned int pi, double & fx, double & fy, double & fz ) const
    {
      const size_t a = m_offset[cell] + pi;
      fx += m_res.f[3*a]; fy += m_res.f[3*a+1]; fz += m_res.f[3*a+2];
    }

    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( size_t cell, unsigned int pi, double & fx, double & fy, double & fz, double & ep, Mat3d & virial ) const
    {
      const size_t a = m_offset[cell] + pi;
      fx += m_res.f[3*a]; fy += m_res.f[3*a+1]; fz += m_res.f[3*a+2];
      ep += m_res.e[a];
      const double* v = m_res.vir + 6*a;
      virial.m11 += v[0]; virial.m22 += v[1]; virial.m33 += v[2];
      virial.m12 += v[3]; virial.m21 += v[3];
      virial.m13 += v[4]; virial.m31 += v[4];
      virial.m23 += v[5]; virial.m32 += v[5];
    }
  };
  
}
}

namespace exanb
{
  template<bool PerAtomCharge> struct ComputeCellParticlesTraits< exaStamp::EwaldTableFunc<PerAtomCharge> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };

  template<> struct ComputeCellParticlesTraits< exaStamp::EwaldScatterFunc >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace onika
{
  namespace parallel
  {
    template<> struct ParallelForFunctorTraits< exaStamp::EwaldRhoFunc > { static inline constexpr bool CudaCompatible = true; };
    template<class T> struct ParallelForFunctorTraits< exaStamp::EwaldZeroFunc<T> > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::EwaldKDataFunc > { static inline constexpr bool CudaCompatible = true; };
    template<bool ComputeEnergy, bool ComputeVirial>
    struct ParallelForFunctorTraits< exaStamp::EwaldLongRangeForceComputeFunc<ComputeEnergy,ComputeVirial> > { static inline constexpr bool CudaCompatible = true; };
  }
}

namespace exaStamp
{
inline namespace coulombic_ewald
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class EwaldLongRangePC : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( EwaldParameters , ewald_config , INPUT , OPTIONAL );
    ADD_SLOT( GridT      , grid         , INPUT_OUTPUT );
    ADD_SLOT( Domain     , domain       , INPUT , REQUIRED );
    ADD_SLOT( double     , rcut_max     , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( EwaldRho   , ewald_rho    , INPUT_OUTPUT , EwaldRho{} );
    ADD_SLOT( ParticleSpecies  , species           , INPUT , REQUIRED );
    ADD_SLOT( bool       , per_atom_charge , INPUT , true , DocString{"read charges from per particle charge field instead of species charges"} );
    ADD_SLOT( MPI_Comm   , mpi          , INPUT );
    ADD_SLOT( bool       , trigger_thermo_state    , INPUT , OPTIONAL );

    // work buffers kept between time steps : compact particle numbering, charges, exp(2.i.pi n.u) tables, results
    onika::memory::CudaMMVector<size_t> m_offset;
    onika::memory::CudaMMVector<double> m_q;
    onika::memory::CudaMMVector<Complexd> m_table;
    onika::memory::CudaMMVector<double> m_f;
    onika::memory::CudaMMVector<double> m_e;
    onika::memory::CudaMMVector<double> m_vir;
    onika::memory::CudaMMVector<EwaldKForce> m_kf;
    onika::memory::CudaMMVector<EwaldKVirial> m_kv;

  public:
    // Operator execution
    inline void execute () override final
    {
      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
      {
        log_energy = *trigger_thermo_state ;
      }
      else
      {
        ldbg << "trigger_thermo_state missing " << std::endl;
      }
      
      if( ! ewald_config.has_value() )
      {
        ldbg << "ewald_config not set, skip ewal_long_range computation" << std::endl;
        return ;
      }

      *rcut_max = std::max( *rcut_max , ewald_config->radius );
      
      if( grid->number_of_cells() == 0 ) return;

      // k vectors depend on the cell, refuse to compute with stale ones
      const Mat3d cell = ewald_cell_matrix( domain->xform() , domain->bounds_size() );
      if( ! ewald_same_cell( cell , ewald_config->cell ) )
      {
        fatal_error() << "coulombic_ewald_long_range : domain cell "<<cell<<" differs from the one used to build k vectors "
                      << ewald_config->cell << ". Call coulombic_ewald_init before force computation when the cell changes." << std::endl;
      }

      const size_t nk = ewald_config->nknz;
      if (ewald_rho->rho.size() < nk) {
        ewald_rho->rho.resize(nk);
      }
      ewald_rho->nk = nk;

      const ReadOnlyEwaldParameters ro_params = *ewald_config;

      // compact numbering of local (non ghost) particles : particle p of cell c has index m_offset[c] + p
      const size_t ncells = grid->number_of_cells();
      m_offset.resize( ncells );
      size_t n = 0;
      for(size_t c=0;c<ncells;c++)
      {
        m_offset[c] = n;
        if( ! grid->is_ghost_cell(c) ) n += grid->cell_number_of_particles(c);
      }
      const int nentries = EwaldTable::entries( ro_params );
      if( m_q.size() < n ) m_q.resize( n );
      if( m_table.size() < n * nentries ) m_table.resize( n * nentries );
      if( m_f.size() < 3*n ) m_f.resize( 3*n );
      if( log_energy && m_e.size() < n ) { m_e.resize( n ); m_vir.resize( 6*n ); }
      const EwaldTable table = { m_table.data() , n , int(ro_params.kxmax) + 1 , int(ro_params.kxmax + ro_params.kymax) + 2 ,
                                 { int(ro_params.kxmax) , int(ro_params.kymax) , int(ro_params.kzmax) } };

      auto rx = grid->field_accessor( field::rx );
      auto ry = grid->field_accessor( field::ry );
      auto rz = grid->field_accessor( field::rz );
      auto fx = grid->field_accessor( field::fx );
      auto fy = grid->field_accessor( field::fy );
      auto fz = grid->field_accessor( field::fz );
      auto ep = grid->field_accessor( field::ep );
      auto virial = grid->field_accessor( field::virial );

      // 1. per particle charges and exp(2.i.pi n.u) tables
      const Vec3d bsize = domain->bounds_size();
      auto compute_tables = [&]( auto per_atom_charge_tag , auto charge_or_type )
      {
        static constexpr bool PerAtomCharge = decltype(per_atom_charge_tag)::value;
        EwaldTableFunc<PerAtomCharge> table_func = { species->data() , m_offset.data() , table ,
                                                     domain->bounds().bmin , Vec3d{ 1.0/bsize.x , 1.0/bsize.y , 1.0/bsize.z } ,
                                                     m_q.data() };
        compute_cell_particles( *grid , false , table_func , onika::make_flat_tuple(rx,ry,rz,charge_or_type) , parallel_execution_context() );
      };
      if( *per_atom_charge ) compute_tables( std::true_type{}  , grid->field_accessor( field::charge ) );
      else                   compute_tables( std::false_type{} , grid->field_accessor( field::type ) );

      // 2. structure factor of local particles, summed over all ranks
      Complexd * rho = ewald_rho->rho.data();
      onika::parallel::parallel_for( nk , EwaldZeroFunc<Complexd>{ rho } , parallel_execution_context() );
      const bool gpu_available = ( global_cuda_ctx() != nullptr ) && global_cuda_ctx()->has_devices();
      if( n > 0 )
      {
        const size_t tile = gpu_available ? EWALD_RHO_GPU_TILE : EWALD_RHO_CPU_TILE;
        const size_t ntiles = ( n + tile - 1 ) / tile;
        const size_t ngroups = ro_params.nkgroups;
        EwaldRhoFunc rho_func = { ro_params.kgroups , ro_params.Gdata , table , ngroups , tile , m_q.data() , rho };
        onika::parallel::parallel_for( ntiles * ngroups , rho_func , parallel_execution_context() );
      }
      int nprocs = 1;
      MPI_Comm_size( *mpi , &nprocs );
      static_assert( sizeof(Complexd) == 2*sizeof(double) );
      if( nprocs > 1 ) MPI_Allreduce( MPI_IN_PLACE , (double*) rho , nk*2 , MPI_DOUBLE , MPI_SUM , *mpi );

      if( n == 0 ) return;

      // per k data for the force pass
      if( m_kf.size() < nk ) m_kf.resize( nk );
      if( log_energy && m_kv.size() < nk ) m_kv.resize( nk );
      onika::parallel::parallel_for( nk , EwaldKDataFunc{ ro_params.Gdata , rho , m_kf.data() , log_energy ? m_kv.data() : nullptr } , parallel_execution_context() );

      // 3. forces, energies and virial in compact arrays.
      // GPU : one thread per (k chunk, particle), k vectors split in chunks so that there are enough threads (~128k),
      // each chunk keeps at least 64 k vectors. CPU : one work item per block of EWALD_FORCE_CPU_BLOCK particles.
      size_t nchunks = 1;
      size_t block = EWALD_FORCE_CPU_BLOCK;
      size_t nitems = ( n + block - 1 ) / block;
      if( gpu_available )
      {
        block = 0;
        nchunks = std::max( size_t(1) , std::min( ( size_t(131072) + n - 1 ) / n , nk / 64 ) );
      }
      const size_t kchunk = ( nk + nchunks - 1 ) / nchunks;
      nchunks = ( nk + kchunk - 1 ) / kchunk;
      if( gpu_available ) nitems = nchunks * n;
      const bool atomic = nchunks > 1;
      const EwaldParticleResults res = { m_f.data() , m_e.data() , m_vir.data() };
      if( atomic )
      {
        onika::parallel::parallel_for( 3*n , EwaldZeroFunc<double>{ res.f } , parallel_execution_context() );
        if( log_energy )
        {
          onika::parallel::parallel_for( n , EwaldZeroFunc<double>{ res.e } , parallel_execution_context() );
          onika::parallel::parallel_for( 6*n , EwaldZeroFunc<double>{ res.vir } , parallel_execution_context() );
        }
      }
      if( log_energy )
      {
        EwaldLongRangeForceComputeFunc<true,true> force_func = { ro_params , table , kchunk , atomic , block , m_q.data() , m_kf.data() , m_kv.data() , res };
        onika::parallel::parallel_for( nitems , force_func , parallel_execution_context() );
      }
      else
      {
        EwaldLongRangeForceComputeFunc<false,false> force_func = { ro_params , table , kchunk , atomic , block , m_q.data() , m_kf.data() , m_kv.data() , res };
        onika::parallel::parallel_for( nitems , force_func , parallel_execution_context() );
      }

      // 4. add to particle fields
      EwaldScatterFunc scatter_func = { m_offset.data() , res };
      if( log_energy ) compute_cell_particles( *grid , false , scatter_func , onika::make_flat_tuple(fx,fy,fz,ep,virial) , parallel_execution_context() );
      else             compute_cell_particles( *grid , false , scatter_func , onika::make_flat_tuple(fx,fy,fz) , parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Reciprocal space part of the Ewald summation (direct sum over k vectors built by coulombic_ewald_init).
Computes forces, and when trigger_thermo_state is true, per particle energy (reciprocal + self + neutralizing background)
and per particle reciprocal virial.
)EOF";
    }

  };
  
  template<class GridT> using EwaldLongRangePCTmpl = EwaldLongRangePC<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_ewald_long_range)
  {  
    OperatorNodeFactory::instance()->register_factory( "coulombic_ewald_long_range" , make_grid_variant_operator<EwaldLongRangePCTmpl> );
  }

}
}
