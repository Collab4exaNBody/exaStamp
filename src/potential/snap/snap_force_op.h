#pragma once

#include <exanb/core/physics_constants.h>
#include <exanb/core/concurent_add_contributions.h>
#include <onika/cuda/cuda.h>

#include <cmath>

// tells if we use Complex arithmetic classes or unroll all to scalar expressions
//#define SNAP_AUTOGEN_COMPLEX_MATH 1

#include "snap_compute_buffer.h"
#include "snap_compute_ui.h"
#include "snap_compute_yi.h"
#include "snap_compute_duidrj.h"
#include "snap_compute_deidrj.h"

#ifdef SNAP_AUTOGEN_COMPLEX_MATH
#include "snap_math.h"
#endif

namespace exaStamp
{
  using namespace exanb;

  template<class SnapXSForceExtStorageT>
  ONIKA_HOST_DEVICE_FUNC
  static inline void add_nbh_contrib_to_force(
                 int twojmax, int idxu_max, int jelem
               , double wj_jj, double rcutij_jj, double sinnerij_jj, double dinnerij_jj
               , double x, double y, double z, double z0, double r, double rsq
               , double const * __restrict__ rootpqarray
               , int const * __restrict__ y_jju_map, int idxu_max_alt
               , double rmin0, double rfac0, bool switch_flag, bool switch_inner_flag, bool chem_flag
               , double * __restrict__ fij
               , SnapXSForceExtStorageT& ext )
  {
//    printf("generic add_nbh_contrib_to_force\n");
    // here we use UiTot as scratch buffer to compute Uarray for jj
    snap_compute_uarray( twojmax, rootpqarray, ext.m_U_array.r(), ext.m_U_array.i(), x, y, z, z0, r );

    // reads ulist, sinner and dinner
    // writes dulist
    snap_compute_duidrj( twojmax, idxu_max, x,y,z, rcutij_jj , wj_jj
                       , ext.m_U_array.r(), ext.m_U_array.i()
                       , rootpqarray
                       , sinnerij_jj, dinnerij_jj
                       , rmin0, rfac0, switch_flag, switch_inner_flag, chem_flag                             
                       , ext.m_DU_array.r(), ext.m_DU_array.i() // OUTPUTS
                       );

    // reads dulist and ylist
    // writes final force
    snap_compute_deidrj_alt( jelem, twojmax, idxu_max
                       , y_jju_map, idxu_max_alt
                       , ext.m_DU_array.r(), ext.m_DU_array.i()
                       , ext.m_Y_array.r(), ext.m_Y_array.i()
                       , fij );
  }

  template<class SnapXSForceExtStorageT, int twojmax>
  ONIKA_HOST_DEVICE_FUNC
  static inline void add_nbh_contrib_to_force(
                 onika::IntConst<twojmax> _twojmax_, int idxu_max, int jelem
               , double wj_jj, double rcutij_jj, double sinnerij_jj, double dinnerij_jj
               , double x, double y, double z, double z0, double r, double rsq
               , double const * __restrict__ rootpqarray
               , int const * __restrict__ y_jju_map, int idxu_max_alt
               , double rmin0, double rfac0, bool switch_flag, bool switch_inner_flag, bool chem_flag
               , double * __restrict__ fij
               , SnapXSForceExtStorageT& ext )
  {
//    printf("specific add_nbh_contrib_to_force\n");

    const double r0inv = 1.0 / sqrt(r * r + z0 * z0);
    const double a_r = r0inv * z0;
    const double a_i = -r0inv * z;
    const double b_r = r0inv * y;
    const double b_i = -r0inv * x;

    const double rscale0 = rfac0 * M_PI / (rcutij_jj - rmin0);
    const double dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

    double rinv = 1.0 / r;
    double ux = x * rinv;
    double uy = y * rinv;
    double uz = z * rinv;

    const double dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

    double dr0inv[3];
    dr0inv[0] = dr0invdr * ux;
    dr0inv[1] = dr0invdr * uy;
    dr0inv[2] = dr0invdr * uz;

    double dz0[3];
    dz0[0] = dz0dr * ux;
    dz0[1] = dz0dr * uy;
    dz0[2] = dz0dr * uz;

    double da_r[3] , da_i[3] , db_r[3] , db_i[3];
    for (int k = 0; k < 3; k++) {
      da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
      da_i[k] = -z * dr0inv[k];
    }
    da_i[2] += -r0inv;
    for (int k = 0; k < 3; k++) {
      db_r[k] = y * dr0inv[k];
      db_i[k] = -x * dr0inv[k];
    }
    db_i[0] += -r0inv;
    db_r[1] += r0inv;

    const double sfac  = wj_jj * snap_compute_sfac ( rmin0, switch_flag, switch_inner_flag, r, rcutij_jj, sinnerij_jj, dinnerij_jj );
    const double dsfac = wj_jj * snap_compute_dsfac( rmin0, switch_flag, switch_inner_flag, r, rcutij_jj, sinnerij_jj, dinnerij_jj );


    const double da_x_r = da_r[0];
    const double da_x_i = da_i[0];
    
    const double da_y_r = da_r[1];
    const double da_y_i = da_i[1];

    const double da_z_r = da_r[2];
    const double da_z_i = da_i[2];


    const double db_x_r = db_r[0];
    const double db_x_i = db_i[0];
    
    const double db_y_r = db_r[1];
    const double db_y_i = db_i[1];

    const double db_z_r = db_r[2];
    const double db_z_i = db_i[2];

    const auto * __restrict__ Y_r = ext.m_Y_array.r() + jelem * idxu_max_alt;
    const auto * __restrict__ Y_i = ext.m_Y_array.i() + jelem * idxu_max_alt;

#   define CONST_DECLARE(var,value) static constexpr double var = value

#   ifdef SNAP_AUTOGEN_COMPLEX_MATH

    using namespace SnapMath;    
    const SnapMath::Complexd a = { a_r , a_i };
    const SnapMath::Complexd b = { b_r , b_i };
    const SnapMath::Complexd da_x = {da_x_r,da_x_i};
    const SnapMath::Complexd da_y = {da_y_r,da_y_i};
    const SnapMath::Complexd da_z = {da_z_r,da_z_i};
    const SnapMath::Complexd db_x = {db_x_r,db_x_i};
    const SnapMath::Complexd db_y = {db_y_r,db_y_i};
    const SnapMath::Complexd db_z = {db_z_r,db_z_i};    
    const SnapMath::Complex3d da = { da_x , da_y , da_z };
    const SnapMath::Complex3d db = { db_x , db_y , db_z };
    const SnapMath::Double3d u = { ux , uy , uz };
    static constexpr SnapMath::Complexd U_ZERO = {0.,0.};
    static constexpr SnapMath::Complexd U_UNIT = {1.,0.};
    static constexpr SnapMath::Complex3d dU_ZERO = { U_ZERO , U_ZERO , U_ZERO };

#   define U_DECLARE(var)          SnapMath::Complexd var
#   define BAKE_U_BLEND(c,var)     const auto c##_##var = ( conj(c) * var ); \
                                   const auto d##c##_##var = ( conj(d##c) * var + conj(c) * d##var )
#   define U_BLEND(c,var)          c##_##var
#   define U_ASSIGN(var,expr)      var = expr

#   define dU_BLEND(c,var) d##c##_##var
#   define dU_POSTPROCESS(var,jju) d##var = dsfac * var * u + sfac * d##var

#   define CONJC3D(var) conj(var)
#   define dU_DECLARE(var) SnapMath::Complex3d d##var
#   define dU_ASSIGN(var,expr) d##var = expr
#   define dU_STORE_FSKIP(var,xxx,jju_map) /**/
#   define dU_STORE(var,xxx,jju_map)       do{ /*int jju_map=jju; if(y_jju_map!=nullptr) jju_map=y_jju_map[jju]; assert(jju_map==cjjumap);*/ \
      fij[0] +=      d##var.x.r * Y_r[jju_map] + d##var.x.i * Y_i[jju_map];  \
      fij[1] +=      d##var.y.r * Y_r[jju_map] + d##var.y.i * Y_i[jju_map];  \
      fij[2] +=      d##var.z.r * Y_r[jju_map] + d##var.z.i * Y_i[jju_map]; }while(0)
#   define dU_STORE_FHALF(var,xxx,jju_map) do{ /*int jju_map=jju; if(y_jju_map!=nullptr) jju_map=y_jju_map[jju]; assert(jju_map==cjjumap);*/ \
      fij[0] += 0.5*(d##var.x.r * Y_r[jju_map] + d##var.x.i * Y_i[jju_map]); \
      fij[1] += 0.5*(d##var.y.r * Y_r[jju_map] + d##var.y.i * Y_i[jju_map]); \
      fij[2] += 0.5*(d##var.z.r * Y_r[jju_map] + d##var.z.i * Y_i[jju_map]); }while(0)

#   else

    static constexpr double U_UNIT_R = 1.;
    static constexpr double U_UNIT_I = 0.;
    static constexpr double dU_ZERO_X_R = 0. ;
    static constexpr double dU_ZERO_X_I = 0. ;
    static constexpr double dU_ZERO_Y_R = 0. ;
    static constexpr double dU_ZERO_Y_I = 0. ;
    static constexpr double dU_ZERO_Z_R = 0. ;
    static constexpr double dU_ZERO_Z_I = 0. ;

#   define U_DECLARE(var)          double var##_r , var##_i
#   define dU_DECLARE(var) double d##var##_x##_r , d##var##_x##_i , d##var##_y##_r , d##var##_y##_i , d##var##_z##_r , d##var##_z##_i

#   define BAKE_U_BLEND(c,var)     const double c##_##var##_r = c##_r * var##_r + c##_i * var##_i , \
                                                c##_##var##_i = c##_r * var##_i - c##_i * var##_r ; \
                                   const auto d##c##_d##var##_x##_r = d##c##_x##_r * var##_r + d##c##_x##_i * var##_i \
                                                                    + c##_r * d##var##_x##_r + c##_i * d##var##_x##_i; \
                                   const auto d##c##_d##var##_x##_i = d##c##_x##_r * var##_i - d##c##_x##_i * var##_r \
                                                                    + c##_r * d##var##_x##_i - c##_i * d##var##_x##_r; \
                                   const auto d##c##_d##var##_y##_r = d##c##_y##_r * var##_r + d##c##_y##_i * var##_i \
                                                                    + c##_r * d##var##_y##_r + c##_i * d##var##_y##_i; \
                                   const auto d##c##_d##var##_y##_i = d##c##_y##_r * var##_i - d##c##_y##_i * var##_r \
                                                                    + c##_r * d##var##_y##_i - c##_i * d##var##_y##_r; \
                                   const auto d##c##_d##var##_z##_r = d##c##_z##_r * var##_r + d##c##_z##_i * var##_i \
                                                                    + c##_r * d##var##_z##_r + c##_i * d##var##_z##_i; \
                                   const auto d##c##_d##var##_z##_i = d##c##_z##_r * var##_i - d##c##_z##_i * var##_r \
                                                                    + c##_r * d##var##_z##_i - c##_i * d##var##_z##_r
#   define U_BLEND_R(c,var)        c##_##var##_r
#   define U_BLEND_I(c,var)        c##_##var##_i
#   define CONJ_R(r_expr) (r_expr)
#   define CONJ_U_R(var) var##_r
#   define CONJ_I(i_expr) -(i_expr)
#   define CONJ_U_I(var) -var##_i

#   define U_ASSIGN_R(var,expr) var##_r = expr
#   define U_ASSIGN_I(var,expr) var##_i = expr

#   define CONJC3D_X_R(r_expr) (r_expr)
#   define CONJC3D_dU_X_R(var) var##_x##_r
#   define CONJC3D_X_I(i_expr) -(i_expr)
#   define CONJC3D_dU_X_I(var) -var##_x##_i
#   define CONJC3D_Y_R(r_expr) (r_expr)
#   define CONJC3D_dU_Y_R(var) var##_y##_r
#   define CONJC3D_Y_I(i_expr) -(i_expr)
#   define CONJC3D_dU_Y_I(var) -var##_y##_i
#   define CONJC3D_Z_R(r_expr) (r_expr)
#   define CONJC3D_dU_Z_R(var) var##_z##_r
#   define CONJC3D_Z_I(i_expr) -(i_expr)
#   define CONJC3D_dU_Z_I(var) -var##_z##_i
 
#   define dU_BLEND_X_R(c,var) d##c##_d##var##_x##_r
#   define dU_BLEND_X_I(c,var) d##c##_d##var##_x##_i
#   define dU_BLEND_Y_R(c,var) d##c##_d##var##_y##_r
#   define dU_BLEND_Y_I(c,var) d##c##_d##var##_y##_i
#   define dU_BLEND_Z_R(c,var) d##c##_d##var##_z##_r
#   define dU_BLEND_Z_I(c,var) d##c##_d##var##_z##_i

#   define dU_POSTPROCESS(var,jju) \
        d##var##_x##_r = dsfac * var##_r * ux + sfac * d##var##_x##_r; \
        d##var##_x##_i = dsfac * var##_i * ux + sfac * d##var##_x##_i; \
        d##var##_y##_r = dsfac * var##_r * uy + sfac * d##var##_y##_r; \
        d##var##_y##_i = dsfac * var##_i * uy + sfac * d##var##_y##_i; \
        d##var##_z##_r = dsfac * var##_r * uz + sfac * d##var##_z##_r; \
        d##var##_z##_i = dsfac * var##_i * uz + sfac * d##var##_z##_i

#   define dU_ASSIGN_X_R(var,expr) d##var##_x##_r = expr
#   define dU_ASSIGN_X_I(var,expr) d##var##_x##_i = expr
#   define dU_ASSIGN_Y_R(var,expr) d##var##_y##_r = expr
#   define dU_ASSIGN_Y_I(var,expr) d##var##_y##_i = expr
#   define dU_ASSIGN_Z_R(var,expr) d##var##_z##_r = expr
#   define dU_ASSIGN_Z_I(var,expr) d##var##_z##_i = expr

#   define dU_STORE_FSKIP(var,xxx,jju_map) /**/
#   define dU_STORE(var,xxx,jju_map) do{ /*int jju_map=jju; if(y_jju_map!=nullptr) jju_map=y_jju_map[jju]; assert(jju_map!=-1);*/ \
      fij[0] +=      d##var##_x##_r * Y_r[jju_map] + d##var##_x##_i * Y_i[jju_map];  \
      fij[1] +=      d##var##_y##_r * Y_r[jju_map] + d##var##_y##_i * Y_i[jju_map];  \
      fij[2] +=      d##var##_z##_r * Y_r[jju_map] + d##var##_z##_i * Y_i[jju_map]; }while(0)
#   define dU_STORE_FHALF(var,xxx,jju_map) do{ /*int jju_map=jju; if(y_jju_map!=nullptr) jju_map=y_jju_map[jju]; assert(jju_map!=-1);*/ \
      fij[0] += 0.5*(d##var##_x##_r * Y_r[jju_map] + d##var##_x##_i * Y_i[jju_map]); \
      fij[1] += 0.5*(d##var##_y##_r * Y_r[jju_map] + d##var##_y##_i * Y_i[jju_map]); \
      fij[2] += 0.5*(d##var##_z##_r * Y_r[jju_map] + d##var##_z##_i * Y_i[jju_map]); }while(0)

#   endif

    static constexpr int jmax = twojmax/2;
    static_assert( jmax==2 || jmax==3 || jmax==4 );

#   define SNAP_AUTOGEN_NO_UNDEF 1
    if constexpr ( jmax == 2 )
    {
#     include "compute_ui_jmax2.hxx"
    }

    if constexpr ( jmax == 3 )
    {
#     include "compute_ui_jmax3.hxx"
    }

#   undef SNAP_AUTOGEN_NO_UNDEF
    if constexpr ( jmax == 4 )
    {
#     include "compute_ui_jmax4.hxx"
    }

    fij[0] *= 2.;
    fij[1] *= 2.;
    fij[2] *= 2.;
  }

  // Force operator
  template<class SnapConfParamT>
  struct SnapXSForceOp 
  {
    const SnapConfParamT snaconf;
    
    const size_t * const __restrict__ cell_particle_offset = nullptr;
    const double * const __restrict__ beta = nullptr;
    const double * const __restrict__ bispectrum = nullptr;
    
    const double * const __restrict__ coeffelem = nullptr;
    const unsigned int coeffelem_size = 0;
    const unsigned int ncoeff = 0;
   
    const double * const __restrict__ wjelem = nullptr; // data of m_factor in snap_ctx
    const double * const __restrict__ radelem = nullptr;
    const double * const __restrict__ sinnerelem = nullptr;
    const double * const __restrict__ dinnerelem = nullptr;
    const double rcutfac = 0.0;
    const bool eflag = false;
    const bool quadraticflag = false;
    
    // exaStamp conversion specific falgs
    const bool conv_energy_units = true;
    const double conv_energy_factor = UnityConverterHelper::convert(1., "eV");

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      CellParticlesT cells
      ) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      Mat3d& virial,
      CellParticlesT cells
      ) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      CellParticlesT cells,
      GridCellLocksT locks,
      ParticleLockT& lock_a
      ) const
    {
      FakeMat3d virial;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class Mat3dT,class GridCellLocksT, class ParticleLockT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator ()
      (
      int jnum ,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      Mat3dT& virial ,
      CellParticlesT cells,
      GridCellLocksT locks,
      ParticleLockT& lock_a
      ) const
    {
      static constexpr bool compute_virial = std::is_same_v< Mat3dT , Mat3d >;
      static constexpr bool CPAA = gpu_device_execution();
      static constexpr bool LOCK = ! gpu_device_execution();

      assert( ncoeff == static_cast<unsigned int>(snaconf.ncoeff) );

      buf.ext.init( snaconf );

      // energy and force contributions to the particle
      Mat3dT _vir; // default constructor defines all elements to 0
      double _fx = 0.;
      double _fy = 0.;
      double _fz = 0.;

      // start of SNA

      const int itype = type;
      const double radi = radelem[itype];

      int ninside = 0;
      for (int jj = 0; jj < jnum; jj++)
      {
//        const double delx = buf.drx[jj];
//        const double dely = buf.dry[jj];
//        const double delz = buf.drz[jj];
        const double rsq = buf.d2[jj];
	      const int jtype = buf.nbh_pt[jj][field::type];
	      const double cut_ij = ( radi + radelem[jtype] ) * rcutfac;
	      const double cutsq_ij = cut_ij * cut_ij;
        if( rsq < cutsq_ij && rsq > 1e-20 )
        {
          if( ninside != jj ) { buf.copy( jj , ninside ); }
          ninside++;
        }
      }

      const double* __restrict__ betaloc = coeffelem + itype * (snaconf.ncoeff + 1 ) + 1;
      //const int idxu_max = snaconf.idxu_max; // used by macro ULIST_J_A
      
      // compute Ui, Yi for atom I
      // reads sinner, dinner, wj
      // writes ulist and ulisttot
      //const unsigned int ulist_size = jnum * snaconf.idxu_max;
      /*
      snap_compute_ui( snaconf.nelements, snaconf.twojmax, snaconf.idxu_max
                     , snabuf.element, buf.drx,buf.dry,buf.drz, snabuf.rcutij, snaconf.rootpqarray, snabuf.sinnerij, snabuf.dinnerij, snabuf.wj
                     , snaconf.wselfall_flag, snaconf.switch_flag, snaconf.switch_inner_flag, snaconf.chem_flag
                     , snaconf.wself, snaconf.rmin0, snaconf.rfac0
                     , snabuf.ulist_r_ij, snabuf.ulist_i_ij, snabuf.ulisttot_r, snabuf.ulisttot_i // OUTPUTS
                     , ninside, snaconf.chem_flag ? itype : 0);
      
      Content of snap_compute_ui unrolled here => */

//      double UiTot_array_r[ snaconf.idxu_max * snaconf.nelements ];
//      double UiTot_array_i[ snaconf.idxu_max * snaconf.nelements ];
      snap_uarraytot_zero( snaconf.nelements, snaconf.idxu_max, buf.ext.m_UTot_array.r(), buf.ext.m_UTot_array.i() );
      snap_uarraytot_init_wself( snaconf.nelements, snaconf.twojmax, snaconf.idxu_max, snaconf.wself, snaconf.wselfall_flag, buf.ext.m_UTot_array.r(), buf.ext.m_UTot_array.i(), snaconf.chem_flag ? itype : 0 );

      for (int jj = 0; jj < ninside; jj++)
      {
        const int jtype = buf.nbh_pt[jj][field::type];
        const int jelem = snaconf.chem_flag ? jtype : 0 ;

        const double x = buf.drx[jj];
        const double y = buf.dry[jj];
        const double z = buf.drz[jj];
        const double rsq = buf.d2[jj];
        const double r = sqrt(rsq);
        const double rcutij_jj = ( radi + radelem[jtype] ) * rcutfac;
        const double theta0 = (r - snaconf.rmin0) * snaconf.rfac0 * M_PI / (rcutij_jj - snaconf.rmin0);
        const double z0 = r / tan(theta0);

        //double Ui_array_r[snaconf.idxu_max];
        //double Ui_array_i[snaconf.idxu_max];
        //snap_compute_uarray( snaconf.twojmax, snaconf.rootpqarray, Ui_array_r, Ui_array_i, x, y, z, z0, r );

        double sinnerij_jj = 0.0;
        double dinnerij_jj = 0.0;
        if (snaconf.switch_inner_flag) {
          sinnerij_jj = 0.5*(sinnerelem[itype]+sinnerelem[jtype]);
          dinnerij_jj = 0.5*(dinnerelem[itype]+dinnerelem[jtype]);
        }
        const double wj_jj = wjelem[jtype];
        const double sfac_jj = snap_compute_sfac( snaconf.rmin0, snaconf.switch_flag, snaconf.switch_inner_flag, r, rcutij_jj, sinnerij_jj, dinnerij_jj );

        //snap_add_uarraytot( snaconf.twojmax, jelem, snaconf.idxu_max, sfac_jj*wj_jj, Ui_array_r, Ui_array_i, UiTot_array_r /*snabuf.ulisttot_r*/, UiTot_array_i /*snabuf.ulisttot_i*/ );
        snap_add_nbh_contrib_to_uarraytot( snaconf.twojmax, sfac_jj*wj_jj, x,y,z,z0,r, snaconf.rootpqarray, buf.ext.m_UTot_array.r() + snaconf.idxu_max * jelem, buf.ext.m_UTot_array.i() + snaconf.idxu_max * jelem, buf.ext );
      }
      /****************** end of Ui sum computation **********************/



      // for neighbors of I within cutoff:
      // compute Fij = dEi/dRj = -dEi/dRi
      // add to Fi, subtract from Fj
      // scaling is that for type I

      // reads ulisttot
      // writes ylist
      /*
      snap_compute_yi( snaconf.nelements, snaconf.twojmax, snaconf.idxu_max, snaconf.idxz_max
                     , snaconf.idxz, snaconf.idxcg_block, snaconf.cglist
                     , snabuf.ulisttot_r, snabuf.ulisttot_i
                     , snaconf.idxb_max, snaconf.idxb_block, snaconf.bnorm_flag
                     , snabuf.ylist_r, snabuf.ylist_i // OUTPUTS
                     , betaloc );
      Content of snap_compute_yi inlined here =>
      */
//      double Yi_array_r[ snaconf.idxu_max * snaconf.nelements ];
//      double Yi_array_i[ snaconf.idxu_max * snaconf.nelements ];
      snap_zero_yi_array( snaconf.nelements, snaconf.idxu_max_alt, buf.ext.m_Y_array.r(), buf.ext.m_Y_array.i() );
/*
      snap_add_yi_contribution( snaconf.nelements, snaconf.twojmax, snaconf.idxu_max, snaconf.idxz_max
                              , snaconf.idxz, snaconf.idxcg_block, snaconf.cglist
                              , buf.ext.m_UTot_array.r(), buf.ext.m_UTot_array.i()
                              , snaconf.idxb_max, snaconf.idxb_block, snaconf.bnorm_flag
                              , buf.ext.m_Y_array.r(), buf.ext.m_Y_array.i()
                              , betaloc );
*/
      snap_add_yi_contribution_alt( snaconf.nelements, snaconf.twojmax, snaconf.idxu_max, snaconf.idxz_max_alt
                              , snaconf.idxz_alt, snaconf.idxcg_block, snaconf.cglist
                              , snaconf.y_jju_map, snaconf.idxu_max_alt
                              , buf.ext.m_UTot_array.r(), buf.ext.m_UTot_array.i()
                              , snaconf.idxb_max, snaconf.idxb_block, snaconf.bnorm_flag
                              , buf.ext.m_Y_array.r(), buf.ext.m_Y_array.i()
                              , betaloc );

      /******************* end of Yi computation ********************/


      for (int jj = 0; jj < ninside; jj++)
      {
        const int jtype = buf.nbh_pt[jj][field::type];
        const int jelem = snaconf.chem_flag ? jtype : 0 ;

        const double x = buf.drx[jj];
        const double y = buf.dry[jj];
        const double z = buf.drz[jj];
        const double rsq = buf.d2[jj];
        const double r = sqrt(rsq);
        const double rcutij_jj = ( radi + radelem[jtype] ) * rcutfac;
        const double theta0 = (r - snaconf.rmin0) * snaconf.rfac0 * M_PI / (rcutij_jj - snaconf.rmin0);
        const double z0 = r / tan(theta0);

        double sinnerij_jj = 0.0;
        double dinnerij_jj = 0.0;
        if (snaconf.switch_inner_flag) {
          sinnerij_jj = 0.5*(sinnerelem[itype]+sinnerelem[jtype]);
          dinnerij_jj = 0.5*(dinnerelem[itype]+dinnerelem[jtype]);
        }
        const double wj_jj = wjelem[jtype];

        double fij[3];
	      fij[0]=0.;
	      fij[1]=0.;
	      fij[2]=0.;

        add_nbh_contrib_to_force( snaconf.twojmax, snaconf.idxu_max, jelem , wj_jj, rcutij_jj, sinnerij_jj, dinnerij_jj , x, y, z, z0, r, rsq
                                , snaconf.rootpqarray, snaconf.y_jju_map, snaconf.idxu_max_alt
                                , snaconf.rmin0, snaconf.rfac0, snaconf.switch_flag, snaconf.switch_inner_flag, snaconf.chem_flag
                                , fij, buf.ext );

	      fij[0] *= conv_energy_factor;
	      fij[1] *= conv_energy_factor;
	      fij[2] *= conv_energy_factor;

	      Mat3d v_contrib = tensor( Vec3d{ fij[0] , fij[1] , fij[2] }, Vec3d{ buf.drx[jj],buf.dry[jj],buf.drz[jj] } );
        if constexpr ( compute_virial ) { _vir += v_contrib * -1.0; }

        _fx += fij[0];
        _fy += fij[1];
        _fz += fij[2];
        
        size_t cell_b=0, p_b=0;
        buf.nbh.get(jj, cell_b, p_b);
        concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double> (
            locks[cell_b][p_b]
          , cells[cell_b][field::fx][p_b], cells[cell_b][field::fy][p_b], cells[cell_b][field::fz][p_b]
          , -fij[0]                      , -fij[1]                      , -fij[2] );
      }

      if constexpr ( compute_virial )
        concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double,Mat3d> ( lock_a, fx, fy, fz, virial, _fx, _fy, _fz, _vir );
      else
        concurent_add_contributions<ParticleLockT,CPAA,LOCK,double,double,double> ( lock_a, fx, fy, fz, _fx, _fy, _fz );

      if (eflag)
      {
        double _en = 0.;

        const long bispectrum_ii_offset = snaconf.ncoeff * ( cell_particle_offset[buf.cell] + buf.part );

        // evdwl = energy of atom I, sum over coeffs_k * Bi_k
        const double * const coeffi = coeffelem /*[itype]*/;
	      //	coeffelem[ itype * (snaconf.ncoeff + 1 ) + icoeff + 1 ]
	      double evdwl = coeffi[itype * (snaconf.ncoeff + 1)];
	      //	std::cout << "TYpe = " << itype << "e0 = " << evdwl << std::endl;

        // snabuf.copy_bi2bvec();

        // E = beta.B + 0.5*B^t.alpha.B

        // linear contributions

        for (int icoeff = 0; icoeff < snaconf.ncoeff; icoeff++)
        {
          evdwl += betaloc[icoeff] * bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;
	      // evdwl += coeffi[itype * (snaconf.ncoeff + 1) + icoeff+1] * bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;
        }

        // quadratic contributions

        if (quadraticflag)
        {
          int k = snaconf.ncoeff+1;
          for (int icoeff = 0; icoeff < snaconf.ncoeff; icoeff++) 
          {
            double bveci = bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;
            evdwl += 0.5*coeffi[k++]*bveci*bveci;
            for (int jcoeff = icoeff+1; jcoeff < snaconf.ncoeff; jcoeff++) 
            {
              double bvecj = bispectrum[ bispectrum_ii_offset + jcoeff ] /*bispectrum[ii][jcoeff]*/ ;
              evdwl += coeffi[k++]*bveci*bvecj;
            }
          }
        }
        evdwl *= 1.;//scale[itype][itype];
        _en = evdwl; // ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
        
        if( conv_energy_units )
        {
          _en *= conv_energy_factor;
        }

        en += _en;
      }
      
/*
      if( buf.ext.m_U_array.m_ptr != nullptr )
      {
        printf("U_array has been dynamically allocated\n");
      }
      if( buf.ext.m_DU_array.m_ptr != nullptr )
      {
        printf("DU_array has been dynamically allocated\n");
      }
*/

    }
  };

}

namespace exanb
{
  template<class SnapConfParamT>
  struct ComputePairTraits< exaStamp::SnapXSForceOp<SnapConfParamT> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible      = true;
    static inline constexpr bool BufferLessCompatible         = false;
    static inline constexpr bool CudaCompatible               = true;
  };

}

