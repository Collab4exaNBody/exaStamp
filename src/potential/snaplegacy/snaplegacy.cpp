#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <exanb/core/cpp_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/file_utils.h>

#include <exaStamp/potential/snap/snap_params.h>

#include <exaStamp/potential/snaplegacy/SnapLegacyCG.h>
#include <exaStamp/potential/snaplegacy/SnapLegacyBS.h>
#include <exaStamp/potential/snaplegacy/snap_legacy_read_lammps.h>
#include <exaStamp/potential/snaplegacy/snap_legacy_config.h>

#include <exanb/core/config.h> // for MAX_PARTICLE_NEIGHBORS constant
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <vector>
#include <memory>
#include <iostream>

//#define XSTAMP_DUMP_BISPECTRUM 1

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type >
    >
  class SnapLegacyForce : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( SnapParms             , parameters        , INPUT , REQUIRED );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT , false );
    ADD_SLOT( bool                  , symmetric_forces  , INPUT , false );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
    ADD_SLOT( GridParticleLocks     , particle_locks    , INPUT        , OPTIONAL , DocString{"particle spin locks"} );

    // additional storage space added to compute buffer created by compute_cell_particle_pairs
    struct alignas(DEFAULT_ALIGNMENT) SnapStorageExt
    {
      alignas(DEFAULT_ALIGNMENT) int species[exanb::MAX_PARTICLE_NEIGHBORS];
    };
    
    // functor that populate compute buffer's extended storage for particle charges
    struct CopyParticleType
    {
      template<typename ComputeBufferT, typename FieldArraysT, class NbhDataT>
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, FieldArraysT cells, size_t cell_b, size_t p_b, const NbhDataT& nbh_data) const noexcept
      {
        assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
        tab.ext.species[tab.count] = cells[cell_b][field::type][p_b];
        DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, nbh_data );
      }
    };

    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    static constexpr bool UseWeights = false;
    static constexpr bool UseNeighbors = true;
    static constexpr bool UseLocks = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors,SnapStorageExt,CopyParticleType>;
    using CellParticles = typename GridT::CellParticles;

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz, field::_type ,field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};

        
    public:

    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      if( m_rcut == 0.0 )
      {
        std::string lammps_param = onika::data_file_path( parameters->lammps_param );
        std::string lammps_coef = onika::data_file_path( parameters->lammps_coef ); 
        ldbg << "Snap: read lammps files "<<lammps_param<<" and "<<lammps_coef<<std::endl;
        snap_legacy_read_lammps(lammps_param, lammps_coef, m_config);
        ldbg <<"rfac0="<<m_config.rfac0() <<", rmin0="<<m_config.rmin0() <<", rcutfac="<<m_config.rcutfac() <<", twojmax="<<m_config.twojmax()<<", nmat="<<m_config.materials().size()<<std::endl;
        m_rcut = m_config.rcutfac(); // because LAMMPS uses angstrom while exastamp uses nm
      }      

      *rcut_max = std::max( *rcut_max , m_rcut );
      
      size_t n_cells = grid->number_of_cells();
      if( n_cells==0 )
      {
        return ;
      }

      if( m_cg == nullptr )
      {
        m_cg_nt = parameters->nt;
        
        for( const auto& mat : m_config.materials() )
        {
          ldbg << '\t' << mat.name() << ": radelem="<<mat.radelem()<<", weight="<<mat.weight()<<", ncoefs="<<mat.number_of_coefficients()<<std::endl;
          for(size_t i=0;i<mat.number_of_coefficients();i++)
          {
            ldbg << "\t\t" << mat.coefficient(i) << std::endl;
          }
        }
        
        double jmax = m_config.twojmax()*0.5;
        int nmat = m_config.materials().size();
        
        if( nmat != 1 )
        {
          lerr << "Snap: ERROR: only 1 material is allowed" << std::endl;
          std::abort();
        }
      
        const SnapMaterial& mat = m_config.materials()[0];

        // temporay, enable mutiple species if they all have weight=1. modifications needed for true multimaterial
        m_factor.assign( MAX_PARTICLE_SPECIES, 1.0 );
        m_factor[0] = mat.weight();

        m_coefs.resize( mat.number_of_coefficients() );
        for(size_t i=0;i<mat.number_of_coefficients();i++)
        {
          m_coefs[i] = mat.coefficient(i);
        }
        
        ldbg<<"compute Cg with jmax="<<jmax<<", NT="<<m_cg_nt<<std::endl;
        m_cg = std::make_shared<SnapLegacyCG>( m_config.twojmax()*0.5 , m_cg_nt );
        m_cg->compute();        
      }

      ldbg << "symmetric_forces = " << *symmetric_forces << std::endl;

      size_t nt = omp_get_max_threads();
      if( nt > m_thread_ctx.size() )
      {
        size_t old_nt = m_thread_ctx.size();
        m_thread_ctx.resize( nt );
        for(size_t i=old_nt;i<nt;i++)
        {
          assert( m_thread_ctx[i].m_snapbs == nullptr );
          m_thread_ctx[i].m_snapbs = std::make_shared<SnapLegacyBS>( m_cg->get_jmax(), m_coefs.data(), m_factor.data() );
        }
      }
      ForceOp force_op { *m_cg , & m_thread_ctx, m_rcut, *symmetric_forces };

      ComputePairOptionalLocks<UseLocks> cp_locks { particle_locks->data() };
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();

      if( domain->xform_is_identity() )
      {
        NullXForm cp_xform;
        auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
        compute_cell_particle_pairs( *grid, m_rcut, *ghost, optional, force_buf, force_op , compute_force_field_set , parallel_execution_context() );
      }
      else
      {
        LinearXForm cp_xform { domain->xform() };
        auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
        compute_cell_particle_pairs( *grid, m_rcut, *ghost, optional, force_buf, force_op , compute_force_field_set , parallel_execution_context() );
      }
    }

    private:

    struct PerThreadContext
    {
      std::shared_ptr<SnapLegacyBS> m_snapbs = nullptr;
    };

    // Reaction Field Compute functor
    struct alignas(DEFAULT_ALIGNMENT) ForceOp 
    {
      const SnapLegacyCG& m_cg;
      std::vector<PerThreadContext> * const m_thread_ctx;
      const double m_rcut;
      const bool m_symetric_forces = false;

      template<class CellsAccessorT, class GridCellLocksT, class ParticleLockT>
      inline void operator ()(
        size_t n, ComputeBuffer& buf, double& en, double& fx, double& fy, double& fz, unsigned int type,
        CellsAccessorT cells, GridCellLocksT locks, ParticleLockT& lock_a
        ) const
      {
        FakeMat3d virial;
        this->operator () ( n,buf,en,fx,fy,fz,type, virial, cells , locks, lock_a);
      }

      template<class CellsAccessorT>
      inline void operator () (
        size_t n, ComputeBuffer& buf, double& en, double& fx, double& fy, double& fz, unsigned int type,
        CellsAccessorT cells
        ) const
      {
        FakeMat3d virial;
        ComputePairOptionalLocks<false> locks = {};
        FakeParticleLock lock_a = {};
        this->operator () ( n,buf,en,fx,fy,fz,type,virial, cells, locks , lock_a );
      }

      template<class CellsAccessorT>
      inline void operator ()
        (
        size_t n,
        ComputeBuffer& buf,
        double& en,
        double& fx,
        double& fy,
        double& fz,
        unsigned int type, // to recover particle type
        Mat3d& virial ,
        CellsAccessorT cells
        ) const
      {
        ComputePairOptionalLocks<false> locks = {};
        FakeParticleLock lock_a = {};
        this->operator () ( n,buf,en,fx,fy,fz,type,virial, cells, locks , lock_a );
      }

      template<class CellsAccessorT, class Mat3dT,class GridCellLocksT, class ParticleLockT>
      inline void operator ()
        (
        size_t n,
        ComputeBuffer& buf,
        double& en,
        double& fx,
        double& fy,
        double& fz,
        unsigned int type, // to recover particle type
        Mat3dT& virial ,
        CellsAccessorT cells,
        GridCellLocksT locks,
        ParticleLockT& lock_a
        ) const
      {
        static constexpr bool compute_virial = true; //std::is_same_v< Mat3dT , Mat3d >;

        // get thread specific compute context
        size_t tid = omp_get_thread_num();
        assert( tid < m_thread_ctx->size() );
        SnapLegacyBS& snap_bs = * (*m_thread_ctx)[tid].m_snapbs;

        // energy and force contributions to the particle
        double _en = 0.;
        double _fx = 0.;
        double _fy = 0.;
        double _fz = 0.;

/*
        std::cout.precision(6);
        std::cout << "part #"<<id<<" : "<<std::endl << std::flush;
*/
        snap_bs.set_neighbours( buf.drx, buf.dry, buf.drz, buf.ext.species, m_rcut , n);
        snap_bs.compute_cmm(m_rcut);
        snap_bs.compute_bs( type, m_rcut, m_cg );

#       ifndef NDEBUG
        {
          int nbs = SnapLegacyBS::n_idx_bs( m_cg.get_jmax()*2 );
          for(int i=0;i<nbs;i++)
          {
            complex<double> bs = snap_bs.bs_val(i);
            complex3d dbs = snap_bs.dbs_val(i);
            exanb::ldbg << "BS("<<i<<") = "<<bs.real()<<"+"<<bs.imag()<<"i , dBS = ( "<< dbs.x.real()<<"+"<<dbs.x.imag()<<"i ,"<< dbs.y.real()<<"+"<<dbs.y.imag()<<"i ,"<< dbs.z.real()<<"+"<<dbs.z.imag()<<"i )"<< std::endl;
          }
        }
#       endif

#       ifdef XSTAMP_DUMP_BISPECTRUM
        {
          int nbs = SnapLegacyBS::n_idx_bs( m_cg.get_jmax()*2 );
          Vec3d pos { cells[buf.cell][field::rx][buf.part], cells[buf.cell][field::ry][buf.part], cells[buf.cell][field::rz][buf.part] };
          std::cout << "BISPECTRUM: "<<pos;
          for(int i=0;i<nbs;i++)
          {
            complex<double> bs = snap_bs.bs_val(i);
            complex3d dbs = snap_bs.dbs_val(i);
#           pragma omp critical(dump_bispectrum)
            std::cout <<" "<<bs.real();
          }
          std::cout<<std::endl;
        }
#       endif

        // SNAP energy of the atom
        const double e_tot = std::real(snap_bs.en_val(type));	
        const double e = e_tot / n;

        Mat3d _vir; // default constructor defines all elements to 0
        // assert( _vir.m11==0 && _vir.m12==0 && _vir.m13==0 && _vir.m21==0 && _vir.m22==0 && _vir.m23==0 && _vir.m31==0 && _vir.m32==0 && _vir.m33==0);

        for(unsigned int i=0;i<n;++i)
        {
          const complex3d fc = snap_bs.force_val(i);
          const double3d F = double3d( std::real(fc.x) , std::real(fc.y) , std::real(fc.z) );

#         ifdef USTAMP_POTENTIAL_WITH_WEIGHTS
          const double weight = tab.nbh_data.get(i); //compute_pair_buffer_get_weight( tab, i );
          e *= weight;
          F.x *= weight;
          F.y *= weight;
          F.z *= weight;
#         endif

          _fx += F.x;
          _fy += F.y;
          _fz += F.z;
	        _en += e;

      	  auto v_contrib = tensor( Vec3d{F.x,F.y,F.z}, Vec3d{buf.drx[i],buf.dry[i],buf.drz[i]} );
          _vir += v_contrib * -0.5;

          size_t cell_b=0, p_b=0;
          buf.nbh.get(i, cell_b, p_b);

          if( m_symetric_forces )
          {
            auto& lock_b = locks[cell_b][p_b];
            lock_b.lock();
            cells[cell_b][field::fx][p_b] += F.x;
            cells[cell_b][field::fy][p_b] += F.y;
            cells[cell_b][field::fz][p_b] += F.z;
            // cells[cell_b][field::en][p_b] += 0.5 * e // ???
            if constexpr ( compute_virial ) { cells[cell_b][field::virial][p_b] += v_contrib * 0.5; }
            lock_b.unlock();
          }

        }

        if( m_symetric_forces ) { lock_a.lock(); }
        en += _en;
        fx -= _fx;
        fy -= _fy;
        fz -= _fz;
        virial += _vir;
        if( m_symetric_forces ) { lock_a.unlock(); }
      }
    };

    SnapConfig m_config;
    std::vector<PerThreadContext> m_thread_ctx;
    std::shared_ptr<SnapLegacyCG> m_cg = nullptr;
    std::vector<double> m_coefs;
    std::vector<double> m_factor;
    double m_rcut = 0.0;
    int m_cg_nt = 2;
  };

  template<class GridT> using SnapLegacyForceTmpl = SnapLegacyForce<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(snaplegacy)
  {
    OperatorNodeFactory::instance()->register_factory( "snaplegacy_force" ,make_grid_variant_operator< SnapLegacyForceTmpl > );
  }

}


