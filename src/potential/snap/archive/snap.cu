



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
#include <onika/cpp_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/file_utils.h>

#include <exaStamp/potential/snap/snap_params.h>
#include <exaStamp/potential/snap/snap_read_lammps.h>
#include <exaStamp/potential/snap/snap_config.h>

#include "snapCg.h"
#include "snapBs.h"

#ifdef XNB_CUDA_VERSION
#include <onika/cuda/cuda_context.h>
#include "snap_gpu.h"
#endif

#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <vector>
#include <memory>
#include <iostream>

// this allows for parallel compilation of templated operator for each available field set



namespace exaStamp
{

  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;
  using namespace SnapExt;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class SnapComputeForce : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( SnapParms             , parameters        , INPUT , REQUIRED );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT , false );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
    ADD_SLOT( GridParticleLocks     , particle_locks    , INPUT        , OPTIONAL , DocString{"particle spin locks"} );

    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    static constexpr bool UseWeights = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors>;
    using CellParticles = typename GridT::CellParticles;
//   using ParticleLock = decltype( ComputePairOptionalLocks<false>{}[0][0] );

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};
        
  public:
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      ldbg << "rcut="<<m_rcut <<std::endl << std::flush;
      if( m_rcut == 0.0 )
      {
        std::string lammps_param = onika::data_file_path( parameters->lammps_param );
        std::string lammps_coef = onika::data_file_path( parameters->lammps_coef ); 
        ldbg << "Snap: read lammps files "<<lammps_param<<" and "<<lammps_coef<<std::endl << std::flush;
        snap_read_lammps(lammps_param, lammps_coef, m_config);
        ldbg <<"rfac0="<<m_config.rfac0() <<", rmin0="<<m_config.rmin0() <<", rcutfac="<<m_config.rcutfac() <<", twojmax="<<m_config.twojmax()<<", bzeroflag="<<m_config.bzeroflag()<<", nmat="<<m_config.materials().size()<<std::endl << std::flush;
        m_rcut = m_config.rcutfac(); // because LAMMPS uses angstrom while exastamp uses nm
	      m_rfac0 = m_config.rfac0();
	      m_rmin0 = m_config.rmin0();	
	      m_bzflag = m_config.bzeroflag();	
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
          //std::cout<<"coef["<<i<<"] = "<<m_coefs[i]<<std::endl;
        }
        
        ldbg<<"compute Cg with jmax="<<jmax<<", NT="<<m_cg_nt<<std::endl;
        m_cg = std::make_shared<snapCg>( m_config.twojmax()*0.5 , m_cg_nt );
        m_cg->compute();        
      }
          
      size_t nt = omp_get_max_threads();
      if( nt > m_thread_ctx.size() )
      {
        size_t old_nt = m_thread_ctx.size();
        m_thread_ctx.resize( nt );
        for(size_t i=old_nt;i<nt;i++)
        {
          assert( m_thread_ctx[i].m_snapbs == nullptr );
          m_thread_ctx[i].m_snapbs = std::make_shared<snapBs>( m_cg->get_jmax(), *m_cg, m_coefs.data(), m_factor[0] );
        }
      }

      ForceOp force_op { m_thread_ctx, m_rcut, m_rfac0, m_rmin0, m_bzflag };

      ComputePairNullWeightIterator cp_weight{};
      
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();

#     ifdef XNB_CUDA_VERSION
      static SnapGPUContext<SnapExt::CUDA_BLOCK_SIZE,3> snap_gpu_jmax3;
      static SnapGPUContext<SnapExt::CUDA_BLOCK_SIZE,4> snap_gpu_jmax4;
      bool go_gpu = false;
      if( global_cuda_ctx() != nullptr ) go_gpu = global_cuda_ctx()->has_devices() && ( m_cg->get_jmax()==3 || m_cg->get_jmax()==4 );
      if( go_gpu )
      {
#       pragma omp critical(cuda_snap_alloc)
        {
          switch( int( m_cg->get_jmax() ) )
          {
            case 3 :
              if( snap_gpu_jmax3.d_bs_fblock == nullptr ) snap_gpu_jmax3.initialize( * (global_cuda_ctx()) , *(m_thread_ctx[0].m_snapbs) );
              break;
            case 4 :
              if( snap_gpu_jmax4.d_bs_fblock == nullptr ) snap_gpu_jmax4.initialize( * (global_cuda_ctx()) , *(m_thread_ctx[0].m_snapbs) );
              break;
            default:
              std::abort();
              break;
          }
        }

        assert( !m_thread_ctx.empty() );
        ldbg<<"going GPU ..."<<std::endl;

        //ProfilingTimer timer;
        //profiling_timer_start(timer);
        
        LinearXForm cp_xform { domain->xform() };
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, ComputePairOptionalLocks<false>{}  );
       // cuda_snap_force( snap_gpu_context, *(m_thread_ctx[0].m_snapbs), *grid, m_rcut, m_factor[0], m_coefs[0], m_rfac0, m_rmin0, m_bzflag, *ghost, optional , compute_force_field_set );
        switch( int( m_cg->get_jmax() ) )
        {
          case 3 :
            cuda_snap_force( snap_gpu_jmax3, *(m_thread_ctx[0].m_snapbs), *grid, m_rcut, m_factor[0], m_coefs[0], m_rfac0, m_rmin0, m_bzflag, *ghost, optional , compute_force_field_set );
            snap_gpu_jmax3.synchronize();
            break;
          case 4 :
            cuda_snap_force( snap_gpu_jmax4, *(m_thread_ctx[0].m_snapbs), *grid, m_rcut, m_factor[0], m_coefs[0], m_rfac0, m_rmin0, m_bzflag, *ghost, optional , compute_force_field_set );
            snap_gpu_jmax4.synchronize();
            break;
          default: std::abort(); break;
        }        

        //parallel_execution_context()->account_gpu_execution_time ( profiling_timer_elapsed_restart(timer) );
      }
      else
#     endif
      {
        LinearXForm cp_xform { domain->xform() };
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, ComputePairOptionalLocks<true>{ particle_locks->data() } );
        compute_cell_particle_pairs( *grid, m_rcut, *ghost, optional, force_buf, force_op, compute_force_field_set, parallel_execution_context()  );
      }

    }

    private:

    struct PerThreadContext
    {
      std::shared_ptr<snapBs> m_snapbs = nullptr;
    };

    // Reaction Field Compute functor
    struct alignas(DEFAULT_ALIGNMENT) ForceOp 
    {
      std::vector<PerThreadContext>& m_thread_ctx;
      const double m_rcut;
      const double m_rfac0;
      const double m_rmin0;
      const size_t m_bzflag;

      template<class CellsAccessorT, class GridCellLocksT, class ParticleLockT>
      inline void operator ()
        (
        size_t n,
        ComputeBuffer& buf,
        double& en,
        double& fx,
        double& fy,
        double& fz,
        CellsAccessorT cells,
        GridCellLocksT locks,
        ParticleLockT& lock_a
        ) const
      {
        FakeMat3d virial;
        this->operator () ( n,buf,en,fx,fy,fz,virial, cells,locks,lock_a );
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
        CellsAccessorT cells
        ) const
      {
        FakeMat3d virial;
        ComputePairOptionalLocks<false> locks = {};
        FakeParticleLock lock_a = {};
        this->operator () ( n,buf,en,fx,fy,fz,virial, cells, locks , lock_a );
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
        Mat3d& virial ,
        CellsAccessorT cells
        ) const
      {
        ComputePairOptionalLocks<false> locks = {};
        FakeParticleLock lock_a = {};
        this->operator () ( n,buf,en,fx,fy,fz,virial, cells, locks , lock_a );
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
        Mat3dT& virial ,
        CellsAccessorT cells,
        GridCellLocksT locks,
        ParticleLockT& lock_a
        ) const
      {
        static constexpr bool compute_virial = std::is_same_v< Mat3dT , Mat3d >;

        size_t tid = omp_get_thread_num();
        assert( tid < m_thread_ctx.size() );
        snapBs& snap_bs = * m_thread_ctx[tid].m_snapbs;

        // energy and force contributions to the particle
        double _en = 0.;
        double _fx = 0.;
        double _fy = 0.;
        double _fz = 0.;

#       pragma omp simd
        for(unsigned int i=0;i<n;++i)
        {
          buf.d2[i] = std::sqrt( buf.d2[i] );
        }

        snap_bs.set_neighbours( buf.drx, buf.dry, buf.drz, buf.d2, m_rcut , n);
        snap_bs.compute_cmm(m_rcut, m_rfac0, m_rmin0);
        snap_bs.compute_bs();
        //snap_bs.compute_bs0();	

        // SNAP energy of the atom
        const double e_tot = snap_bs.en_val();
        double e = 0.;

	      if (m_bzflag == 1)
	      {
	        // SNAP energy zero of the atom
	        const double e_zero = snap_bs.en_zero_val();
	        // SNAP energy of the atom (retrieve e0 = sum b_k B0_k) where B0 is the bispectrum of an isolated atom without neighbors
	        e = (e_tot-e_zero) / n;
	      }
	      else
	      {
	        e = e_tot / n;		
	      }
	
        Mat3dT _vir; // default constructor defines all elements to 0
        //assert( _vir.m11==0 && _vir.m12==0 && _vir.m13==0 && _vir.m21==0 && _vir.m22==0 && _vir.m23==0 && _vir.m31==0 && _vir.m32==0 && _vir.m33==0);

        for(unsigned int i=0;i<n;++i)
        {
          const double3d F = snap_bs.force_val(i);

          auto v_contrib = tensor( Vec3d{F.x,F.y,F.z}, Vec3d{buf.drx[i],buf.dry[i],buf.drz[i]} );

          _fx += F.x;
          _fy += F.y;
          _fz += F.z;
          _en += e;
          _vir += v_contrib * -1.0;

          size_t cell_b=0, p_b=0;
          buf.nbh.get(i, cell_b, p_b);

          auto& lock_b = locks[cell_b][p_b];
          lock_b.lock();
          cells[cell_b][field::fx][p_b] += F.x;
          cells[cell_b][field::fy][p_b] += F.y;
          cells[cell_b][field::fz][p_b] += F.z;
          // cells[cell_b][field::en][p_b] += 0.5 * e // ???
          if constexpr ( compute_virial ) { cells[cell_b][field::virial][p_b] += v_contrib ; }
          lock_b.unlock();
        }

        en += _en; // written only by central atom, no symmetrical contribution
        lock_a.lock();
        fx -= _fx;
        fy -= _fy;
        fz -= _fz;
        virial += _vir;
        lock_a.unlock();
      }
    };

    SnapConfig m_config;
    std::vector<PerThreadContext> m_thread_ctx;
    std::shared_ptr<snapCg> m_cg = nullptr;
    std::vector<double> m_coefs;
    std::vector<double> m_factor;
    double m_rcut = 0.0;
    double m_rfac0 = 0.99363;
    double m_rmin0 = 0.;
    int m_cg_nt = 2;
    bool m_bzflag = true;
  };

  template<class GridT> using SnapComputeForceTmpl = SnapComputeForce<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(snap)
  {
    OperatorNodeFactory::instance()->register_factory( "snap_force" ,make_grid_variant_operator< SnapComputeForceTmpl > );
  }

}


