#pragma xstamp_cuda_enable

#pragma xstamp_grid_variant

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/log.h>
#include <exanb/core/cpp_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/file_utils.h>

#include "n2p2.h"

#include "InterfaceExastamp.h"
#include <exanb/core/physics_constants.h>

#ifdef XNB_CUDA_VERSION
#include <onika/cuda/cuda_context.h>
#include <exanb/core/cpu_gpu_partition.h>
#endif

#include <vector>
#include <memory>
#include <iostream>

// this allows for parallel compilation of templated operator for each available field set



namespace exaStamp
{
  using onika::memory::DEFAULT_ALIGNMENT;
  
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type ,field::_id >
    >
  class N2P2ComputeForce : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( NNPParams             , parameters        , INPUT        , REQUIRED );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0      );
    ADD_SLOT( ParticleSpecies       , species           , INPUT        , REQUIRED );
    ADD_SLOT( int64_t               , timestep          , INPUT        , REQUIRED );
    ADD_SLOT( GridChunkNeighbors    , chunk_neighbors   , INPUT        , GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT        , false    );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT            );
    ADD_SLOT( Domain                , domain            , INPUT        , REQUIRED );



    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    using ComputeBuffer = ComputePairBuffer2<false,false>;
    using CellParticles = typename GridT::CellParticles;
    //    using ParticleLock  = decltype( ComputePairOptionalLocks<false>{}[0][0] );

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type ,field::_id >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type ,field::_id, field::_virial >;
    using ComputeFields              = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};

    /**
    *   N2P2 Interface class
    *   TODO: For later, include directly the n2p2 root folder to exastamp git repository 
    *
    *   Note: Si jamais, la classe c++ se situe dans mon home: 
    *   /ccc/home/cont001/ocre/lasvenesr/src_md/n2p2/src/libnnpif/ExaStamp/InterfaceExastamp.cpp
    */
    std::vector< std::shared_ptr<nnp::InterfaceExastamp> > m_interfaces;
    
  public:
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      if( m_rcut == 0.0 )
      {
        m_rcut = parameters->cutoff; // ang
      }

      *rcut_max = std::max( *rcut_max , m_rcut );
      
      size_t n_cells = grid->number_of_cells();
      if( n_cells == 0 )
      {
        return ;
      }
       
      uint64_t n_particles = grid->number_of_particles() - grid->number_of_ghost_particles(); 
      size_t n_threads = omp_get_max_threads();

      if( n_threads > m_interfaces.size() )
      {
        size_t old_n_threads = m_interfaces.size();
        m_interfaces.resize( n_threads );
        for(size_t i = old_n_threads; i < n_threads; i++)
        {
          assert( m_interfaces[i] == nullptr );

          m_interfaces[i] = std::make_shared<nnp::InterfaceExastamp> ();

#         pragma omp master
          {
            ldbg << NNP_FLAG_INFO << " : n2p2 input directory is "          << parameters->dir          << std::endl;
            ldbg << NNP_FLAG_INFO << " : number of species : "              << species->size()          << std::endl;
            ldbg << NNP_FLAG_INFO << " : show extrapolation warning : "     << parameters->showew       << std::endl;
            ldbg << NNP_FLAG_INFO << " : reset extrapolation warning : "    << parameters->resetew      << std::endl;
            ldbg << NNP_FLAG_INFO << " : show extrapolation warning sum : " << parameters->showewsum    << std::endl;
            ldbg << NNP_FLAG_INFO << " : max extrapolation warning : "      << parameters->maxew        << std::endl;
            ldbg << NNP_FLAG_INFO << " : conversion factor length : "       << parameters->cflength     << std::endl;
            ldbg << NNP_FLAG_INFO << " : conversion factor energy : "       << parameters->cfenergy     << std::endl;
            ldbg << NNP_FLAG_INFO << " : cutoff (n2p2) : "                  << parameters->cutoff       << std::endl;
#           ifndef N2P2_NO_SF_GROUPS            
            ldbg << NNP_FLAG_INFO << " : Using symmetry functions groups"                               << std::endl;
#           else
            ldbg << NNP_FLAG_INFO << " : Not using symmetry functions groups"                           << std::endl;
#           endif
          }                      

          if (!m_interfaces[i]->isInitialized()) {
                m_interfaces[i]->initialize(parameters->dir,
                                            parameters->showew,
                                            parameters->resetew,
                                            parameters->showewsum,
                                            parameters->maxew,
                                            parameters->cflength,
                                            parameters->cfenergy,
                                            parameters->cutoff,
                                            species->size(),
											parameters->dumpout);
          }
          m_interfaces[i]->resizeStructure(n_particles, *timestep);
        }
      }

      size_t tid = omp_get_thread_num();     
      if (*rcut_max < m_interfaces[tid]->getMaxCutoffRadius()) {
        lerr << "ERROR: stamp max radius is lower than n2p2 max radius" << std::endl;
        std::abort();
      }
		
      ForceOp force_op { *rcut_max, m_interfaces, *timestep };

      ComputePairNullWeightIterator          cp_weight{};
      ComputePairOptionalLocks<false>        cp_locks {};
      GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();
      ComputePairTrivialCellFiltering cpu_cell_filter = {};

      if( domain->xform_is_identity() )
      {
        NullXForm cp_xform;
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks, cpu_cell_filter );
        compute_cell_particle_pairs( *grid, m_rcut, *ghost, optional, force_buf, force_op , compute_force_field_set );
      }
      else
      {
        LinearXForm cp_xform { domain->xform() };
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks, cpu_cell_filter );
        compute_cell_particle_pairs( *grid, m_rcut, *ghost, optional, force_buf, force_op , compute_force_field_set );
      }
    } 
    
    private:
    
    struct alignas(DEFAULT_ALIGNMENT) ForceOp 
    {
      const double m_rcut;      
      // One interface per thread
      std::vector< std::shared_ptr<nnp::InterfaceExastamp> > & m_interfaces;
      int64_t timestep;
      
      inline void operator ()
        (
        size_t n,
        ComputeBuffer& buf,
        double& en,
        double& fx,
        double& fy,
        double& fz,
        unsigned int type, // On a besoin du type de l'atome courant
        unsigned int id, // idem pour l'identifiant de l'atome courant
        CellParticles* unused
	//        ComputePairOptionalLocks<false> unused2,
	//        ParticleLock& unused3
        ) const
      {
        Mat3d virial;
        this->operator () ( n,buf,en,fx,fy,fz,type,id,virial, unused );
      }

      inline void operator ()
        (
        size_t n,
        ComputeBuffer& buf,
        double& en,
        double& fx,
        double& fy,
        double& fz,
        unsigned int type,
        unsigned int id,
        Mat3d& virial ,
        CellParticles*
	//        ComputePairOptionalLocks<false>,
	//        ParticleLock&
        ) const
      {
        size_t tid = omp_get_thread_num();
        assert(tid < m_interfaces.size());
        nnp::InterfaceExastamp& interface = *m_interfaces[tid];

        // energy and force contributions to the particle
        double _fx = 0.;	
        double _fy = 0.;
        double _fz = 0.;
	







#       pragma omp simd
        for(unsigned int i=0;i<n;++i)
        {
          buf.d2[i] = std::sqrt( buf.d2[i] );
        }

	// On définit l'atome courant comme central (cf: articles/fonctions de symétrie)
	// en lui donnant le nombre de voisins, ainsi que leurs positions relatives 
	// et distances.
        interface.setCentralAtom(id, type, n, buf.drx, buf.dry, buf.drz, buf.d2);

	// Une fois l'atome central configuré (voisins, positions, etc...)
	// On va calculer pour cet atome, les fonctions de symétrie (groupes),
	// l'énergie de la structure en propageant le vecteur descripteur
	// dans le réseau, ce qui va nous donner l'énergie.
        interface.process();
	
	std::vector<nnp::Vec3D> atomFbis;
	atomFbis.resize(n+1);	
	interface.getParticleForces(n, atomFbis);
        en = interface.getEnergy();	

        Mat3d _vir; // default constructor defines all elements to 0
        assert( _vir.m11==0 && _vir.m12==0 && _vir.m13==0 && _vir.m21==0 && _vir.m22==0 && _vir.m23==0 && _vir.m31==0 && _vir.m32==0 && _vir.m33==0);

        for(unsigned int i=0;i<n+1;++i)
        {
	  const nnp::Vec3D F = atomFbis[i];
	  _vir += tensor( Vec3d{F[0], F[1], F[2]}, Vec3d{buf.drx[i],buf.dry[i],buf.drz[i]} ) * -1.0;
	}
	const nnp::Vec3D FcentralAtom = atomFbis[n];
	fx = -1. * FcentralAtom[0];
	fy = -1. * FcentralAtom[1];
	fz = -1. * FcentralAtom[2];
        virial += _vir;
	  
	// Enfin, on doit réinitialiser certains champs, comme les champs
	// structure.hasDerivatives{Functions, FunctionGroups} pour indiquer qu'il faudra
	// recalculer les fonctions de symétrie au prochain pas de temps.        
        interface.resetCentralAtom();



	
      }    
    };
 
    double m_rcut = 0.0;









  };

  template<class GridT> using N2P2ComputeForceTmpl = N2P2ComputeForce<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "n2p2_force" ,make_grid_variant_operator< N2P2ComputeForceTmpl > );
  }

}


