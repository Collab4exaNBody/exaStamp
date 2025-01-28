// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

#include <yaml-cpp/yaml.h>
#include <memory>
#include <utility>// std::pair

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/math/basic_types_yaml.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/log.h>
#include <exanb/core/thread.h>  // GridParticleLocks

#include <exaStamp/molecule/impropers_potentials_parameters.h>
#include <exaStamp/molecule/periodic_r_delta.h>
#include <exaStamp/molecule/impropers_force_functor.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_ep >
    >
  class ComputeForcesImpropersNode : public OperatorNode
  {    
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( Domain                       , domain                   , INPUT );
    ADD_SLOT( ChemicalImpropers            , chemical_impropers       , INPUT, OPTIONAL );
    ADD_SLOT( ParticleSpecies              , species                  , INPUT, REQUIRED );
    ADD_SLOT( ImpropersPotentialParameters , potentials_for_impropers , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( GridT                        , grid                     , INPUT_OUTPUT );
    ADD_SLOT( GridParticleLocks            , particle_locks           , INPUT_OUTPUT);

    ADD_SLOT( MoleculeComputeParameterSet  , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );
    ADD_SLOT( IntramolecularParameterIndexLists, intramolecular_parameters , INPUT_OUTPUT, DocString{"Intramolecular functional parmater index lists"} );

    ADD_SLOT( bool                      , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( bool                      , compute_virial      , INPUT , false );

  public:
    inline void execute ()  override final
    {
      if( ! grid.has_value() || grid->number_of_cells()==0 )
      {
        return;
      }

      if( ! chemical_impropers.has_value() )
      {
        fatal_error() << "chemical_impropers input missing" << std::endl;
      }

      // do we need to compute energy and virial
      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
      {
        log_energy = *trigger_thermo_state ;
      }
      else
      {
        ldbg << "trigger_thermo_state missing " << std::endl;
      }
      const bool need_virial = log_energy && *compute_virial;
      using VirialFieldT = decltype( grid->field_accessor( field::virial ) );
      VirialFieldT virial_field = {};
      if( need_virial ) virial_field = grid->field_accessor( field::virial );

      ldbg<<"n_impropers = "<<chemical_impropers->size()<<std::endl;

      const Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 

#     ifndef NDEBUG
      if( ! domain->xform_is_identity() )
      {
        Vec3d tmp = domain->xform() * size_box;
        if( fabs(tmp.x) <= 1.5 || fabs(tmp.y) <= 1.5 || fabs(tmp.z) <= 1.5 )
        {
          lerr<<"xform="<<domain->xform()<<", size_box="<<size_box<<", tmp="<<tmp<<std::endl;
          std::abort();
        }
      }
#     endif

      auto compute_impropers_force_opt_virial = [&]( auto cpvir )
      {
        auto cells = grid->cells_accessor();
        ImproperForceOp<decltype(cells),VirialFieldT,cpvir.value> improper_force_op =
          { cells
          , particle_locks->data()
          , molecule_compute_parameters->m_func_params.data()
          , intramolecular_parameters->m_improper_param_idx.data()
          , chemical_impropers->data()
          , size_box , half_min_size_box
          , domain->xform() , domain->xform_is_identity()
          , virial_field };

        /* auto parallel_op = */ parallel_for( chemical_impropers->size(), improper_force_op, parallel_execution_context("improper_force") );
        // auto exec_ctrl_obj = parallel_execution_stream() << std::move(parallel_op) ;
      };

      if( need_virial ) compute_impropers_force_opt_virial( onika::TrueType{} );
      else              compute_impropers_force_opt_virial( onika::FalseType{} );

      ldbg<<"compute_force_improper done"<<std::endl;
    }

  };

  template<class GridT> using ComputeForcesImpropersNodeTmpl = ComputeForcesImpropersNode<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_forces_molecules_impropers)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_force_improper", make_grid_variant_operator< ComputeForcesImpropersNodeTmpl > );
  }

}

