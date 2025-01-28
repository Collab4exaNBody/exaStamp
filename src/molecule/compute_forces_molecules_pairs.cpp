// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>

#include <exaStamp/molecule/pairs_force_functor.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/math/basic_types_yaml.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/log.h>

#include <string>
#include <iostream>
#include <sstream>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class IntramolecularPairForce : public OperatorNode
  {

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( Domain                  , domain                , INPUT );
    ADD_SLOT( ChemicalPairs           , chemical_pairs        , INPUT_OUTPUT, OPTIONAL );
    ADD_SLOT( ParticleSpecies         , species               , INPUT, REQUIRED );
    ADD_SLOT( BondsPotentialParameters, potentials_for_bonds  , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( GridT                   , grid                  , INPUT_OUTPUT );
    ADD_SLOT( GridParticleLocks       , particle_locks        , INPUT_OUTPUT);

    ADD_SLOT( MoleculeComputeParameterSet  , molecule_compute_parameters , INPUT_OUTPUT, DocString{"Intramolecular functionals' parameters"} );
    ADD_SLOT( IntramolecularParameterIndexLists, intramolecular_parameters , INPUT_OUTPUT, DocString{"Intramolecular functional parmater index lists"} );

    ADD_SLOT( bool                      , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( bool                      , compute_virial      , INPUT , false );
    ADD_SLOT( bool                      , per_atom_charge     , INPUT , true );

  public:
    inline void execute() override final
    {
      using onika::parallel::parallel_for;

      if( ! grid.has_value() || grid->number_of_cells()==0 )
      {
        return;
      }

      if( ! chemical_pairs.has_value() )
      {
        lerr << "chemical_pairs input missing" << std::endl;
        std::abort();
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
      using ChargeFieldT = decltype( grid->field_accessor( field::charge ) );
      using VirialFieldT = decltype( grid->field_accessor( field::virial ) );
      ChargeFieldT charge_field = {};
      VirialFieldT virial_field = {};
      if( *per_atom_charge )
      {
        charge_field = grid->field_accessor( field::charge );
      }
      if( need_virial )
      {
        virial_field = grid->field_accessor( field::virial );
      }
      
      const Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 

#     ifndef NDEBUG
      Vec3d tmp = domain->xform() * size_box;
      if( fabs(tmp.x) <= 1.5 || fabs(tmp.y) <= 1.5 || fabs(tmp.z) <= 1.5 )
      {
        fatal_error()<<"xform="<<domain->xform()<<", size_box="<<size_box<<", tmp="<<tmp<<std::endl;
      }
#     endif

      ldbg<<"n_pairs = "<<chemical_pairs->size()<<" , pair params = "<< molecule_compute_parameters->m_pair_params.size()<<std::endl;

      auto compute_pair_force_opt_virial = [&]( auto cpvir )
      {
        auto cells = grid->cells_accessor();
        IntramolecularPairForceOp<decltype(cells),ChargeFieldT,VirialFieldT,cpvir.value> pair_force_op =
          { cells
          , particle_locks->data()
          , molecule_compute_parameters->m_pair_params.data()
          , intramolecular_parameters->m_pair_param_idx.data()
          , chemical_pairs->data()
          , species->data()
          , size_box , half_min_size_box
          , domain->xform() , domain->xform_is_identity()
          , *per_atom_charge
          , charge_field, virial_field };

        /* auto parallel_op = */ parallel_for( chemical_pairs->size(), pair_force_op, parallel_execution_context("intramol_pair_force") );
        // auto exec_ctrl_obj = parallel_execution_stream() << std::move(parallel_op) ;
      };

      if( need_virial ) compute_pair_force_opt_virial( onika::TrueType{} );
      else              compute_pair_force_opt_virial( onika::FalseType{} );

      ldbg<<"intramolecular_pair_force done"<<std::endl;
    }

  };

  template<class GridT> using IntramolecularPairForceTmpl = IntramolecularPairForce<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_forces_molecules_pairs)
  {
    OperatorNodeFactory::instance()->register_factory( "intramolecular_pair_force", make_grid_variant_operator< IntramolecularPairForceTmpl > );
  }

}

