// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

// #pragma xstamp_grid_variant // DO NOT REMOVE THIS LINE

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/domain.h>
#include <exaStamp/molecule/molecule_list.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/molecule_compute_buffer.h>
#include <exaStamp/molecule/periodic_r_delta.h>
#include <exaStamp/molecule/potential_functional.h>
#include <exaStamp/molecule/intramolecular_compute.h>

#include <onika/flat_tuple.h>
#include <onika/parallel/parallel_for.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_ep >
    >
  class IntramolecularForceCompute : public OperatorNode
  {
    using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
    static constexpr bool has_virial_field = has_virial_field_t::value;
    using ComputeBufferT = std::conditional_t< has_virial_field , MoleculeVirialComputeBuffer, MoleculeComputeBuffer >;
  
    ADD_SLOT( GridT                    , grid                        , INPUT_OUTPUT);
    ADD_SLOT( Domain                   , domain                      , INPUT , REQUIRED );
    ADD_SLOT( MoleculeSpeciesVector    , molecules                   , INPUT , MoleculeSpeciesVector{}    , DocString{"Molecule descriptions"} );
    ADD_SLOT( MoleculeLists            , molecule_list               , INPUT , MoleculeLists{}            , DocString{"List of atoms for each molecule"} );
    ADD_SLOT( MoleculeSetComputeParams , molecule_compute_parameters , INPUT , MoleculeSetComputeParams{} , DocString{"Intramolecular functionals' parameters"} );

  public:
    inline void execute ()  override final
    {
      using onika::parallel::parallel_for;
    
      if( grid->number_of_cells() == 0 || molecules->m_molecules.empty() || molecule_list->number_of_molecules() == 0 ) return;
    
      auto cells = grid->cells();
      const size_t nmol = molecule_list->number_of_molecules();
      const Vec3d size_box { std::abs(domain->extent().x - domain->origin().x)
                           , std::abs(domain->extent().y - domain->origin().y)
                           , std::abs(domain->extent().z - domain->origin().z) };
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 
      const auto xform = domain->xform();
      
      size_t n_cells = grid->number_of_cells();
      size_t nghosts = 0;
      for(size_t i=0;i<n_cells;i++) if(grid->is_ghost_cell(i)) nghosts += grid->cell(i).size();

      ldbg << "ghost particles = "<< nghosts <<std::endl;
      ldbg << "number of molecules = " << nmol << std::endl;
      ldbg << "virial = " << std::boolalpha << has_virial_field << std::endl;
      ldbg << "molecule buffer size = " << sizeof(ComputeBufferT) << std::endl;
      
      auto fields = onika::make_flat_tuple(
        grid->field_const_accessor(field::rx) ,
        grid->field_const_accessor(field::ry) , 
        grid->field_const_accessor(field::rz) ,
        grid->field_accessor(field::fx) , 
        grid->field_accessor(field::fy) , 
        grid->field_accessor(field::fz) , 
        grid->field_accessor(field::ep) );

      auto intramol_op = make_intramolecular_functor( *molecule_list, *molecules, *molecule_compute_parameters , xform, size_box, half_min_size_box, grid->cells_accessor() , fields );
      
      parallel_for( nmol , intramol_op , parallel_execution_context() );
/*
#     pragma omp parallel
      {              
#       pragma omp for schedule(dynamic)        
        for(size_t m=0;m<nmol;m++)
        {
          func( m );
        } // end of loop on molecules       
      } // end of parallel section
*/

    }

  };

  template<class GridT> using IntramolecularForceComputeTmpl = IntramolecularForceCompute<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("intramolecular_compute", make_grid_variant_operator< IntramolecularForceComputeTmpl > );
  }

}
