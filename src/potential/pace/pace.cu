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
#include <onika/file_utils.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/core/particle_type_id.h>

#include <vector>
#include <memory>
#include <iostream>
//#include <fmt/core.h>

#include <mpi.h>
// #include "cnpy.h"
// #include "wigner/wigner.hpp"
// #include "wigner/wigner_3nj.hpp"


#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"
#include "ace/ace_b_basis.h"

#include "pace_params.h"
#include "pace_config.h"
#include "pace_force_op.h"

namespace exaStamp
{

  bool hasExtension(const std::string& filename, const std::string& extension) {
    if (filename.length() >= extension.length()) {
      return std::equal(extension.rbegin(), extension.rend(), filename.rbegin());
    }
    return false;
  }

  static char const *const elements_pace[] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};
  static constexpr int elements_num_pace = sizeof(elements_pace) / sizeof(const char *);
  
  static int AtomicNumberByName_pace(const char *elname)
  {
    for (int i = 1; i < elements_num_pace; i++)
      if (strcmp(elname, elements_pace[i]) == 0) return i;
    return -1;
  }
  
  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class PaceForce : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( MPI_Comm              , mpi               , INPUT , REQUIRED);
    ADD_SLOT( PaceParams            , parameters        , INPUT , REQUIRED );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT , false );
    ADD_SLOT( bool                  , conv_coef_units   , INPUT , true );
    ADD_SLOT( bool                  , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
    ADD_SLOT( GridParticleLocks     , particle_locks    , INPUT , OPTIONAL , DocString{"particle spin locks"} );

    ADD_SLOT( long                  , timestep          , INPUT , REQUIRED , DocString{"Iteration number"} );
    //    ADD_SLOT( std::string           , bispectrumchkfile , INPUT , OPTIONAL , DocString{"file with reference values to check bispectrum correctness"} );
    ADD_SLOT( ParticleSpecies       , species           , INPUT        , REQUIRED );
    ADD_SLOT( ParticleTypeMap       , particle_type_map , INPUT        , REQUIRED );    

    ADD_SLOT( PaceContext        , pace_ctx          , PRIVATE );
    ADD_SLOT( bool               , pace_init          , INPUT_OUTPUT, false );

    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    static constexpr bool UseWeights = false;
    static constexpr bool UseNeighbors = true;
    //static constexpr bool UseLocks = true;
    //    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors>;
    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors,PaceComputeBuffer,CopyParticleType>;

    using CellParticles = typename GridT::CellParticles;

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    // using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz >;
    // using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_virial>;
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type ,field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};
    static constexpr FieldSet< field::_type> compute_bispectrum_field_set{};
        
  public:
    
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      
      bool recursive = ( (*parameters).pace_algo == "recursive" );
      bool cTildeBasis = false;
      PaceContext PaceCtx = *pace_ctx;
      PaceCtx.aceimpl = new ACEImpl;
      PaceCtx.aceimpl->basis_set = new ACECTildeBasisSet;
      PaceCtx.aceimpl->ace = new ACERecursiveEvaluator();

      ACEBBasisSet bBasisSet;
      ACECTildeBasisSet cTildeBasisSet;

      auto potential_file_name = (*parameters).pace_coef;
      
      if (hasExtension(potential_file_name, ".yaml")) {
        bBasisSet = ACEBBasisSet(potential_file_name);
        cTildeBasisSet = bBasisSet.to_ACECTildeBasisSet();
        *PaceCtx.aceimpl->basis_set = cTildeBasisSet;
        cTildeBasis = true;
      } else {
        cTildeBasis = false;          
        *PaceCtx.aceimpl->basis_set = ACECTildeBasisSet(potential_file_name);
      }
      
      PaceCtx.aceimpl->ace->set_recursive(recursive);
      PaceCtx.aceimpl->ace->element_type_mapping.init((*parameters).nt + 1);
      
      const auto& sp = *species;
      const int n = (*parameters).nt;
      for (int i = 1; i <= n; i++) {
        const char *elemname = sp[i-1].m_name;
        int atomic_number = AtomicNumberByName_pace(elemname);
        if (atomic_number == -1) std::cout << elemname << "is not a valid element" << std::endl;
        SPECIES_TYPE mu = PaceCtx.aceimpl->basis_set->get_species_index_by_name(elemname);
        if (mu != -1) {
          std::cout << "Mapping LAMMPS atom type #"<< i << "("<<elemname<<") -> ACE species type #"<< mu << std::endl;
          PaceCtx.aceimpl->ace->element_type_mapping(i) = mu;
        } else {
          std::cout << "Element "<< elemname << " is not supported by ACE-potential from file " << potential_file_name << std::endl;
        }
      }
      PaceCtx.aceimpl->ace->set_basis(*PaceCtx.aceimpl->basis_set, 1);
      //               *pace_init = true;
      
      double cutoff=0.;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          *rcut_max = std::max( cutoff, PaceCtx.aceimpl->basis_set->radial_functions->cut(i,j) );
        }
      }

      std::cout << "Total number of basis functions" << std::endl;
      
      for (SPECIES_TYPE mu = 0; mu < PaceCtx.aceimpl->basis_set->nelements; mu++) {
        int n_r1 = PaceCtx.aceimpl->basis_set->total_basis_size_rank1[mu];
        int n = PaceCtx.aceimpl->basis_set->total_basis_size[mu];
        std::cout <<"\t"<< PaceCtx.aceimpl->basis_set->elements_name[mu] << ": "<< n_r1 << " (r=1) " << n << " (r>1)" << std::endl;
      }
      
      size_t n_cells = grid->number_of_cells();
      if( n_cells==0 )
        {
          return ;
        }
      
      if( ! particle_locks.has_value() )
        {
          fatal_error() << "No particle locks" << std::endl;
        }
      
      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
        {
          ldbg << "trigger_thermo_state = " << *trigger_thermo_state << std::endl;
          log_energy = *trigger_thermo_state ;
        }
      else
        {
          ldbg << "trigger_thermo_state missing " << std::endl;
        }
      
      // exanb objects to perform computations with neighbors      
      ComputePairNullWeightIterator cp_weight{};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();      
      LinearXForm cp_xform { domain->xform() };
      
      auto compute_opt_locks = [&](auto cp_locks)
      {
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
        PaceForceOp force_op { PaceCtx,
                               ! (*conv_coef_units) };
        compute_cell_particle_pairs( *grid, *rcut_max, *ghost, optional, force_buf, force_op , compute_force_field_set , parallel_execution_context() );
      };
      if( omp_get_max_threads() > 1 ) compute_opt_locks( ComputePairOptionalLocks<true>{ particle_locks->data() } );
      else                            compute_opt_locks( ComputePairOptionalLocks<false>{} );
      
    }

  };

  template<class GridT> using PaceForceTmpl = PaceForce<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(pace)
  {
    OperatorNodeFactory::instance()->register_factory( "pace_force" ,make_grid_variant_operator< PaceForceTmpl > );
  }

}


