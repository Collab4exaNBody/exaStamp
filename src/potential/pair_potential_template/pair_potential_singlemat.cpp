#pragma xstamp_cuda_enable

#pragma xstamp_grid_variant

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/log.h>
#include <exanb/core/compact_grid_pair_weights.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <iostream>
#include <type_traits>

// this allows for parallel compilation of templated operator for different variant functor : singe pair pot, multi-param pair pot, and rigid molecule pair pot
#define USTAMP_POTENTIAL_KIND_PAIR 0
#define USTAMP_POTENTIAL_KIND_RIGIDMOL 1
#define USTAMP_POTENTIAL_KIND_MULTI_PAIR 2

${VARIANT:
#include "pair_potential_singlemat_simple.h"
#include "pair_potential_singlemat_multiparam.h"
#include "pair_potential_singlemat_rigidmol.h"
}

// generate object and namespace names
#define OPERATOR_NAME_STR USTAMP_STR(OPERATOR_NAME)
#define PRIV_NAMESPACE_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_details)

//----------------------- debugging configuration --------------------------
/*
// example of debugging particle ids
#define DEBUG_ADDITIONAL_FIELDS ,field::_id
#define DEBUG_ADDITIONAL_PARAMETERS int64_t id,
#define DEBUG_ADDITIONAL_PARAMETER_NAMES id,
#define DEBUG_ADDITIONAL_CODE if(id==0) { std::cout<<"id="<<id<<" : nnbh="<<n<<" ep="<<ep<<" _ep="<<_ep<< std::endl; }
*/

// exemple of debugging particle dr norm with cell/particle indices
/*#define DEBUG_PER_NBH_ADDITIONAL_CODE if(tab.cell==41 && tab.part==733) { \
_Pragma("omp critical(dbg_mesg)") std::cout<<tab.cell<<" "<<tab.part<<" "<<tab.count<<" "<<tab.nbh_count<<" "<<i<<" "<<r<< std::endl; }  
#define DEBUG_ADDITIONAL_CODE if(tab.cell==41 && tab.part==733) { \
_Pragma("omp critical(dbg_mesg)") \
std::cout<<tab.cell*10000+tab.part<<" "<<tab.cell<<" "<<tab.part<<" "<<tab.count<<" "<<tab.nbh_count<<" "<<rmin<<" "<<rmax<< std::endl; }
*/

//--------------------------- end of debug configuration ---------------------------------


// enable potential variant code generation only if it is meaningful for this particular potential
#if defined(USTAMP_POTENTIAL_ENABLE_RIGIDMOL) || !defined(USTAMP_POTENTIAL_RIGIDMOL)

namespace exaStamp
{
  using namespace exanb;

  // ------------------------ Operator implementation --------------------------------------  
  template<
    class GridT,
    class = AssertGridHasFields< GridT ,field::_fx ,field::_fy ,field::_fz POTENTIAL_ADDITIONAL_FIELDS DEBUG_ADDITIONAL_FIELDS >
    >
  class OPERATOR_NAME : public OperatorNode
  {
    
#   if defined(USTAMP_POTENTIAL_RIGIDMOL) // Rigid molecule for single molecule type pair, multiple atom parameters

    using ForceOp = PRIV_NAMESPACE_NAME::RigidMolForceOp;
    using PotParameters = PRIV_NAMESPACE_NAME::RigidMolPotentialParameters;
    static inline constexpr bool ComputeBufferStoreNeighbors = true;
    using ComputeBufferExtendedStorage = RigidMoleculePairContext;
    
#   elif defined(USTAMP_POTENTIAL_MULTI_PARAM) // single potential, multiple pair parameters

    using ForceOp = PRIV_NAMESPACE_NAME::PairMultiForceOp;
    using PotParameters = PRIV_NAMESPACE_NAME::PotentialMultiParameters;
    static inline constexpr bool ComputeBufferStoreNeighbors = false;
    using ComputeBufferExtendedStorage = NoExtraStorage;
    
#   else  // single potential, single parameters

    using ForceOp = PRIV_NAMESPACE_NAME::ForceOp;
    using PotParameters = USTAMP_POTENTIAL_PARAMS;
    static inline constexpr bool ComputeBufferStoreNeighbors = false; // in case of debugging, force this one to true
    using ComputeBufferExtendedStorage = NoExtraStorage;
    
#   endif

    // compile time constant indicating if grid has virial field
#   if defined(USTAMP_POTENTIAL_RIGIDMOL)
    using ComputeFields = FieldSet< /* field::_fx ,field::_fy ,field::_fz, field::_ep, field::_couple */ >;
#   else
    static inline constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;
    static inline constexpr bool has_ep_field = GridHasField<GridT,field::_ep>::value;
    using _ComputeFields = FieldSet< field::_fx ,field::_fy ,field::_fz POTENTIAL_ADDITIONAL_FIELDS DEBUG_ADDITIONAL_FIELDS >;
    using _ComputeFieldsVirial    = FieldSet< field::_fx ,field::_fy ,field::_fz ,field::_virial POTENTIAL_ADDITIONAL_FIELDS DEBUG_ADDITIONAL_FIELDS >;
    using _ComputeFieldsEp = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz POTENTIAL_ADDITIONAL_FIELDS DEBUG_ADDITIONAL_FIELDS >;
    using _ComputeFieldsEpVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_virial POTENTIAL_ADDITIONAL_FIELDS DEBUG_ADDITIONAL_FIELDS >;
    using ComputeFields = std::conditional_t< has_ep_field ,
                            std::conditional_t< has_virial_field , _ComputeFieldsEpVirial , _ComputeFieldsEp > ,
                            std::conditional_t< has_virial_field , _ComputeFieldsVirial , _ComputeFields > >;
#   endif

    using CellParticles = typename GridT::CellParticles;

    // attributes processed during computation

    struct ComputePairScratch
    {
      Mat3d xform;
      PairPotentialMinimalParameters pair_params;
      ForceOp cp_force;
      double rcut_max = 0.0; // max neighbor search distance
      USTAMP_POTENTIAL_PARAMS common_parameters = {};
#     if defined(USTAMP_POTENTIAL_RIGIDMOL) || defined(USTAMP_POTENTIAL_MULTI_PARAM)
      onika::memory::CudaMMVector<PotParameters> parameters_storage;
#     endif
      unsigned int species_index = 0;
      bool has_weight = false;
      bool xform_is_identity = false;
      bool compute_ghosts = false;
    };

    // ============ operator slots ====================    
#   if defined(USTAMP_POTENTIAL_RIGIDMOL) || defined(USTAMP_POTENTIAL_MULTI_PARAM)
    ADD_SLOT( USTAMP_POTENTIAL_PARAMS   , common_parameters    , INPUT , OPTIONAL , DocString{"Potential parameters applied to type pairs without user specified parameters"} );
#   endif
    ADD_SLOT( PotParameters             , parameters           , INPUT , PotParameters{} );
    ADD_SLOT( double                    , rcut                 , INPUT , REQUIRED );
    ADD_SLOT( double                    , rcut_max             , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors      , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( ParticleSpecies           , species              , INPUT , OPTIONAL );

    ADD_SLOT( std::string               , type                 , INPUT , OPTIONAL );
    
    ADD_SLOT( bool                      , ghost                , INPUT , false );
    ADD_SLOT( GridT                     , grid                 , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain               , INPUT , REQUIRED );
    ADD_SLOT( CompactGridPairWeights    , compact_nbh_weight   , INPUT , OPTIONAL );
    ADD_SLOT( bool                      , enable_pair_weights  , INPUT, true );
    ADD_SLOT( ComputePairScratch        , compute_pair_scratch , PRIVATE );
    // ===============================================    

  public:

    inline void task_1( const exanb::GridChunkNeighbors& chunk_neighbors, const GridT& grid )
    {
        assert( chunk_neighbors.number_of_cells() == grid.number_of_cells() );
    }

    inline void task_2( double rcut , const ParticleSpecies& species , const std::string& type , PotParameters& parameters , ComputePairScratch& compute_pair_scratch , double& rcut_max, const USTAMP_POTENTIAL_PARAMS & common_parameters
             /*copy-captured*/, const bool species_has_value, const bool type_has_value )
    {
        // inside a task construct, it is legit to access local variables, they're copy captured (i.e. firstprivate) during lambda construction
        if( ! species_has_value )
        {
          lerr_stream() << "no species input, cannot continue" << std::endl;
          std::abort();
        }
        
        // we need to specify a particular species only for mono-material pair potential
        size_t specy_index = 0;
#       if !defined(USTAMP_POTENTIAL_RIGIDMOL) && !defined(USTAMP_POTENTIAL_MULTI_PARAM)
        if( species.size() > 1 && ! type_has_value )
        {
          fatal_error() <<"Error: exactly 1 atom type expected (single material) but found "<<species.size()<<std::endl;
        }
        if( type_has_value )
        {
          std::string specy_name = type;
          for(size_t s=0;s<species.size();s++)
          {
            if( species.at(s).m_name == specy_name ) { specy_index = s; }
          }
          ldbg << "specy_name="<<specy_name<< ", specy_index = "<<specy_index<<std::endl;
        }
#       endif        
        compute_pair_scratch.species_index = specy_index;

        //compute_pair_scratch.rcut = rcut;        
        double param_max_rcut = rcut;
        double rigid_mol_max_radius = 0.0;

#       if defined(USTAMP_POTENTIAL_RIGIDMOL) || defined(USTAMP_POTENTIAL_MULTI_PARAM)
        // this is the default value for potential parameters that will be used when 
        compute_pair_scratch.common_parameters = common_parameters;

        // compute max rcut accross pair potential parameters
        if( parameters.m_user_pot_parameters != nullptr )
        {
          for( const auto& p : *(parameters.m_user_pot_parameters) )
          {
            param_max_rcut = std::max( param_max_rcut , p.second.rcut );
          }
          ldbg<<"param_max_rcut = "<< param_max_rcut<<std::endl;
        }

        for(size_t i=0;i<species.size();i++)
        {
          const auto & sp = species.at(i);
          if( sp.m_rigid_atom_count > 1 )
          {
            for(size_t j=0;j<sp.m_rigid_atom_count;j++)
            {
              rigid_mol_max_radius = std::max( rigid_mol_max_radius , norm( sp.m_rigid_atoms[j].m_pos ) );
            }
          }
        }
        ldbg << "rigidmol max radius = " << rigid_mol_max_radius << std::endl;

#       endif // defined(USTAMP_POTENTIAL_RIGIDMOL) || defined(USTAMP_POTENTIAL_MULTI_PARAM)

        compute_pair_scratch.rcut_max = param_max_rcut + 2 * rigid_mol_max_radius ;

        rcut_max = std::max( rcut_max , compute_pair_scratch.rcut_max );
        ldbg << "compute_pair default rcut=" << rcut << " / rcut_max=" << compute_pair_scratch.rcut_max << " / global rcut_max="<< rcut_max << std::endl;
    }

    inline void task_3( const ParticleSpecies& species , const std::string& type , GridT& grid , PotParameters& parameters, const double rcut, ComputePairScratch& compute_pair_scratch
            /*copy-captured*/, const bool chunk_neighbors_has_value )
    {
        if( grid.number_of_particles() > 0 )
        {
          if( ! chunk_neighbors_has_value )
          {
            lerr_stream() << "no neighbors, cannot continue" << std::endl;
            std::abort();
          }

#         if defined(USTAMP_POTENTIAL_RIGIDMOL) || defined(USTAMP_POTENTIAL_MULTI_PARAM)

#         ifdef USTAMP_POTENTIAL_RIGIDMOL
          using MultiPairParams = PRIV_NAMESPACE_NAME::RigidMolPotentialParameters;
          using SinglePairParam = PRIV_NAMESPACE_NAME::RigidMolPotentialPairParam;
#         endif

#         ifdef USTAMP_POTENTIAL_MULTI_PARAM
          using MultiPairParams = PRIV_NAMESPACE_NAME::PotentialMultiParameters;
          using SinglePairParam = PRIV_NAMESPACE_NAME::PotentialPairParam;
#         endif

          const bool rebuild_force_op_parameters = ( compute_pair_scratch.cp_force.p==nullptr );

          if( rebuild_force_op_parameters )
          {
            ldbg << "rebuild_force_op_parameters" <<std::endl;

#           if defined(USTAMP_POTENTIAL_RIGIDMOL) || defined(USTAMP_POTENTIAL_MULTI_PARAM)
            compute_pair_scratch.parameters_storage.resize(1);
            compute_pair_scratch.parameters_storage[0].m_user_pot_parameters = nullptr;
            compute_pair_scratch.cp_force.p = compute_pair_scratch.parameters_storage.data();
            auto & pot_params = compute_pair_scratch.parameters_storage[0];
#           endif

            // special check for optimizations :
            // all singla atom species must declared first, then all rigid molecule.
            // this way we can account only for the first Ns single atom species to build set of type pairs
            size_t first_rigid_molecule_index = species.size();
            for(size_t i=0;i<species.size();i++)
            {
              if( species[i].m_rigid_atom_count < 1 )
              {
                fatal_error()<<"rigid_atom_count must be 1 or more"<<std::endl;
              }
              if( species[i].m_rigid_atom_count > 1 )
              {
                first_rigid_molecule_index = std::min( first_rigid_molecule_index , i );
              }
              else if( i > first_rigid_molecule_index )
              {
                fatal_error()<<"All single atom species must be declared first (before rigid molecules)"<<std::endl;
              }
            }
            ldbg << "nb species = "<<species.size()<<" , nb single atoms = "<<first_rigid_molecule_index<<std::endl;

#           ifdef USTAMP_POTENTIAL_RIGIDMOL
            size_t total_rigidmol_atoms = 0;
            for(size_t i=0;i<species.size();i++)
            {
              const size_t rma = species.at(i).m_rigid_atom_count;
              pot_params.m_n_atoms[i] = rma;
              pot_params.m_atoms_start[i] = total_rigidmol_atoms;
              ldbg << "Type "<<i<<" : atoms="<<pot_params.m_n_atoms[i]<<" , start="<<pot_params.m_atoms_start[i]<<std::endl;
              if( rma > 1 )
              {
                if( (total_rigidmol_atoms+rma) > MultiPairParams::MAX_ALL_MOLECULE_ATOMS )
                {
                  fatal_error() << "Too many rigid molecule atoms (max="<<MultiPairParams::MAX_ALL_MOLECULE_ATOMS<<")"<<std::endl;
                }
                for(size_t a=0;a<rma;a++)
                {
                  pot_params.m_atoms[total_rigidmol_atoms+a] = species.at(i).m_rigid_atoms[a];
                }
                total_rigidmol_atoms += rma;
              }
            }
            ldbg << "total_rigidmol_atoms = "<<total_rigidmol_atoms<<std::endl;
#           endif

            const unsigned int n_type_pairs = unique_pair_count( first_rigid_molecule_index );
            if( n_type_pairs > MultiPairParams::MAX_TYPE_PAIR_IDS )
            {
              fatal_error() << "Too many atom types for singlemat rigid molecule implementation. Increase MAX_TYPE_PAIR_IDS"<<std::endl;
            }
 
            ldbg << "n_type_pairs = "<<n_type_pairs<<" , m_user_pot_parameters="<< (void*)parameters.m_user_pot_parameters <<std::endl;
            pot_params.m_nb_pair_params = n_type_pairs;
            for(unsigned int pair_id=0;pair_id<n_type_pairs;pair_id++)
            {
              unsigned int type_a=0, type_b=0;
              pair_id_to_type_pair( pair_id, type_a, type_b );
              if( type_a >= first_rigid_molecule_index || type_b >= first_rigid_molecule_index )
              {
                fatal_error() << "Inconsistent atom type or not a single atom type : type_a="<<type_a<<", type_b="<<type_b<<", Nspecies="<<species.size()<<", Nsingle="<<first_rigid_molecule_index<<std::endl;
              }
              pot_params.m_pair_params[pair_id] = SinglePairParam {};
              pot_params.m_pair_params[pair_id].p = compute_pair_scratch.common_parameters;
              pot_params.m_pair_params[pair_id].rcut = rcut; 
              pot_params.m_pair_params[pair_id].pair_params.m_atom_a = species.at(type_a);
              pot_params.m_pair_params[pair_id].pair_params.m_atom_b = species.at(type_b);
              pot_params.m_pair_params[pair_id].ecut = 0.0;
              
              bool pair_pot_params_found = false;
              if( parameters.m_user_pot_parameters != nullptr )
              {
                auto it = parameters.m_user_pot_parameters->find( std::make_pair( species.at(type_a).name() , species.at(type_b).name() ) );
                if( it == parameters.m_user_pot_parameters->end() )
                {
                  it = parameters.m_user_pot_parameters->find( std::make_pair( species.at(type_b).name() , species.at(type_a).name() ) );
                }
                if( it != parameters.m_user_pot_parameters->end() )
                {
                  pair_pot_params_found = true;
                  pot_params.m_pair_params[pair_id].p = it->second.p;
                  pot_params.m_pair_params[pair_id].rcut = it->second.rcut;
                  ldbg << "found parameter for pair "<<species.at(type_a).name() <<"/"<< species.at(type_b).name()<<" , pair_id="<<pair_id<<" , rcut="<<pot_params.m_pair_params[pair_id].rcut<<std::endl;
                }
              }
              if( ! pair_pot_params_found )
              {
                ldbg << "NO parameter found for pair "<<species.at(type_a).name() <<"/"<< species.at(type_b).name()<<" , pair_id="<<pair_id<<" , rcut="<<pot_params.m_pair_params[pair_id].rcut<<" , parameters.m_user_pot_parameters="<<(const void*)parameters.m_user_pot_parameters <<std::endl;
              }

              pot_params.m_pair_params[pair_id].ecut = energy_cutoff( pot_params.m_pair_params[pair_id].p , pot_params.m_pair_params[pair_id].pair_params , pot_params.m_pair_params[pair_id].rcut ); ;
            }
            
            //if( parameters.m_user_pot_parameters != nullptr ) delete parameters.m_user_pot_parameters;
            parameters.m_user_pot_parameters = nullptr;
            
            
            // reduce number of parameters if possible
            std::vector< std::vector<char> > serialized_param_array( n_type_pairs );
            std::vector< int > param_index( n_type_pairs );
            for(unsigned int pair_id=0;pair_id<n_type_pairs;pair_id++)
            {
              param_index[pair_id] = pair_id;
              serialized_param_array[pair_id] = serialize_pair_potential( pot_params.m_pair_params[pair_id].rcut
                                                                        , pot_params.m_pair_params[pair_id].p
                                                                        , USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(pot_params.m_pair_params[pair_id].pair_params) );
            }
            std::sort(param_index.begin(),param_index.end(),
                      [&serialized_param_array](int a,int b)->bool
                      {
                        return serialized_param_array[a] < serialized_param_array[b];
                      });
            std::vector< int > param_map( n_type_pairs , -1 );
            int n_unique_params = -1;
            int plast = -1;
            if(n_type_pairs>0)
            {
              param_map.push_back(0);
              n_unique_params = 1;
              plast = param_index[0];
            }
            for(unsigned int pair_id=1;pair_id<n_type_pairs;pair_id++)
            {
              int pcur = param_index[pair_id];
              if( serialized_param_array[pcur] == serialized_param_array[ plast ] )
              {
                param_map[pair_id] = n_unique_params-1;
              }
              else
              {
                param_map[pair_id] = n_unique_params++;
                plast = pcur;
              }
            }
            ldbg << "n_type_pairs="<<n_type_pairs<<" , n_unique_params="<<n_unique_params<<std::endl;
          }
          ldbg << "compute_pair_scratch.rcut_max="<< compute_pair_scratch.rcut_max << " , rcut="<<rcut << std::endl;
          
#         else

          const size_t specy_index = compute_pair_scratch.species_index;
          compute_pair_scratch.pair_params.m_atom_a = species.at(specy_index);
          compute_pair_scratch.pair_params.m_atom_b = species.at(specy_index);
          const double ecut = energy_cutoff( parameters, compute_pair_scratch.pair_params , rcut );
          compute_pair_scratch.cp_force = ForceOp { parameters , compute_pair_scratch.pair_params , ecut };
          
#         endif
        }
    }

    inline void task_4( const Domain& domain , const exanb::GridChunkNeighbors& chunk_neighbors , const bool ghost , const bool enable_pair_weights , const CompactGridPairWeights& compact_nbh_weight , GridT& grid , ComputePairScratch& compute_pair_scratch
            /*copy-captured*/, OPERATOR_NAME* self, const bool nbh_weight_has_value )
    {
        if( grid.number_of_particles() > 0 )
        {
          compute_pair_scratch.xform_is_identity = domain.xform_is_identity();
          compute_pair_scratch.xform = domain.xform();
          compute_pair_scratch.has_weight = nbh_weight_has_value && enable_pair_weights;
          if( compute_pair_scratch.has_weight ) { compute_pair_scratch.has_weight = ! compact_nbh_weight.empty(); }
          compute_pair_scratch.compute_ghosts = ghost;

          int compute_case = ( (!compute_pair_scratch.xform_is_identity) ? 2 : 0 ) | ( compute_pair_scratch.has_weight ? 1 : 0 ) ;
          switch( compute_case )
          {
            case 0 : compute_force(grid,chunk_neighbors,&compact_nbh_weight,compute_pair_scratch,self->parallel_execution_context(), std::false_type{} , std::false_type{} ); break;
            case 1 : compute_force(grid,chunk_neighbors,&compact_nbh_weight,compute_pair_scratch,self->parallel_execution_context(), std::false_type{} , std::true_type {} ); break;
            case 2 : compute_force(grid,chunk_neighbors,&compact_nbh_weight,compute_pair_scratch,self->parallel_execution_context(), std::true_type {} , std::false_type{} ); break;
            case 3 : compute_force(grid,chunk_neighbors,&compact_nbh_weight,compute_pair_scratch,self->parallel_execution_context(), std::true_type {} , std::true_type {} ); break;
          }
        }
    }

    inline void execute() override final
    {
      const bool species_has_value = species.has_value();      
      const bool type_has_value = type.has_value();
      const bool chunk_neighbors_has_value = chunk_neighbors.has_value();
      const bool nbh_weight_has_value = compact_nbh_weight.has_value();      

      ldbg << "enable_pair_weights="<< (*enable_pair_weights) <<" , compact_nbh_weight.has_value() = "<<compact_nbh_weight.has_value()<<std::endl;

      USTAMP_POTENTIAL_PARAMS def_pot_params = {};
#     if defined(USTAMP_POTENTIAL_RIGIDMOL) || defined(USTAMP_POTENTIAL_MULTI_PARAM)      
      if( common_parameters.has_value() )
      {
        ldbg << "common_parameters has value" << std::endl;
        def_pot_params = *common_parameters;
      }
      else
      {
        ldbg << "common_parameters has NO value" << std::endl;
      }
#     endif
      
#     ifndef NDEBUG
      task_1( *chunk_neighbors , *grid );
#     endif
      task_2( *rcut , *(species.get_pointer()) , *(type.get_pointer()) , *parameters , *compute_pair_scratch , *rcut_max, def_pot_params /*copy-captured*/ ,species_has_value, type_has_value );
      task_3( *(species.get_pointer()) , *(type.get_pointer()) , *grid , *parameters, *rcut, *compute_pair_scratch /*copy-captured*/ ,chunk_neighbors_has_value );
      auto* self = this;
      task_4( *domain , *chunk_neighbors , *ghost , *enable_pair_weights , * compact_nbh_weight.get_pointer() , *grid , *compute_pair_scratch /*copy-captured*/ , self , nbh_weight_has_value );
    }

    static inline LinearXForm make_cp_xform( const Mat3d& xform , std::true_type ) { return LinearXForm { xform }; }
    static inline NullXForm make_cp_xform( const Mat3d& , std::false_type ) { return NullXForm { }; }

    static inline CompactPairWeightIterator make_cp_weight( const CompactGridPairWeights* cgpw , std::true_type ) { return CompactPairWeightIterator { cgpw->m_cell_weights.data() }; }
    static inline ComputePairNullWeightIterator make_cp_weight( const CompactGridPairWeights* , std::false_type ) { return ComputePairNullWeightIterator { }; }

    template<bool HasWeight>
    static inline auto make_cp_buf_factory( std::integral_constant<bool,HasWeight> )
    {
      return ComputePairBufferFactory< ComputePairBuffer2<HasWeight,ComputeBufferStoreNeighbors,ComputeBufferExtendedStorage> > {};
    }

    template<bool HasXForm, bool HasWeight>
    static inline void compute_force(
                               GridT & grid
                             , const exanb::GridChunkNeighbors& chunk_neighbors
                             , const CompactGridPairWeights* cgpw
                             , const ComputePairScratch& compute_pair_scratch
                             , onika::parallel::ParallelExecutionContext * gpu_exec_ctx
                             , std::integral_constant<bool,HasXForm> has_xform
                             , std::integral_constant<bool,HasWeight> has_weights )
    {

      auto optional = make_compute_pair_optional_args( exanb::GridChunkNeighborsLightWeightIt<false> { chunk_neighbors }
                                           , make_cp_weight(cgpw,has_weights)
                                           , make_cp_xform(compute_pair_scratch.xform,has_xform)
                                           , ComputePairOptionalLocks<false> {} );

      compute_cell_particle_pairs(
          grid
        , compute_pair_scratch.rcut_max
        , compute_pair_scratch.compute_ghosts
        , optional
        , make_cp_buf_factory( has_weights )
        , compute_pair_scratch.cp_force
        , ComputeFields{}
        , gpu_exec_ctx );
    }

    // used only once to compute energy cutoff
    template<class PotParamsT>
    static inline double energy_cutoff(const PotParamsT& p, const PairPotentialMinimalParameters & pair_params, double rcut)
    {
      double tmp_e = 0.0;
      if( rcut > 0.0 )
      {
        double tmp_de = 0.0;
        USTAMP_POTENTIAL_COMPUTE(p,pair_params,rcut,tmp_e,tmp_de);
      }
      return tmp_e;
    }

  };


  namespace TemplateHelper
  {
    template<class GridT> using OPERATOR_NAME = ::exaStamp::OPERATOR_NAME < GridT >;
  }
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( OPERATOR_NAME_STR , make_grid_variant_operator< TemplateHelper::OPERATOR_NAME > );
  }

}

#endif // defined(USTAMP_POTENTIAL_ENABLE_RIGIDMOL) || !defined(USTAMP_POTENTIAL_RIGIDMOL)


#undef OPERATOR_NAME
#undef OPERATOR_NAME_STR

