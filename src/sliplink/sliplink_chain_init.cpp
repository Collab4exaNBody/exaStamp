#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_random.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/particle_id_constants.h>
#include <exaStamp/sliplink/sliplink.h>
#include <exanb/core/particle_id_translation.h>
#include <exanb/core/integer_range_algorithm.h>

#include <mpi.h>

#include <vector>
#include <algorithm>
#include <random>


namespace exaStamp
{

  template<
    class GridT
    , class = AssertGridHasFields< GridT, field::_id >
    >
  class SlipLinkChainInitOperator : public OperatorNode
  {

    using ParticleTuple = typename GridT::CellParticles::TupleValueType;
    using ParticleTupleInput = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_id >;

    ADD_SLOT(MPI_Comm           , mpi             , INPUT , MPI_COMM_WORLD );
    ADD_SLOT(SlipLinkParameters , sliplink_config , INPUT, REQUIRED ); 
    ADD_SLOT(Domain             , domain          , INPUT );
    ADD_SLOT(GridT              , grid            , INPUT_OUTPUT );

  public:

    inline void execute () override final
    {
      MPI_Comm comm = *mpi;
      size_t nc = sliplink_config->number_of_chains;
      size_t n_beads = sliplink_config->beads_per_chain;
      
      double sigma2 = sliplink_config->sigma2;
      double bond_max_dist = sliplink_config->bond_max_dist;

      int nprocs = 1;
      int rank = 0;
      MPI_Comm_size(comm,&nprocs);
      MPI_Comm_rank(comm,&rank);

      // at initialization, all MPI processes cover the whole grid, but they have a different subset of chains
      IJK dims = domain->grid_dimension();

      grid->clear_particles();
      grid->set_max_neighbor_distance( 0.0 );
      grid->set_offset( IJK{0,0,0} );
      grid->set_origin( domain->bounds().bmin );
      grid->set_cell_size( domain->cell_size() );
      grid->set_dimension( dims );

      // number of chains for this MPI process
      size_t chain_start = sub_range_begin(rank,nprocs,0,nc);
      size_t chain_end = sub_range_end(rank,nprocs,0,nc);
      assert( chain_start == (rank*nc)/nprocs );
      assert( chain_end == ((rank+1)*nc)/nprocs );
      size_t chain_count = chain_end - chain_start;

      // random distributions for beads
      std::uniform_real_distribution<double> uniform_box_x( domain->bounds().bmin.x , domain->bounds().bmax.x );
      std::uniform_real_distribution<double> uniform_box_y( domain->bounds().bmin.y , domain->bounds().bmax.y );
      std::uniform_real_distribution<double> uniform_box_z( domain->bounds().bmin.z , domain->bounds().bmax.z );
      std::normal_distribution<double> gaussian_bead_displ(0.0, sigma2);

      // grid cells
      auto cells = grid->cells();
#     ifndef NDEBUG
      size_t n_cut_bond_dist = 0;
#     endif

      uint64_t bead_id_start = chain_start * n_beads;
      uint64_t bead_id = bead_id_start;

      auto& re = rand::random_engine();

      for(size_t c=0;c<chain_count;c++)
      {
        Vec3d r { uniform_box_x(re) , uniform_box_y(re) , uniform_box_z(re) };        

        for(size_t b=0;b<n_beads;b++)
        {
          assert( bead_id == ( (chain_start+c)*n_beads + b ) );

          IJK loc = domain_periodic_location( *domain, r );
          ParticleTuple t = ParticleTupleInput(r.x,r.y,r.z, bead_id);

          // add particle to grid
          assert( grid->contains(loc) );
          size_t cell_index = grid_ijk_to_index( dims , loc );
          cells[cell_index].push_back( t );

          Vec3d displ { gaussian_bead_displ(re),
                        gaussian_bead_displ(re),
                        gaussian_bead_displ(re) };

          // FIXME: need pseudo-gaussian distribution with compact support instead
          if( norm(displ) > bond_max_dist )
          {
            displ = displ / norm(displ) * bond_max_dist;
#           ifndef NDEBUG
            ++ n_cut_bond_dist;
#           endif
          }

          r = r + displ;
          ++ bead_id;
        }
      }
      assert( bead_id == ( chain_start + chain_count ) * n_beads );
      grid->rebuild_particle_offsets();

#     ifndef NDEBUG
      if( n_cut_bond_dist )
      {
        lout << n_cut_bond_dist << " bead bonds distances have been truncated to "<<bond_max_dist<<std::endl;
      }
#     endif
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Written by Michelin & CEA/DIF
Credits to Claire Lemarchand, Ioanis Tanis, Thierry Carrard

Initializes the bead chains.
Number of chains is controlled by number_of_chains input.
Every chain has beads_per_chain beads.
)EOF";
    }

  };

  template<class GridT> using SlipLinkChainInit = SlipLinkChainInitOperator<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_chain_init", make_grid_variant_operator< SlipLinkChainInit > );
  }

}

