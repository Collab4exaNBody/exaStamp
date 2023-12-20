#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/compact_grid_pair_weights.h>
#include <exanb/compute/compute_pair_optional_args.h>

#include <memory>

#include <exanb/particle_neighbors/grid_particle_neighbors.h>

namespace exaStamp
{

  // =====================================================================
  // ========================== TestWeight ========================
  // =====================================================================

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_idmol>
    >
  struct TestWeightNode : public OperatorNode
  {
    using WeightType = std::map<std::string, std::vector<double> >;
    using NeighborPairWeight = std::vector< std::vector< double > >;
    using GridParticleNeighbors = exanb::GridParticleNeighbors;

    ADD_SLOT( GridT                , grid          , INPUT  );
    ADD_SLOT( GridParticleNeighbors, primary_neighbors , INPUT  );
    //ADD_SLOT( NeighborPairWeight   , nbh_weight    , INPUT  );
    ADD_SLOT( CompactGridPairWeights , compact_nbh_weight      , INPUT , OPTIONAL );

    void execute() override final
    {
      const GridT& grid = *(this->grid);
      const GridParticleNeighbors& nbh = *primary_neighbors;
      //const NeighborPairWeight& nbh_weight = *(this->nbh_weight);

      auto cells = grid.cells();
      size_t n_cells = nbh.size();

      CompactPairWeightIterator cp_weight { compact_nbh_weight->m_cell_weights.data() };
      auto pw_ctx = cp_weight.make_ctx();
      
      for(size_t cell_a=0;cell_a<n_cells;cell_a++)
      {
        size_t n_particles = nbh[cell_a].nbh_start.size();

        const uint32_t* __restrict__ pair_start   = nbh       [cell_a].nbh_start.data();
        const uint64_t* __restrict__ pairs        = nbh       [cell_a].neighbors.data();
        //const double *  __restrict__ pair_weights = nbh_weight[cell_a]           .data();

        for(size_t p_a=0;p_a<n_particles;p_a++)
        {
          uint64_t idmol_a = cells[cell_a][field::idmol][p_a];
          uint64_t id_a = cells[cell_a][field::id][p_a];
          uint8_t type_a = cells[cell_a][field::type][p_a];

          size_t pair_index = 0;
          if( p_a > 0 ) { pair_index = pair_start[p_a-1]; }
          size_t pair_end = pair_start[p_a];

          for(size_t nbh_index=0;pair_index<pair_end;pair_index++ , nbh_index++)
          {
            size_t cell_b=0, p_b=0;
            exanb::decode_cell_particle( pairs[pair_index] , cell_b , p_b);

            uint64_t idmol_b = cells[cell_b][field::idmol][p_b];
            uint64_t id_b = cells[cell_b][field::id][p_b];
            uint8_t type_b = cells[cell_b][field::type][p_b];

            double weight = cp_weight.get(cell_a,p_a,nbh_index,pw_ctx);
            ldbg << "pair id : " << pairs[pair_index] << " ; id a : " << id_a << " ; id b : " << id_b << " type_a : " << static_cast<int>(type_a) << " type_b : " << static_cast<int>(type_b) << ";  weight : " << weight << std::endl;

            if( weight < 1.0 )
            {
              if(idmol_b!=idmol_a) { std::abort(); }
            }
          }
        }
      }
      lout << "Pair weights are OK !" << std::endl;
    }

  };

  template<class GridT> using TestWeightNodeTmpl = TestWeightNode<GridT>;

   // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "debug_pair_weight", make_grid_variant_operator< TestWeightNodeTmpl > );
  }

}
