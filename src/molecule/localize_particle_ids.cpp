#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/particle_id_translation.h>

#include <exaStamp/molecule/mol_connectivity.h>

#include <random>
#include <mpi.h>

#include <vector>
#include <utility>

namespace exaStamp
{

  template<
      class GridT
    , class = AssertGridHasFields< GridT, field::_cmol >
    >
  class LocalizeParticleIdsOperator : public OperatorNode
  {
    ADD_SLOT(GridT  , grid          , INPUT_OUTPUT );
    ADD_SLOT(double , bond_max_dist , INPUT, REQUIRED ); // molecule bond max distance
    ADD_SLOT(ParticleIdMap , id_map , INPUT, REQUIRED );

  public:

    // TODO: remove SLs that not connected to at least one particle in the central area (not ghost)

    inline void execute () override final
    {
      // grid cells
      auto cells = grid->cells();
      
      IJK dims = grid->dimension();
      
      // resolve ids to local ids (encoded cell index/particle/index)      
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,i,_)
        {
          MoleculeConnectivity * __restrict__ cmol = cells[i][field::cmol];
          double const * __restrict__ rx = cells[i][field::rx];
          double const * __restrict__ ry = cells[i][field::ry];
          double const * __restrict__ rz = cells[i][field::rz];
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            Vec3d r = { rx[j], ry[j], rz[j] };
            for(size_t k=0;k<cmol[j].size();k++)
            {
              cmol[j][k] = global_to_nearest_local_id( cmol[j][k], *id_map, *grid, r, *bond_max_dist );
            }
          }
        }
        GRID_OMP_FOR_END
      }
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return
R"EOF(
localize particle ids in cmol. it means find the particle instance location given its Id and a maximum distance
)EOF";
    }

  };

  template<class GridT> using LocalizeParticleIds = LocalizeParticleIdsOperator<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "localize_particle_ids", make_grid_variant_operator< LocalizeParticleIds > );
  }

}

