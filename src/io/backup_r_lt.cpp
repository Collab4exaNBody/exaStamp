#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/position_long_term_backup.h>

#include <mpi.h>
#include <memory>

namespace exaStamp
{
  using namespace exanb;

  template<class GridT,  class = AssertGridHasFields< GridT, field::_id> >
  class PositionBackupLongTerm : public OperatorNode
  {
    ADD_SLOT( MPI_Comm               , mpi         , INPUT , REQUIRED );
    ADD_SLOT( GridT  , grid     , INPUT );
    ADD_SLOT( Domain , domain            , INPUT , REQUIRED );  // 
    ADD_SLOT( PositionLongTermBackup , backup_r_lt , INPUT_OUTPUT);

  public:
    inline void execute ()  override final
    {    
      GridT& grid = *(this->grid);
      size_t n_cells = grid.number_of_cells();
      IJK dims = grid.dimension();
      auto cells = grid.cells();
      ssize_t gl = grid.ghost_layers();

      Mat3d xform = domain->xform();

      backup_r_lt->m_cell_offset.assign( n_cells+1 , 0 );
      unsigned long long total_particles = 0;
      GRID_FOR_BEGIN(dims,i,loc)
      {
        size_t n_particles = 0;
        if( ! grid.is_ghost_cell(loc) ) { n_particles = cells[i].size(); }
        backup_r_lt->m_cell_offset[i] = total_particles;
        total_particles += n_particles;
      }
      GRID_FOR_END

      const unsigned long long uninitialized_id = std::numeric_limits<unsigned long long>::max() ;

      backup_r_lt->m_cell_offset[n_cells] = total_particles;
      backup_r_lt->m_ids.assign( total_particles , uninitialized_id );
      backup_r_lt->m_positions.assign( total_particles , Vec3d{0.,0.,0.} );
      backup_r_lt->m_xform = xform;

      unsigned long long idMin = uninitialized_id;
      unsigned long long idMax = 0;
      uint64_t inner_particles = 0;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) reduction(min:idMin) reduction(max:idMax) reduction(+:inner_particles) )
        {
          size_t i = grid_ijk_to_index( dims , loc + gl );
          const size_t n_particles = cells[i].size();
          size_t start = backup_r_lt->m_cell_offset[i];
          const auto* __restrict__ id = cells[i][field::id];
          const auto* __restrict__ rx = cells[i][field::rx];
          const auto* __restrict__ ry = cells[i][field::ry];
          const auto* __restrict__ rz = cells[i][field::rz];

          for(size_t j=0;j<n_particles;j++)
          {
            unsigned long long pid = id[j];
            assert( pid != uninitialized_id );
            idMin = std::min( idMin , pid );
            assert( idMin != uninitialized_id );
            idMax = std::max( idMax , pid );
            backup_r_lt->m_ids[start+j] = pid;
            backup_r_lt->m_positions[start+j] = Vec3d { rx[j] , ry[j] , rz[j] };
          }
          inner_particles += n_particles;
        }
        GRID_OMP_FOR_END
      }
      
      assert( inner_particles == total_particles );
      assert( inner_particles == grid.number_of_particles()-grid.number_of_ghost_particles() );

      MPI_Allreduce(MPI_IN_PLACE , &total_particles , 1 , MPI_UNSIGNED_LONG_LONG , MPI_SUM , *mpi);
      MPI_Allreduce(MPI_IN_PLACE , &idMin , 1 , MPI_UNSIGNED_LONG_LONG , MPI_MIN , *mpi);
      MPI_Allreduce(MPI_IN_PLACE , &idMax , 1 , MPI_UNSIGNED_LONG_LONG , MPI_MAX , *mpi);     

      assert( total_particles==0 || idMin<idMax );

      backup_r_lt->m_idMin = idMin;
      backup_r_lt->m_idMax = idMax+1; // exclusive max range is [idMin;idMax+1[

#     ifndef NDEBUG
      for(auto x : backup_r_lt->m_ids)
      {
        assert( x != uninitialized_id );
        assert( x>=idMin && x<=idMax );
      }
#     endif
    }

  };

  template<class GridT> using PositionBackupLongTermTmpl = PositionBackupLongTerm<GridT>;

 // === register factories ===  
  ONIKA_AUTORUN_INIT(backup_r_lt)
  {
   OperatorNodeFactory::instance()->register_factory( "backup_r_lt", make_grid_variant_operator< PositionBackupLongTermTmpl > );
  }

}

