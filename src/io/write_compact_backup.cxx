#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/fields.h>
#include <exanb/core/backup_r.h>
#include <exanb/core/domain.h>


namespace exaStamp
{
/*
*/

  struct CompactBackupHeader
  {
      AABB bounds;
      double cell_size;
      IJK dom_dims;
      IJK subdom_offset;
      IJK subdom_dims;    
  };

  template<typename GridT>
  struct WriteCompactBackup : public OperatorNode
  {
    ADD_SLOT( MPI_Comm    , mpi      , INPUT , MPI_COMM_WORLD, DocString{"MPI communicator"} );
    ADD_SLOT( GridT       , grid     , INPUT , REQUIRED );
    ADD_SLOT( Domain      , domain   , INPUT , REQUIRED );
    ADD_SLOT( std::string , file     , INPUT , std::string("compact.dump") );

    inline void execute ()  override final
    {
      MPI_Comm comm = *mpi;
      int np=1, rank=0;
      MPI_Comm_size(comm,&np);
      MPI_Comm_rank(comm,&rank);

      const IJK dims = grid->dimension();
      const auto cells = grid->cells();
      const double cell_size = domain->cell_size();
      ssize_t gl = grid.ghost_layers();

      const auto xform = domain->xform();

      std::vector<uint32_t> buffer;
      buffer.reserve(65536);
      
      std::string filename = ( *file ) + std::to_string( rank );
      std::ofstream fout( filename );
      
      const AABB bounds = domain->bounds();
      const double cell_size = grid->cell_size();
      const IJK dom_dims = domain->grid_dimension();
      const IJK subdom_offset = grid->offset();
      const IJK subdom_dims = grid->dimension();
      
      fout.write( (const char*) &bounds , sizeof(bounds) );
      fout.write( (const char*) &cell_size , sizeof(cell_size) );
      fout.write( (const char*) &dom_dims , sizeof(dom_dims) );
      fout.write( (const char*) &subdom_offset , sizeof(subdom_offset) );
      fout.write( (const char*) &subdom_dims , sizeof(subdom_dims) );

      GRID_FOR_BEGIN(dims-2*gl,_,loc /*, schedule(dynamic)*/ )
      {
        const size_t i = grid_ijk_to_index( dims , loc + gl );
        const Vec3d cell_origin = grid->cell_position( loc );

        const auto* __restrict__ rx = cells[i][field::rx];
        const auto* __restrict__ ry = cells[i][field::ry];
        const auto* __restrict__ rz = cells[i][field::rz];

        const auto* __restrict__ vx = cells[i].field_pointer_or_null(field::vx);
        const auto* __restrict__ vy = cells[i].field_pointer_or_null(field::vy);
        const auto* __restrict__ vz = cells[i].field_pointer_or_null(field::vz);
        const bool has_velocity = ( vx!=nullptr && vy!=nullptr && vz!=nullptr );

        const uint32_t n_particles = cells[i].size();
        fout.write( (const char*) &n_particles , sizeof(n_particles) );

        buffer.resize( n_particles * 6 );
        for(uint32_t j=0;j<n_particles;j++)
        {
          buffer[j*6+0] = encode_double_u32( rx[j] , cell_origin.x , cell_size );
          buffer[j*6+1] = encode_double_u32( ry[j] , cell_origin.y , cell_size );
          buffer[j*6+2] = encode_double_u32( rz[j] , cell_origin.z , cell_size );
          if( has_velocity )
          {
             * (float*) ( & buffer[j*6+3] ) = vx[j];
             * (float*) ( & buffer[j*6+4] ) = vy[j];
             * (float*) ( & buffer[j*6+5] ) = vz[j];
          }
        }
        fout.write( (const char*) buffer.data() , buffer.size() * sizeof(uint32_t) );

      }
      GRID_FOR_END

      fout.close();
    }

  };
    
 // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "write_compact_backup", make_grid_variant_operator< WriteCompactBackup > );
  }

}

