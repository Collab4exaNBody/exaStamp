#include <exanb/core/domain.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/memory/allocator.h> // for ONIKA_ASSUME_ALIGNED macro
#include <exanb/core/thread.h> // for ONIKA_ASSUME_ALIGNED macro

#include <exanb/core/basic_types_operators.h>
#include <exanb/core/domain.h>
#include <exanb/core/string_utils.h>

#include <exaStamp/sliplink/sliplink.h>
#include <exanb/core/particle_id_translation.h>
#include <exanb/core/particle_id_codec.h>


#include <random>
#include <mpi.h>
#include <vector>
#include <iomanip>

namespace exaStamp
{

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_id >
    >
  class SlipLinkAnalyticsOperator : public OperatorNode
  {
    ADD_SLOT(MPI_Comm           , mpi             , INPUT , MPI_COMM_WORLD );
    ADD_SLOT(SlipLinkParameters , sliplink_config , INPUT, REQUIRED );
    ADD_SLOT(double             , bond_max_stretch , INPUT , 0.5 ); // fraction of bond_max_dist.
    ADD_SLOT(Domain             , domain              , INPUT , REQUIRED);
    ADD_SLOT(GridT              , grid            , INPUT_OUTPUT );
    ADD_SLOT(ParticleIdMap      , id_map          , INPUT , REQUIRED );
    ADD_SLOT(long               , sl_regen_count  , INPUT , REQUIRED );

    ADD_SLOT( long               , timestep            , INPUT, REQUIRED );
    ADD_SLOT( double             , physical_time       , INPUT, REQUIRED );
    ADD_SLOT( bool               , lb_flag             , INPUT, false );
    ADD_SLOT( bool               , move_flag           , INPUT, false );
    ADD_SLOT( bool               , print_header        , INPUT, false );
    ADD_SLOT( double             , lb_inbalance        , INPUT, 0.0 );

  public:

    inline void execute () override final
    {
      MPI_Comm comm = *mpi;
      int nprocs = 1;
      int rank = 0;
      MPI_Comm_size(comm,&nprocs);
      MPI_Comm_rank(comm,&rank);

      const size_t nc = sliplink_config->number_of_chains;
      const size_t n_beads = sliplink_config->beads_per_chain;
      const double bond_max_search_dist = sliplink_config->bond_max_dist * ( 2. + *bond_max_stretch );

      // ldbg << "sliplinks analytics : bond_max_search_dist = " << bond_max_search_dist << std::endl;

      auto cells = grid->cells();
      size_t gl = grid->ghost_layers();          
      IJK dims = grid->dimension();

      std::vector< Vec3d > ree( nc , Vec3d{0.,0.,0.} );
      std::vector< Vec3d > rcom_cos( nc , Vec3d{0.,0.,0.} );
      std::vector< Vec3d > rcom_sin( nc , Vec3d{0.,0.,0.} );
      spin_mutex_array chain_lock;
      chain_lock.resize( nc );

      double xmin = domain->bounds().bmin.x;
      double xmax = domain->bounds().bmax.x;
      double ymin = domain->bounds().bmin.y;
      double ymax = domain->bounds().bmax.y;
      double zmin = domain->bounds().bmin.z;
      double zmax = domain->bounds().bmax.z;
      
//      ldbg << "Domain = "<< *domain << std::endl;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc)
        {
          size_t i = grid_ijk_to_index(dims,loc+gl);
          uint64_t const * __restrict__ ids = cells[i][field::id]; ONIKA_ASSUME_ALIGNED(ids);
          double const * __restrict__ rx = cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
          double const * __restrict__ ry = cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
          double const * __restrict__ rz = cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
          double xsi, zeta, theta;
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            uint64_t id = ids[j];
            size_t index_in_chain = id % n_beads;
            size_t chain_index = id / n_beads;
            Vec3d r { rx[j] , ry[j] , rz[j] };
            
//#           pragma omp critical(debug_sl_analytics)
//            ldbg << "id="<<id<<" r="<<r<<std::endl;
            
            chain_lock[chain_index].lock();

            theta = (rx[j]-xmin)/(xmax-xmin)*2.0 *M_PI;
            xsi  = cos(theta);
            zeta = sin(theta);
            rcom_cos[chain_index].x += xsi;
            rcom_sin[chain_index].x += zeta;

            theta = (ry[j]-ymin)/(ymax-ymin)*2.0 *M_PI;
            xsi  = cos(theta);
            zeta = sin(theta);
            rcom_cos[chain_index].y += xsi;
            rcom_sin[chain_index].y += zeta;

            theta = (rz[j]-zmin)/(zmax-zmin)*2.0 *M_PI;
            xsi  = cos(theta);
            zeta = sin(theta);
            rcom_cos[chain_index].z += xsi;
            rcom_sin[chain_index].z += zeta;

            chain_lock[chain_index].unlock();

            if( index_in_chain > 0 )
            {
              uint64_t left_bead = global_to_nearest_local_id( id-1, *id_map, *grid, r, bond_max_search_dist );
              assert( is_particle_id_valid(left_bead) );
              size_t c = 0, p = 0;
              decode_cell_particle( left_bead , c , p );
              assert( grid->is_valid_cell_particle(c,p) );
              Vec3d v = r - Vec3d{ cells[c][field::rx][p] , cells[c][field::ry][p] , cells[c][field::rz][p] };
              
              chain_lock[chain_index].lock();
              ree[chain_index] += v;
              chain_lock[chain_index].unlock();              
            }
          }
        }
        GRID_OMP_FOR_END
      }
      
      MPI_Allreduce(MPI_IN_PLACE,ree.data(),nc*3,MPI_DOUBLE,MPI_SUM,comm);
      MPI_Allreduce(MPI_IN_PLACE,rcom_cos.data(),nc*3,MPI_DOUBLE,MPI_SUM,comm);
      MPI_Allreduce(MPI_IN_PLACE,rcom_sin.data(),nc*3,MPI_DOUBLE,MPI_SUM,comm);

      std::vector< Vec3d > rcom( nc , Vec3d{0.,0.,0.} );
      double theta;
      
      if( rank == 0 )
      {

#       pragma omp parallel for
        for(size_t i=0;i<nc;i++)
        {
          rcom_cos[i].x /= n_beads;
          rcom_sin[i].x /= n_beads;
          theta = atan2(-rcom_sin[i].x, -rcom_cos[i].x) + M_PI;
          rcom[i].x = theta/(2.0*M_PI)*(xmax-xmin)+xmin;

          rcom_cos[i].y /= n_beads;
          rcom_sin[i].y /= n_beads;
          theta = atan2(-rcom_sin[i].y, -rcom_cos[i].y) + M_PI;
          rcom[i].y = theta/(2.0*M_PI)*(ymax-ymin) + ymin;

          rcom_cos[i].z /= n_beads;
          rcom_sin[i].z /= n_beads;
          theta = atan2(-rcom_sin[i].z, -rcom_cos[i].z) + M_PI;
          rcom[i].z = theta/(2.0*M_PI)*(zmax-zmin) + zmin;

        }
        
        if( m_first_run )
        {
          m_ree0 = ree;
          double sum_ree0_norm2 = 0.0;
#         pragma omp parallel for reduction(+:sum_ree0_norm2)
          for(size_t i=0;i<nc;i++)
          {
            sum_ree0_norm2 += norm2( m_ree0[i] );
          }
          m_ree0_norm2 = sum_ree0_norm2 / nc;
          //m_rcom0 = rcom;
          m_prev_rcom = rcom;
          m_displ_rcom.assign( rcom.size() , Vec3d{0.,0.,0.} );
        }

        double ac = 0.0;
        double msd = 0.0;

#       pragma omp parallel for reduction(+:ac,msd)
        for(size_t i=0;i<nc;i++)
        {
          m_displ_rcom[i] += find_periodic_closest_point( rcom[i] , m_prev_rcom[i] , domain->bounds() ) - m_prev_rcom[i];
          m_prev_rcom[i] = rcom[i];

          ac += dot( ree[i] , m_ree0[i] );
          msd += norm2( m_displ_rcom[i] /*rcom[i] - m_rcom0[i]*/ );
        }

        double ree_avg_norm = 0.0;
        for(size_t i=0;i<nc;i++)
        {
          ree_avg_norm += norm( ree[i] );
        }
        ree_avg_norm /= nc;
 
        /*std::cout << std::setprecision(20) << ac
         << " , " << std::setprecision(20) << m_ree0_norm2
         << " , " << std::setprecision(20) << msd << std::endl;*/
         
        msd /= nc;
        ac /= nc;
        ac /= m_ree0_norm2;
        
        //lout << "msd = "<<msd <<", ac="<<ac<<std::endl;
        static const std::string header = "     Step        Time             Regen SL   Imb./Stat.          ACD                 MSD              Volume      ";

        char lb_move_char = ' ';
        if( *move_flag ) { lb_move_char = '*'; }
        std::string lb_value;
        if( *lb_flag )
        {
          if( *lb_inbalance == 0.0 )
          {
            lb_value = "   N/A   ";
          }
          else
          {
            lb_value = format_string("%.3e", *lb_inbalance);
          }
        }

        if( m_first_run && (*print_header) )
        {
          lout << header << std::endl;
        }
        
        double volume = 1.0;
        volume *= bounds_volume( domain->bounds() );
        
        lout<<format_string("%9ld % .8e  %15ld  %c%10s  % .10e  % .10e  % .10e",
          *timestep, *physical_time,
          *sl_regen_count,
          lb_move_char,lb_value,
          ac,
          msd,
          volume)
          <<std::endl;
      }
      m_first_run = false;
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return
R"EOF(
Written by Michelin & CEA/DIF
Credits to Claire Lemarchand, Ioanis Tanis, Thierry Carrard

Compute analytics summury values.
)EOF";
    }

  private:
    std::vector<Vec3d> m_ree0;
    //std::vector<Vec3d> m_rcom0;
    std::vector<Vec3d> m_displ_rcom;
    std::vector<Vec3d> m_prev_rcom;
    double m_ree0_norm2 = 0.0;
    bool m_first_run = true;
  };

  template<class GridT> using SlipLinkAnalyticsOperatorTmpl = SlipLinkAnalyticsOperator<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_analytics", make_grid_variant_operator< SlipLinkAnalyticsOperatorTmpl > );
  }

}


