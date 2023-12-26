#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/print_particle.h>

#include <limits>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  struct CellParticleMinMaxFinder
  {
    const char * m_name = nullptr;
    std::function<double(size_t,size_t)> m_func = nullptr;
    ssize_t m_min_cell = -1;
    ssize_t m_max_cell = -1;
    ssize_t m_min_part = -1;
    ssize_t m_max_part = -1;
    double m_min = std::numeric_limits<double>::max();
    double m_max = std::numeric_limits<double>::lowest();
    int m_min_rank = -1;
    int m_max_rank = -1;
    
    inline void update(size_t cell, size_t part)
    {
      if( m_func != nullptr )
      {
        double x = m_func( cell, part );
        if( m_min_cell==-1 || m_min_part==-1 || m_max_cell==-1 || m_max_part==-1 )
        {
          m_min = m_max = x;
          m_min_cell = m_max_cell = cell;
          m_min_part = m_max_part = part;
        }
        else
        {
          if( x < m_min )
          {
            m_min_cell = cell; m_min_part = part; m_min = x;
          }
          if( x > m_max )
          {
            m_max_cell = cell; m_max_part = part; m_max = x;
          }
        }
      }
    }
  };

  // =====================================================================
  // ========================== Min max debugger ========================
  // =====================================================================

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx , field::_fy , field::_fz , field::_vx , field::_vy , field::_vz >
    >
  class DebugMinMaxParticle : public OperatorNode
  {
    ADD_SLOT( MPI_Comm           , mpi             , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( bool , abort_on_out_of_range , false, INPUT);
    ADD_SLOT( GridT , grid , INPUT);

  public:
    void execute() override final
    {
      auto cells = grid->cells();
      size_t n_cells = grid->number_of_cells();

      std::vector<CellParticleMinMaxFinder> min_max_finders = {
          { "rx" , [&cells](size_t c, size_t p) { return cells[c][field::rx][p]; } }
        , { "ry" , [&cells](size_t c, size_t p) { return cells[c][field::ry][p]; } }
        , { "rz" , [&cells](size_t c, size_t p) { return cells[c][field::rz][p]; } }
        , { "velocity" , [&cells](size_t c, size_t p) { return norm( Vec3d{cells[c][field::vx][p],cells[c][field::vy][p],cells[c][field::vz][p]} ); } }
        , { "force" , [&cells](size_t c, size_t p) { return norm( Vec3d{cells[c][field::fx][p],cells[c][field::fy][p],cells[c][field::fz][p]} ); } }
        };

      unsigned long long total_particles = 0;
      for(size_t c=0;c<n_cells;c++)
      {
        if(! grid->is_ghost_cell(c) )
        {
          size_t n_particles = cells[c].size();
          for(size_t p=0;p<n_particles;p++)
          {
            for(auto & f:min_max_finders) f.update(c,p);
            ++ total_particles;
          }
        }
      }

      int rank = 0;
      int np = 1;
      MPI_Comm_rank(*mpi,&rank);
      MPI_Comm_size(*mpi,&np);

      MPI_Allreduce(MPI_IN_PLACE,&total_particles,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,*mpi);
      
      for(auto & f:min_max_finders)
      {
        double all_min = f.m_min;
        MPI_Allreduce(MPI_IN_PLACE,&all_min,1,MPI_DOUBLE,MPI_MIN,*mpi);
        if(all_min != f.m_min) 
        {
          f.m_min_cell = -1;
          f.m_min_part = -1;
          f.m_min_rank = np;
        }
        else
        {
          f.m_min_rank = rank;
        }
        f.m_min = all_min;
        MPI_Allreduce(MPI_IN_PLACE,&(f.m_min_rank),1,MPI_INT,MPI_MIN,*mpi);

        double all_max = f.m_max;
        MPI_Allreduce(MPI_IN_PLACE,&all_max,1,MPI_DOUBLE,MPI_MAX,*mpi);
        if(all_max != f.m_max)
        {
          f.m_max_cell = -1;
          f.m_max_part = -1;
          f.m_max_rank = np;
        }
        else
        {
          f.m_max_rank = rank;
        }
        f.m_max = all_max;
        MPI_Allreduce(MPI_IN_PLACE,&(f.m_max_rank),1,MPI_INT,MPI_MIN,*mpi);
      }

      lout << "total particles = "<<total_particles<<std::endl;
      for(const auto & f:min_max_finders)
      {
        if( ( f.m_min_cell>=0 && f.m_min_part>=0 && f.m_min_rank == rank ) || ( f.m_max_cell>=0 && f.m_max_part>=0 && f.m_max_rank == rank ) )
        {
          std::ostringstream oss;
          oss<< default_stream_format;
          if( f.m_min_cell>=0 && f.m_min_part>=0 && f.m_min_rank == rank )
          {
            oss << f.m_name << " MIN : P"<<rank<<", cell="<<f.m_min_cell<<", part="<<f.m_min_part<<" : ";
            print_particle( oss , cells[f.m_min_cell][f.m_min_part] );
//            oss << "\n";
          }
          if( f.m_max_cell>=0 && f.m_max_part>=0  && f.m_max_rank == rank )
          {
            oss << f.m_name << " MAX : P"<<rank<<", cell="<<f.m_max_cell<<", part="<<f.m_max_part<<" : ";
            print_particle( oss , cells[f.m_max_cell][f.m_max_part] );
//            oss << "\n";
          }
          std::cout << oss.str() << std::flush;
        }
        if( *abort_on_out_of_range && ( std::isnan(f.m_min) || std::isinf(f.m_min) || std::isnan(f.m_max) || std::isinf(f.m_max) ) )
        {
          std::abort();
        }
      }
    }
  };

  template<class GridT> using DebugMinMaxParticleTmpl = DebugMinMaxParticle<GridT>;
  
  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "debug_min_max_particle", make_grid_variant_operator< DebugMinMaxParticleTmpl > );
  }

}
