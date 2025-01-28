#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/memory/allocator.h> // for ONIKA_ASSUME_ALIGNED macro
#include <exanb/compute/compute_pair_optional_args.h>
#include <vector>
#include <iomanip>

#include <mpi.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz >
    >
  class WallOperator : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi      , INPUT , MPI_COMM_WORLD);
    ADD_SLOT( bool     , collective_stats, INPUT , false );    
  
    ADD_SLOT( GridT  , grid    , INPUT_OUTPUT );
    ADD_SLOT( Vec3d  , normal  , INPUT , Vec3d{1.0,0.0,0.0} );
    ADD_SLOT( double , offset  , INPUT , 0.0 );
    ADD_SLOT( double , cutoff  , INPUT , 5.0 );
    ADD_SLOT( Domain , domain  , INPUT , REQUIRED );
    ADD_SLOT( long   , exponent, INPUT , 12 );
    ADD_SLOT( double , epsilon , INPUT , make_quantity(1.0e-19,"J").convert() );

    using has_ep_field_t = typename GridT::CellParticles::template HasField < field::_ep > ;
    static constexpr bool has_ep_field = has_ep_field_t::value;

  public:

    inline void execute () override final
    {
        ldbg<<"WallOperator: N="<<(*normal)<<", D="<<(*offset)<<", R="<<(*cutoff)<<", E="<<*epsilon <<std::endl;
//std::cout<<std::endl<<"WallOperator: N="<<(*normal)<<", D="<<(*offset)<<", R="<<(*cutoff)<<", E="<<*epsilon <<std::endl;
        if( ! domain->xform_is_identity() )
          {
            LinearXForm cp_xform { domain->xform() };
            apply_wall( *grid, *normal, - (*offset), *cutoff, *exponent, *epsilon, cp_xform );
          }
        else
          {
            NullXForm cp_xform{};
            apply_wall( *grid, *normal, - (*offset), *cutoff, *exponent, *epsilon, cp_xform );
          }
    }

    inline void yaml_initialize(const YAML::Node& node) override final
    {
      YAML::Node tmp;
      if( node.IsSequence() && node.size()==3 )
      {
        tmp["normal"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize( tmp );
    }

  private:
  
    template<class OptionalXformT>
    inline void apply_wall(GridT& grid, const Vec3d& N, double D, double R, int exposant, double epsilon, const OptionalXformT& xform )
    {

      // D : -offset
      // R : cutoff

      //std::cout << "In src/compute/wall.cpp entree dans apply_wall" << std::endl << std::flush;
      const double cell_size = grid.cell_size();
      const Vec3d grid_origin = grid.origin() + ( grid.offset() * cell_size );
            
      auto cells = grid.cells();
      IJK dims = grid.dimension();
      ssize_t gl = grid.ghost_layers();

      //double PosMur=-D;
      //lout<<"PosMur="<<PosMur<<std::endl << std::flush ;

      unsigned long n_wall_particles = 0; // number of particle interacting with the wall
      double minwalld = 2*R;

//double rxmin=1e10;
//double rxmax=-1e10;

      //std::cout << "In src/compute/wall.cpp apply_wall avant boucle omp for" << std::endl << std::flush;
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,gloc, schedule(dynamic) reduction(+:n_wall_particles) reduction(min:minwalld) )
//GRID_OMP_FOR_BEGIN(dims-2*gl,_,gloc, schedule(dynamic) reduction(+:n_wall_particles) reduction(min:minwalld) reduction(min:rxmin) reduction(max:rxmax))
        {
          IJK loc = gloc + gl;
          size_t cell_i = grid_ijk_to_index( dims , loc );

          Vec3d cell_corner[8] =
            { grid_origin + (loc+IJK{0,0,0}) *cell_size
            , grid_origin + (loc+IJK{0,0,1}) *cell_size
            , grid_origin + (loc+IJK{0,1,0}) *cell_size
            , grid_origin + (loc+IJK{0,1,1}) *cell_size
            , grid_origin + (loc+IJK{1,0,0}) *cell_size
            , grid_origin + (loc+IJK{1,0,1}) *cell_size
            , grid_origin + (loc+IJK{1,1,0}) *cell_size
            , grid_origin + (loc+IJK{1,1,1}) *cell_size };

#         pragma omp simd
          for(int j=0;j<8;j++) { cell_corner[j] = xform.transformCoord(cell_corner[j]); }

          double mind = dot( cell_corner[0] , N ) + D;
          double maxd = mind;

#         pragma omp simd reduction(min:mind) reduction(max:maxd)
          for(int j=0;j<8;j++)
          {
            double d = dot( cell_corner[j] , N ) + D;
            mind = std::min( mind , d );
            maxd = std::max( maxd , d );
          }
          
          if( mind < R && maxd > (-R) )
          {  
            const double* __restrict__ rx = cells[cell_i][field::rx];
            const double* __restrict__ ry = cells[cell_i][field::ry];
            const double* __restrict__ rz = cells[cell_i][field::rz];
            double* __restrict__ fx = cells[cell_i][field::fx];
            double* __restrict__ fy = cells[cell_i][field::fy];
            double* __restrict__ fz = cells[cell_i][field::fz];
            
            double* __restrict__ ep = nullptr;
            if constexpr (has_ep_field) { ep = cells[cell_i][field::ep]; }

#           ifndef NDEBUG
            const uint64_t* __restrict__ id = cells[cell_i].field_pointer_or_null(field::id);
#           endif

            ONIKA_ASSUME_ALIGNED( rx );
            ONIKA_ASSUME_ALIGNED( ry );
            ONIKA_ASSUME_ALIGNED( rz );
            ONIKA_ASSUME_ALIGNED( fx );
            ONIKA_ASSUME_ALIGNED( fy );
            ONIKA_ASSUME_ALIGNED( fz );
            ONIKA_ASSUME_ALIGNED( ep );

            size_t n = cells[cell_i].size();
                       
#           ifdef NDEBUG
#           pragma omp simd reduction(+:n_wall_particles) reduction(min:minwalld)
#           endif
            for(size_t j=0;j<n;j++)
            {
              Vec3d r { rx[j], ry[j], rz[j] };
              r = xform.transformCoord(r);
              //rxmin=std::min(rxmin,rx[j]);
              //rxmax=std::max(rxmax,rx[j]);
              
              // unsigned distance => bidirectional wall
              double d_sign = dot( r , N ) + D;
              double d = std::abs( d_sign );

              // numerical issues arises if too close to wall
#             ifndef NDEBUG
              if( d < R*1.e-3 )
              {
#               pragma omp critical
                {
                  // lerr <<"in "<< pathname() << std::endl;
                  lerr<<"part #"<<id[j]<<": r="<<rx[j]<<", rWall="<<D<<", d="<<d<<", maxRcut="<<R<<std::endl;
                  lerr<<"particle to close to wall"<<std::endl;
                  std::abort();
                }
              }
#             endif


              if( d <= R )
              {
                ++ n_wall_particles;
                minwalld = std::min( minwalld , d );
                //std::cout << "In src/compute/wall.cpp : distance to wall = " << d << " - " << n_wall_particles <<" particles near wall" <<std::endl;
                
                const double ratio = 1.0 - R / d;
                const double ratio_puis_exposant   = pow(ratio,exposant);

                double f_contrib =  -epsilon * exposant * ( R / (d*d) ) * ratio_puis_exposant / ratio ;
                Vec3d F = N * f_contrib;
                if (d_sign<0.) F = -F;
                //std::cout << "F_wall_d= " << F << std::endl;

                // Energy
                if constexpr (has_ep_field) { ep[j] += epsilon * ratio_puis_exposant; }

                // Forces
                fx[j] += F.x;
                fy[j] += F.y;
                fz[j] += F.z;
              }
            }
          }
        }
        GRID_OMP_FOR_END
      }
      //std::cout << "In src/compute/wall.cpp apply_wall apres boucle omp for" << std::endl << std::flush;

      if( *collective_stats )
      {
        MPI_Allreduce(MPI_IN_PLACE,&n_wall_particles,1,MPI_UNSIGNED_LONG,MPI_SUM,*mpi);
        MPI_Allreduce(MPI_IN_PLACE,&minwalld,1,MPI_DOUBLE,MPI_MIN,*mpi);
//MPI_Allreduce(MPI_IN_PLACE,&rxmin,1,MPI_DOUBLE,MPI_MIN,*mpi);
//MPI_Allreduce(MPI_IN_PLACE,&rxmax,1,MPI_DOUBLE,MPI_MAX,*mpi);
      }

//      std::cout << "In src/compute/wall.cpp : " <<n_wall_particles<<" particles near wall, min distance = "<<minwalld <<std::endl;
//std::cout << "In src/compute/wall.cpp : " <<n_wall_particles<<" particles near wall, min distance = "<<minwalld <<std::endl<<std::endl;
//std::cout << std::endl << "x_min= " << rxmin << " -> d_wall= " << rxmin+D << std::endl << "x_max= " << rxmax << " -> d_wall= " << rxmax+D << " D= " << D << std::endl << std::endl;
//if (abs(rxmin+D)<abs(rxmax+D)) {
//	std::cout << "x_min= " << rxmin << " -> d_wall= " << rxmin+D << std::endl;
//} else {
//	std::cout << "x_max= " << rxmax << " -> d_wall= " << -(rxmax+D) << std::endl;
//}

      if( n_wall_particles > 0 )
      {
        ldbg<<n_wall_particles<<" particles near wall, min distance = "<<minwalld <<std::endl;
      }
    }
  };
  
  // this helps older versions of gcc handle the unnamed default second template parameter
  template <class GridT> using WallOperatorTemplate = WallOperator<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(wall)
  {
    OperatorNodeFactory::instance()->register_factory( "wall", make_grid_variant_operator< WallOperatorTemplate > );
  }

}

