#include <memory>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/quantity.h>
#include <exanb/core/physics_constants.h>
#include <exanb/core/unityConverterHelper.h>
#include <onika/memory/allocator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exaStamp/ttm/source_term.h>

#include <mpi.h>
#include <iomanip>

namespace exaStamp
{

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz , field::_fx, field::_fy, field::_fz >
    >
  class IonicElectronicHeatTransfer : public OperatorNode
  {
    ADD_SLOT( MPI_Comm       , mpi          , INPUT , MPI_COMM_WORLD );
    
    ADD_SLOT( GridT          , grid         , INPUT , REQUIRED );
    ADD_SLOT( Domain         , domain       , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies, species      , INPUT , REQUIRED );
    ADD_SLOT( double         , rcut_max     , INPUT_OUTPUT , 0.0 );  // neighborhood distance, in grid space
    
    ADD_SLOT( double         , dt           , INPUT , REQUIRED );
    ADD_SLOT( double         , physical_time, INPUT , REQUIRED );

    ADD_SLOT( std::shared_ptr<ScalarSourceTerm> , te_source  , INPUT , std::make_shared<ScalarSourceTerm>() );
    ADD_SLOT( std::shared_ptr<ScalarSourceTerm> , ti_source  , INPUT , std::make_shared<ScalarSourceTerm>() );
    ADD_SLOT( double         , g            , INPUT , 0.0 );
    ADD_SLOT( double         , Ke           , INPUT , 1.0 );
    ADD_SLOT( double         , splat_size   , INPUT , 1.0 );

    ADD_SLOT( bool           , copy_ti_te   , INPUT, false );
    ADD_SLOT( long           , grid_subdiv  , INPUT , 3 );
    ADD_SLOT( GridCellValues , grid_cell_values      , INPUT_OUTPUT );
    ADD_SLOT( double         , electronic_energy , INPUT_OUTPUT , 0.0 );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {
      // compile time constant indicating if grid has type field
      using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
      static constexpr has_type_field_t has_type_field{};

      static constexpr double weight_sum_epsilon = 1.e-13; if constexpr (weight_sum_epsilon==0.0){}
      static constexpr double te_deviation_epsilon = 1.e-13;

      const double cell_size = domain->cell_size();
      const ssize_t subdiv = *grid_subdiv;
      const ssize_t n_subcells = subdiv * subdiv * subdiv;
      const double subcell_size = cell_size / subdiv;
      const double subcell_volume = subcell_size * subcell_size * subcell_size;

      // retreive field data accessor. create data field if needed
      if( ! grid_cell_values->has_field("te") )
      {
        grid_cell_values->add_field("te",subdiv,1);
      }
      assert( size_t(subdiv) == grid_cell_values->field("te").m_subdiv );
      assert( size_t(subdiv * subdiv * subdiv) == grid_cell_values->field("te").m_components );
      auto cell_te_data = grid_cell_values->field_data("te");
      
      // particle splatting size
      const double sp_size = *splat_size;
      // const double sp_volume = sp_size * sp_size * sp_size;

      if( sp_size > subcell_size )
      {
        lerr << "in " << pathname() << std::endl
             << "splat_size = "<<sp_size << " is larger than subcell size = "<<subcell_size<< std::endl
             << "Choose smaller splat_size or coarser grid_subdiv value" << std::endl;
        std::abort();
      }

      const double rcut_max_grid = (*rcut_max) / domain->xform_min_scale();
      if( subcell_size+sp_size/2.0 > rcut_max_grid )
      {
        ldbg << "in " << pathname() << std::endl
             << "subcell_size+sp_size/2.0 = "<<subcell_size<<'+'<<sp_size<<"/2.0 = "<<(subcell_size+sp_size/2.0) << std::endl
             << "is larger than" << std::endl
             << "rcut_max/min_scale = "<< *(rcut_max) <<'/'<<domain->xform_min_scale()<<" = "<<rcut_max_grid<<std::endl
             << "adjust rcut_max "<<rcut_max_grid<<" -> "<< (subcell_size+sp_size/2.0) * domain->xform_min_scale() <<std::endl ;
        *rcut_max = std::max( *rcut_max , (subcell_size+sp_size/2.0) * domain->xform_min_scale() );
      }

      if( grid->number_of_cells() == 0 )
      {
        return;
      }
  
//      static const double k = UnityConverterHelper::convert(legacy_constant::boltzmann, "J/K");
      //ldbg << "cell_heat: dt="<<(*dt)<<std::endl;

      size_t nSpecies = species->size();
      if( nSpecies != 1 && !has_type_field )
      {
        lerr << pathname() << std::endl;
        lerr << "no type information, can't retreive masses" << std::endl;
        std::abort();
      }
      double masses[nSpecies];
      for(size_t i=0;i<nSpecies;i++) { masses[i] = species->at(i).m_mass; }

      auto cells = grid->cells();
      const ssize_t n_cells = grid->number_of_cells();
      const IJK dims = grid->dimension();
      const ssize_t gl = grid->ghost_layers();      

      // Temporary storage for per cell Ti and per cell Laplacian(Te)
      std::vector<double> tmp_storage_ti_LapTe( n_cells * n_subcells * 2 , 0.0 );
      double* Ti = tmp_storage_ti_LapTe.data();
      double* cell_L_Te = tmp_storage_ti_LapTe.data() + n_cells * n_subcells;

      // 1. computes per cell Ti
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,i,cell_loc, schedule(dynamic) )
        {
	      const Vec3d cell_origin = grid->cell_position( cell_loc );
	        
	      GridFieldSetPointerTuple< GridT, FieldSet<field::_rx,field::_ry,field::_rz,field::_vx,field::_vy,field::_vz> > ptrs;
	      cells[i].capture_pointers( ptrs );

          const auto* __restrict__ rx = ptrs[field::rx];
          const auto* __restrict__ ry = ptrs[field::ry];
          const auto* __restrict__ rz = ptrs[field::rz];

          const auto* __restrict__ vx = ptrs[field::vx];
          const auto* __restrict__ vy = ptrs[field::vy];
          const auto* __restrict__ vz = ptrs[field::vz];

          const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);

          const unsigned int n = cells[i].size();
          for(unsigned int j=0;j<n;j++)
          {
            Vec3d r { rx[j] , ry[j] , rz[j] };
            Vec3d v { vx[j] , vy[j] , vz[j] };

            const double mass = get_mass( j, atom_type, masses, has_type_field );
            const double v2 = norm2(v);
 
            IJK center_cell_loc;
            IJK center_subcell_loc;
            Vec3d rco = r - cell_origin;
            localize_subcell( rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc );
            center_cell_loc += cell_loc;

            [[maybe_unused]] double sum_w = 0.0;
            [[maybe_unused]] int nb_contribs = 0;
            [[maybe_unused]] int nb_neighbors = 0;

            for(int ck=-1;ck<=1;ck++)
            for(int cj=-1;cj<=1;cj++)
            for(int ci=-1;ci<=1;ci++)
            {
              ++ nb_neighbors;
              IJK nbh_cell_loc;
              IJK nbh_subcell_loc;
              subcell_neighbor( center_cell_loc, center_subcell_loc, subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc );
              if( grid->contains(nbh_cell_loc) )
              {
                ++ nb_contribs;
                ssize_t nbh_cell_i = grid_ijk_to_index( dims , nbh_cell_loc );
                ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
                assert( nbh_cell_i>=0 && nbh_cell_i<n_cells );
                assert( nbh_subcell_i>=0 && nbh_subcell_i<n_subcells );

                // compute weighted contribution of particle to sub cell
                Vec3d nbh_cell_origin = grid->cell_position(nbh_cell_loc);
                AABB subcell_box = { nbh_cell_origin + nbh_subcell_loc*subcell_size , nbh_cell_origin + (nbh_subcell_loc+1)*subcell_size };
                const double w = particle_smoothing(r, sp_size, subcell_box);
                
                sum_w += w;
                size_t scindex = nbh_cell_i * n_subcells + nbh_subcell_i;
                double Ti_contrib = ( mass * w * v2 ) / subcell_volume;
                
#               pragma omp atomic update
                Ti[ scindex ] += Ti_contrib;
              }

            }

            // particle's contribution accurately distributed
            assert( (nb_contribs<nb_neighbors) || std::fabs(sum_w-1.0) < weight_sum_epsilon );
          }

        }
        GRID_OMP_FOR_END
      }

      // if a simple copy Te <- Ti is requested, stop here
      if( *copy_ti_te )
      {
        ldbg << "copy_ti_te : dims=" <<dims<< std::endl;
        for(ssize_t i=0;i<n_cells;i++)
        {
          for(ssize_t j=0;j<n_subcells; j++)
          {
            cell_te_data.m_data_ptr[ i*cell_te_data.m_stride + j ] = Ti[i*n_subcells+j];
          }
        }
        return;
      }

      // Laplacian computation of Te will need some ghost
      if( gl <= 0 )
      {
        lerr << pathname() << std::endl;
        lerr << "No ghost layers, can't continue" << std::endl;
        std::abort();
      }
            
      // 2) Te<->Te and Te<->Ti heat transfer
      //const Vec3d grid_origin = grid->grid_bounds().bmin;
      const Mat3d xform = domain->xform();
      const double coupling_g = *g;
      const double Te_cond = *Ke;
      const double delta_t = *dt;
      auto* Te = cell_te_data.m_data_ptr;
      
      ldbg << "Ke="<<Te_cond<<", g="<<coupling_g<< std::endl;

      // inspired from https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
      static constexpr double Lap27Norm = 26.0;
      static constexpr double Lap27_compact[4] = { -88/Lap27Norm , 6/Lap27Norm, 3/Lap27Norm, 2/Lap27Norm };

#     ifndef NDEBUG
      static constexpr double Lap27 [3][3][3] = {
         { { 2, 3, 2 } ,
           { 3, 6, 3 } ,
           { 2, 3, 2 } } ,
         { { 3, 6, 3 } ,
           { 6,-88,6 } ,
           { 3, 6, 3 } } ,
         { { 2, 3, 2 } ,
           { 3, 6, 3 } ,
           { 2, 3, 2 } } 
        };
#     endif

      double sum_Te = 0.0;
      double norm_dTe = 0.0;
      double entropy_Te = 0.0;

      // compute Te laplace operator
#     pragma omp parallel num_threads(1)
      {
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:sum_Te,norm_dTe,entropy_Te) )
        {
          //const IJK cell_loc = loc + gl;
          //const size_t cell_i = grid_ijk_to_index( dims , cell_loc );

          for(int ck=0;ck<subdiv;ck++)
          for(int cj=0;cj<subdiv;cj++)
          for(int ci=0;ci<subdiv;ci++)
          {
            IJK sc { ci, cj, ck };
            size_t j = cell_i*n_subcells +  grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , sc );

            // Te discrete Laplacian operator
            double L_Te = 0.0;
            for(int nk=-1;nk<=1;nk++)
            for(int nj=-1;nj<=1;nj++)
            for(int ni=-1;ni<=1;ni++)
            {
              IJK nbh { ni, nj, nk };
              IJK nbh_cell_loc;
              IJK nbh_subcell_loc;
              subcell_neighbor( cell_loc, sc, subdiv, nbh, nbh_cell_loc, nbh_subcell_loc );
              if( grid->contains(nbh_cell_loc) )
              {
                ssize_t nbh_cell_i = grid_ijk_to_index( dims , nbh_cell_loc );
                ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
                assert( nbh_cell_i>=0 && nbh_cell_i<n_cells );
                assert( nbh_subcell_i>=0 && nbh_subcell_i<n_subcells );
                size_t nbh_j = nbh_cell_i*n_subcells + nbh_subcell_i;

                int lap_compact_index = std::abs(ni) + std::abs(nj) + std::abs(nk);
                assert( Lap27_compact[lap_compact_index] == Lap27[ni+1][nj+1][nk+1]/Lap27Norm );
//              L_Te += Te[nbh_j] * Lap27[ni+1][nj+1][nk+1] / Lap27Norm;
                L_Te += Te[nbh_j] * Lap27_compact[lap_compact_index];
              }
            }
            cell_L_Te[j] = L_Te;
            if( ! grid->is_ghost_cell(cell_loc) )
            {
              sum_Te += Te[j];
              const double Ce_Te = Ce(Te[j]);
              const double dTe = ( (Te_cond*L_Te) / Ce_Te ) * delta_t;
              norm_dTe += std::fabs( dTe );
              entropy_Te += Te[j] * std::log(Te[j]) * subcell_volume ;
            }
          }
        }
        GRID_OMP_FOR_END
      }

      // check global absolute electronic energy variation (from dissipation)
//#     ifndef NDEBUG
      {
        double tmp[3] = { sum_Te, norm_dTe, entropy_Te };
        MPI_Allreduce(MPI_IN_PLACE,tmp,3,MPI_DOUBLE,MPI_SUM,*mpi);
        sum_Te = tmp[0];
        norm_dTe = tmp[1];
        entropy_Te = tmp[2];
        ldbg <<"sum_Te="<<sum_Te <<" norm_dTe/sum_Te="<<norm_dTe/sum_Te<<" , entropy_Te="<<entropy_Te << std::endl;
      }
//#     endif


      // 3. Compute Te dissipation, Te<->Ti transfer & source terms
      double old_sum_Te = sum_Te;
      double sum_dTe = 0.0;
      double sum_Se = 0.0;
      double sum_Si = 0.0;
      sum_Te = 0.0;
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(static) reduction(+:sum_Te,sum_dTe,sum_Se,sum_Si) )
        {
          //const IJK cell_loc = loc + gl;
          //const size_t cell_i = grid_ijk_to_index( dims , cell_loc );
          const Vec3d cell_origin = grid->cell_position( cell_loc );

          for(int ck=0;ck<subdiv;ck++)
          for(int cj=0;cj<subdiv;cj++)
          for(int ci=0;ci<subdiv;ci++)
          {
            IJK sc { ci, cj, ck };
            Vec3d scr { ci+0.5, cj+0.5, ck+0.5 };
            const size_t scindex = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , sc );
            const size_t idx_ti = cell_i*n_subcells + scindex ;
            const size_t idx_te = cell_i*cell_te_data.m_stride + scindex;

            // source terms
            const Vec3d center = xform * ( cell_origin + scr*subcell_size );
            const double Si = (*ti_source)->S( center, *physical_time );
            const double Se = (*te_source)->S( center, *physical_time );
            
            // sum external contributions (source terms)
            sum_Se += Se;
            sum_Si += Si;
            
            // Te/Ti coupling
            const double g_Te_Ti = coupling_g * ( Te[idx_te] - Ti[idx_ti] );
            
            // Te dependency
            const double Ce_Te = Ce(Te[idx_te]);
            
            // cell temperature increments
            const double dTe = ( Te_cond*cell_L_Te[idx_ti] - g_Te_Ti + Se ) / Ce_Te;
            Te[idx_te] += dTe * delta_t;
            if( ! grid->is_ghost_cell(cell_loc) )
            {
              sum_Te += Te[idx_te];
              sum_dTe += dTe * delta_t;
            }

            //assert( Si == 0.0 );

            // as an output, Ti array stores Xi
            double ksi = 0.0;
            if( Ti[idx_ti] > 0.0 )
            {
              // ksi = ( ( g_Te_Ti + Si ) * subcell_volume ) / ( Ti[j] * subcell_volume );
              ksi = ( g_Te_Ti + Si ) / Ti[idx_ti];
            }
            Ti[idx_ti] = ksi;

            if(coupling_g==0.0) { assert( ksi==0.0 ); }
          }
        }
        GRID_OMP_FOR_END
      }

      {
        double tmp[4] = { sum_Te, sum_dTe, sum_Se, sum_Si };
        MPI_Allreduce(MPI_IN_PLACE,tmp,4,MPI_DOUBLE,MPI_SUM,*mpi);
        sum_Te = tmp[0];
        sum_dTe = tmp[1];
        sum_Se = tmp[2];
        sum_Si = tmp[3];
      }
      
      if(sum_Te > 0.)
      {
        double te_dev = (sum_Te - old_sum_Te - sum_Se) / sum_Te;
        ldbg << "Te dev="<< te_dev << " sum_dTe="<<sum_dTe << " sum_Se="<<sum_Se<<" sum_Si="<<sum_Si<< std::endl;
        if( te_dev > te_deviation_epsilon )
        {
          lerr << "Te deviation too big : "<<old_sum_Te<<" -> "<<sum_Te<<" , dev="<<te_dev<<std::endl;
        }
        // else { lout << "Te deviation Ok : "<<old_sum_Te<<" -> "<<sum_Te<<" , dev="<<te_dev<<std::endl; }
      }
      *electronic_energy = sum_Te * subcell_volume;

      // 4. projects back Ti variation to particles
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
        {
          const IJK cell_loc = loc + gl;
          const size_t i = grid_ijk_to_index( dims , cell_loc );
	      const Vec3d cell_origin = grid->cell_position( cell_loc );

	      GridFieldSetPointerTuple< GridT, FieldSet<field::_rx,field::_ry,field::_rz,field::_vx,field::_vy,field::_vz,field::_fx,field::_fy,field::_fz> > ptrs;
	      cells[i].capture_pointers( ptrs );

          const auto* __restrict__ rx = ptrs[field::rx];
          const auto* __restrict__ ry = ptrs[field::ry];
          const auto* __restrict__ rz = ptrs[field::rz];

          const auto* __restrict__ vx = ptrs[field::vx];
          const auto* __restrict__ vy = ptrs[field::vy];
          const auto* __restrict__ vz = ptrs[field::vz];

          auto* __restrict__ fx = ptrs[field::fx];
          auto* __restrict__ fy = ptrs[field::fy];
          auto* __restrict__ fz = ptrs[field::fz];

          const auto* __restrict__ atom_type = cells[i].field_pointer_or_null(field::type);

          const unsigned int n = cells[i].size();
          for(unsigned int j=0;j<n;j++)
          {
            const double mass = get_mass( j, atom_type, masses, has_type_field );
            const Vec3d r { rx[j] , ry[j] , rz[j] };
            const Vec3d v { vx[j] , vy[j] , vz[j] };
            Vec3d f { 0. , 0. , 0. };
            
            IJK center_cell_loc;
            IJK center_subcell_loc;
            Vec3d rco = r - cell_origin;
            localize_subcell( rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc );
            center_cell_loc += cell_loc;
  
            for(int ck=-1;ck<=1;ck++)
            for(int cj=-1;cj<=1;cj++)
            for(int ci=-1;ci<=1;ci++)
            {
              IJK nbh_cell_loc;
              IJK nbh_subcell_loc;
              subcell_neighbor( center_cell_loc, center_subcell_loc, subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc );
              ssize_t nbh_cell_i = grid_ijk_to_index( dims , nbh_cell_loc );
              ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{subdiv,subdiv,subdiv} , nbh_subcell_loc );
              assert( nbh_cell_i>=0 && nbh_cell_i<n_cells );
              assert( nbh_subcell_i>=0 && nbh_subcell_i<n_subcells );
              
              // compute weighted mass of particle in this sub cell
              Vec3d nbh_cell_origin = grid->cell_position(nbh_cell_loc);
              AABB subcell_box = { nbh_cell_origin + nbh_subcell_loc*subcell_size , nbh_cell_origin + (nbh_subcell_loc+1)*subcell_size };
              const double w = particle_smoothing(r, sp_size, subcell_box);

              const double ksi = Ti[ nbh_cell_i * n_subcells + nbh_subcell_i ];
              f += ksi * mass * w * v;
            }

            if(coupling_g==0.0) { assert( f == (Vec3d{0.,0.,0.}) ); }

            fx[j] += f.x;
            fy[j] += f.y;
            fz[j] += f.z;
          }

        }
        GRID_OMP_FOR_END
      }

    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Handles heat transfer between ionic and electronic temperatures. Electronic temperature (Te) is held by a rectilinear grid,
while ionic temperature (Ti) commes from particles kinetic energy.
1. compute per cell Ti
2. transfer heat between Te and Ti
3. solve heat equation on the rectilinear grid for Te
4. project back Ti to particles through speed adjustment
)EOF";
    }


  private:
    static inline double get_mass( unsigned int j, const uint8_t* type_ptr, const double* masses, std::true_type )
    {
      return masses[type_ptr[j]];
    }

    static inline double get_mass( unsigned int j, const uint8_t* type_ptr, const double* masses, std::false_type )
    {
      return masses[0];
    }

    static inline void localize_subcell( const Vec3d& r, double cell_size, double sub_cellsize, ssize_t subdiv, IJK& cell_loc, IJK& subcell_loc )
    {
      cell_loc = make_ijk( r / cell_size );
      Vec3d ro = r - (cell_loc*cell_size);
      subcell_loc = vclamp( make_ijk(ro / sub_cellsize) , 0 , subdiv-1 );
    }

    static inline void subcell_neighbor( const IJK& cell_loc, const IJK& subcell_loc, ssize_t subdiv, IJK ninc, IJK& nbh_cell_loc, IJK& nbh_subcell_loc )
    {
      nbh_cell_loc = cell_loc;
      nbh_subcell_loc = subcell_loc + ninc;
      if(nbh_subcell_loc.i<0) { -- nbh_cell_loc.i; } else if(nbh_subcell_loc.i>=subdiv) { ++ nbh_cell_loc.i; }
      if(nbh_subcell_loc.j<0) { -- nbh_cell_loc.j; } else if(nbh_subcell_loc.j>=subdiv) { ++ nbh_cell_loc.j; }
      if(nbh_subcell_loc.k<0) { -- nbh_cell_loc.k; } else if(nbh_subcell_loc.k>=subdiv) { ++ nbh_cell_loc.k; }
      nbh_subcell_loc.i = ( nbh_subcell_loc.i + subdiv ) % subdiv;
      nbh_subcell_loc.j = ( nbh_subcell_loc.j + subdiv ) % subdiv;
      nbh_subcell_loc.k = ( nbh_subcell_loc.k + subdiv ) % subdiv;      
    }

    static inline double Ce( double Te )
    {
      return 1.0;
    }

    // @return how much of this particle contributes to region cell_box.
    // sum of contributions for all disjoint cell_box paving the domain is guaranteed to be 1.0
    static inline double particle_smoothing(const Vec3d& r, double sp_size, const AABB& cell_box)
    {
      AABB contrib_box = { r - sp_size*0.5 , r + sp_size*0.5 };
      AABB sub_contrib_box = intersection( contrib_box , cell_box );
      double w = 0.0;
      if( ! is_nil(sub_contrib_box) ) { w = bounds_volume(sub_contrib_box) / (sp_size*sp_size*sp_size); }
      assert( w>=0. && w<=(1.0+1.e-13) );
      return w;
    }

  };

  template<class GridT> using IonicElectronicHeatTransferTmpl = IonicElectronicHeatTransfer<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("ionic_eletronic_heat_transfer", make_grid_variant_operator< IonicElectronicHeatTransferTmpl > );
  }

}
