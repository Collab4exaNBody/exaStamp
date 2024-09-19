//#include <chrono>
#include <memory>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/molecule/id_map.h>
#include <exaStamp/molecule/periodic_r_delta.h>

#include <mpi.h>
#include <exanb/mpi/grid_update_ghosts.h>

namespace exaStamp
{

  template<typename GridT
    , class = AssertGridHasFields< GridT, field::_id >
    >
  class MoleculeUnfold : public OperatorNode
  {
    using UpdateGhostsScratch = typename UpdateGhostsUtils::UpdateGhostsScratch;

    ADD_SLOT( MPI_Comm                 , mpi               , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( long                     , mpi_tag           , INPUT , 0 );

    ADD_SLOT( Domain                   , domain            , INPUT );
    ADD_SLOT( GridT       , grid          , INPUT_OUTPUT );

    ADD_SLOT( IdMap       , id_map        , INPUT );
    ADD_SLOT( IdMapGhosts , id_map_ghosts , INPUT );

    ADD_SLOT( UpdateGhostsScratch      , ghost_comm_buffers, PRIVATE );
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme , INPUT_OUTPUT , REQUIRED );

  public:
    inline void execute ()  override final
    {
      // clear previous content
      id_map->clear();
      id_map_ghosts->clear();

      auto cells = grid->cells_accessor();
      size_t n_cells = grid->number_of_cells();
      const auto idmol_field = grid->field_accessor( field::idmol );
      const auto cmol_field = grid->field_accessor( field::cmol );
//      static constexpr uint64_t null_id = std::numeric_limits<uint64_t>::max();
//      static constexpr uint64_t null_loc = std::numeric_limits<uint64_t>::max();
      static constexpr size_t null_index = std::numeric_limits<size_t>::max();

      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      {
        size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++) cells[cell_i][idmol_field][p_i] = cells[cell_i][field::id][p_i];
      }

      auto pecfunc = [self=this](auto ... args) { return self->parallel_execution_context(args ...); };
      auto pesfunc = [self=this](unsigned int i) { return self->parallel_execution_stream(i); }; 

      std::unordered_map<uint64_t , uint64_t> molecule_owner;

      bool stable;
      int pass = 0;
      do
      {
        const auto ghost_update_fields = onika::make_flat_tuple( idmol_field ); 
        grid_update_ghosts( ldbg, *mpi, *ghost_comm_scheme, *grid, *domain, nullptr,
                            *ghost_comm_buffers, pecfunc,pesfunc, ghost_update_fields,
                            *mpi_tag, true, true, true, true, false, std::false_type{} );
        stable = true;
        ++ pass;
        ldbg << "molecule_unfold: connect molecules, pass="<<pass<<std::endl;

        for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          const auto * __restrict__ ids = cells[cell_i][field::id];
          const auto *  __restrict__ cmol = cells[cell_i][cmol_field];
          auto * __restrict__ mol_ids   = cells[cell_i][idmol_field];
          size_t n = cells[cell_i].size();
          for(size_t p_i=0;p_i<n;p_i++)
          { 
            uint64_t molidmin = std::min( ids[p_i] , mol_ids[p_i] );
            for(size_t j=0;j<cmol[p_i].size();j++)
            {
              uint64_t locations[32];
              const size_t nloc = all_atoms_from_idmap( cmol[p_i][j] , *id_map , *id_map_ghosts, locations, 32 );
              for(size_t k=0;k<nloc;k++)                
              {
                size_t cell=null_index, pos=null_index;
                decode_cell_particle( locations[k] , cell, pos );
                assert( cell!=null_index && pos!=null_index );
                molidmin = std::min( molidmin ,  cells[cell][idmol_field][pos] );
              }
            }
            
            stable = stable && ( molidmin == mol_ids[p_i] );
            mol_ids[p_i] = molidmin;
            for(size_t j=0;j<cmol[p_i].size();j++)
            {
              uint64_t locations[32];
              const size_t nloc = all_atoms_from_idmap( cmol[p_i][j] , *id_map , *id_map_ghosts, locations, 32 );
              for(size_t k=0;k<nloc;k++)                
              {
                size_t cell=null_index, pos=null_index;
                decode_cell_particle( locations[k] , cell, pos );
                assert( cell!=null_index && pos!=null_index );
                if( cells[cell][idmol_field][pos] != molidmin )
                {
                  cells[cell][idmol_field][pos] = molidmin;
                  stable = false;
                }
              }
            }
            
          }
        }
        
      } while( ! stable );

      const auto rxf_field = grid->field_accessor( field::rxf );
      const auto ryf_field = grid->field_accessor( field::ryf );
      const auto rzf_field = grid->field_accessor( field::rzf );

      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      {
        const auto * __restrict__ ids = cells[cell_i][field::id];
        const auto * __restrict__ mol_ids   = cells[cell_i][idmol_field];
        const auto * __restrict__ rx   = cells[cell_i][field::rx];
        const auto * __restrict__ ry   = cells[cell_i][field::ry];
        const auto * __restrict__ rz   = cells[cell_i][field::rz];
        auto * __restrict__ rxf   = cells[cell_i][rxf_field];
        auto * __restrict__ ryf   = cells[cell_i][ryf_field];
        auto * __restrict__ rzf   = cells[cell_i][rzf_field];
        const size_t n = cells[cell_i].size();
        const bool is_ghost = grid->is_ghost_cell(cell_i);
        for(size_t p_i=0;p_i<n;p_i++)
        {
          if( ids[p_i] == mol_ids[p_i] && !is_ghost )
          {
            const auto ploc = encode_cell_particle(cell_i,p_i,0);
            if( molecule_owner.find(ploc) != molecule_owner.end() )
            {
              fatal_error() << "duplicate owner particle for molecule id #"<<mol_ids[p_i]<<std::endl;
            }
            molecule_owner.insert( { mol_ids[p_i] , ploc } );
            rxf[p_i] = rx[p_i];
            ryf[p_i] = ry[p_i];
            rzf[p_i] = rz[p_i];
          }
          else
          {
            rxf[p_i] = std::numeric_limits<double>::quiet_NaN();
            ryf[p_i] = std::numeric_limits<double>::quiet_NaN();
            rzf[p_i] = std::numeric_limits<double>::quiet_NaN();            
          }
        }
      }

      const Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 

      ldbg << "Number of molecules = "<< molecule_owner.size()<<" , size_box = "<<size_box<<" , half_min_size_box = "<<half_min_size_box<<std::endl;

      pass = 0;
      do
      {
        const auto ghost_update_fields = onika::make_flat_tuple( rxf_field, ryf_field, rzf_field ); 
        grid_update_ghosts( ldbg, *mpi, *ghost_comm_scheme, *grid, *domain, nullptr,
                            *ghost_comm_buffers, pecfunc,pesfunc, ghost_update_fields,
                            *mpi_tag, true, true, true, true, false, std::false_type{} );
        stable = true;
        ++ pass;
        ldbg << "molecule_unfold: unfold positions, pass="<<pass<<std::endl;

        for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          const auto *  __restrict__ cmol = cells[cell_i][cmol_field];
          auto * __restrict__ rxf   = cells[cell_i][rxf_field];
          auto * __restrict__ ryf   = cells[cell_i][ryf_field];
          auto * __restrict__ rzf   = cells[cell_i][rzf_field];

          size_t n = cells[cell_i].size();
          for(size_t p_i=0;p_i<n;p_i++)
          {
            const Vec3d ri = { rxf[p_i], ryf[p_i], rzf[p_i] };
            if( !std::isnan(ri.x) && !std::isnan(ri.y) && !std::isnan(ri.z) )
            {
              for(size_t j=0;j<cmol[p_i].size();j++)
              {
                uint64_t locations[32];
                const size_t nloc = all_atoms_from_idmap( cmol[p_i][j] , *id_map , *id_map_ghosts, locations, 32 );
                for(size_t k=0;k<nloc;k++)                
                {
                  size_t cell=null_index, pos=null_index;
                  decode_cell_particle( locations[k] , cell, pos );
                  assert( cell!=null_index && pos!=null_index );
                  if( std::isnan(cells[cell][rxf_field][pos]) || std::isnan(cells[cell][ryf_field][pos]) || std::isnan(cells[cell][rzf_field][pos]) )
                  {
                    const Vec3d rj = { cells[cell][field::rx][pos], cells[cell][field::ry][pos], cells[cell][field::rz][pos] };
                    const Vec3d rij = periodic_r_delta( ri , rj , size_box , half_min_size_box );
                    if( norm(rij)>=half_min_size_box || rij==Vec3d{0,0,0} )
                    {
                      fatal_error() << "in cell #"<<cell_i<<" , particle #"<<p_i<<" , bond #"<<j<<" bond distance too long : "<<norm(rij)<<" >= "<<half_min_size_box<< std::endl;
                    }
                    const Vec3d rj_unfolded = ri + rij;
                    cells[cell][rxf_field][pos] = rj_unfolded.x;
                    cells[cell][ryf_field][pos] = rj_unfolded.y;
                    cells[cell][rzf_field][pos] = rj_unfolded.z;
                    stable = false;
                  }
                }
              }
              
            }
            
          }
        }
        
      } while( ! stable );

      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      {
        auto * __restrict__ rx   = cells[cell_i][field::rx];
        auto * __restrict__ ry   = cells[cell_i][field::ry];
        auto * __restrict__ rz   = cells[cell_i][field::rz];
        const auto * __restrict__ rxf   = cells[cell_i][rxf_field];
        const auto * __restrict__ ryf   = cells[cell_i][ryf_field];
        const auto * __restrict__ rzf   = cells[cell_i][rzf_field];
        size_t n = cells[cell_i].size();
        for(size_t p_i=0;p_i<n;p_i++)
        {
          if( ! std::isnan(rxf[p_i]) && ! std::isnan(ryf[p_i]) && ! std::isnan(rzf[p_i]) )
          {
            rx[p_i] = rxf[p_i];
            ry[p_i] = ryf[p_i];
            rz[p_i] = rzf[p_i];
          }
          else
          {
            fatal_error() << "in cell #"<<cell_i<<" , particle #"<<p_i<<" not unfolded"<<std::endl;
          }
        }
      }

      grid->remove_flat_array( "rxf" );
      grid->remove_flat_array( "ryf" );
      grid->remove_flat_array( "rzf" );
    }

  };

  template<class GridT> using MoleculeUnfoldTmpl = MoleculeUnfold<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "molecule_unfold", make_grid_variant_operator<MoleculeUnfoldTmpl> );
  }

}
