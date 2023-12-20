#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/domain.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <mpi.h>
#include <vector>

namespace exanb
{

  // a simple vector of grid blocks
  struct SubDomainNeighborBlock
  {
    GridBlock m_grid_block = { {0,0,0} , {0,0,0} };
    long m_rank = 0;
  };

  using SubDomainNeighborVector = std::vector<SubDomainNeighborBlock>;

  // simple cost model where the cost of a cell is the number of particles in it
  // 
  template<class GridT>
  struct SubDomainNeighbors : public OperatorNode
  {  
    ADD_SLOT( MPI_Comm  , mpi         , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GridT     , grid        , INPUT , REQUIRED );
    ADD_SLOT( GridBlock , lb_block    , INPUT , OPTIONAL );
    ADD_SLOT( SubDomainNeighborVector , subdomain_neighbors , OUTPUT , SubDomainNeighborVector{} );

    inline bool is_sink() const override final { return true; }

    inline void execute () override final
    {
      if( ! lb_block.has_value() ) return;
    
      MPI_Comm comm = *mpi;
      int np=1, rank=0;
      MPI_Comm_size(comm,&np);
      MPI_Comm_rank(comm,&rank);
      
      int gl = grid->ghost_layers();

      subdomain_neighbors->resize( np , SubDomainNeighborBlock{} );
      SubDomainNeighborBlock myself = { enlarge_block(*lb_block,gl) , rank };
      MPI_Allgather( (char*) &myself, sizeof(SubDomainNeighborBlock), MPI_CHAR, (char*) subdomain_neighbors->data(), sizeof(SubDomainNeighborBlock), MPI_CHAR, comm);
      assert( subdomain_neighbors->at(rank).m_rank == myself.m_rank );
      
      int j = 0;
      for(int i=0;i<np;i++)
      {
        if( i != rank && ! is_empty( intersection( myself.m_grid_block , subdomain_neighbors->at(i).m_grid_block ) ) )
        {
          subdomain_neighbors->at(i).m_rank = i;
          subdomain_neighbors->at( j++ ) = subdomain_neighbors->at(i);
        }
      }
      
      subdomain_neighbors->resize(j);
      
      ldbg << "Subdomain "<< *lb_block << " neighbors :" << std::endl;
      for(const auto& b : *subdomain_neighbors)
      {
        ldbg << "\t P"<<b.m_rank<<" : nbh block "<< enlarge_block(b.m_grid_block,-gl) << " : ghost " << intersection( myself.m_grid_block , b.m_grid_block ) << std::endl;
      }
      
    }

  };

  // === register factory ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "subdomain_neighbors", make_grid_variant_operator< SubDomainNeighbors > );
  }

}

