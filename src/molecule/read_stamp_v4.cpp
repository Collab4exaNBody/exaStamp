#include <onika/math/basic_types.h>
#include <onika/math/basic_types_stream.h>
#include <onika/math/basic_types_yaml.h>
#include <exanb/core/file_utils.h>
#include <exanb/core/grid.h>
#include <onika/log.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/molecule/molecule_species.h>

#include <exanb/core/quantity.h>

#include <exaStamp/molecule/read_stamp_v4.h>

#include <mpi.h>
#include <fstream>

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT>
  struct ReadStampV4 : public OperatorNode
  {
    ADD_SLOT(MPI_Comm               , mpi           , INPUT , MPI_COMM_WORLD                       , DocString{"MPI communicator"} );
    ADD_SLOT(std::string            , file          , INPUT                                        , DocString{"File name"} );
    ADD_SLOT(ReadBoundsSelectionMode, bounds_mode   , INPUT , ReadBoundsSelectionMode::FILE_BOUNDS , DocString{"Bounds mode"} );
    ADD_SLOT(double                 , enlarge_bounds, INPUT , 0.0                                  , DocString{"Bounds enlargement"} );
    ADD_SLOT(Domain                 , domain        , INPUT_OUTPUT                                 , DocString{"Domain"} );
    ADD_SLOT(GridT                  , grid          , INPUT_OUTPUT                                 , DocString{"Particle grid"} );
    ADD_SLOT(long                   , timestep      , OUTPUT                                       , DocString{"Iteration number"} );
    ADD_SLOT(double                 , dt            , OUTPUT                                       , DocString{"Time step"} );
    ADD_SLOT(double                 , physical_time , OUTPUT                                       , DocString{"Physical time"} );
    ADD_SLOT(ParticleSpecies        , species       , INPUT                                        , DocString{"Particle species data block"} );
    ADD_SLOT(bool                   , enable_xform_scale, INPUT, false                             , DocString{"enable scaling of reduced domain to obtain identity xform when possible"} );
    ADD_SLOT(bool                   , pbc_adjust_xform , INPUT , false                             , DocString{"Adjust xform to preserve sizes in directions with periodic boundary conditions"} );
    ADD_SLOT(long                   , version       , INPUT , -1                                   , DocString{"force file version (ignore version number in file) ex: 41 for 4.1, 42 for 4.2"} );
    ADD_SLOT(bool                   , molrig42      , INPUT , true                                 , DocString{"force short (without orientation) 4.2 rigid molecule structures"} );

    // -----------------------------------------------
    // ----------- Operator documentation ------------
    inline std::string documentation() const override final
    {
      return R"EOF(
        Reads a protection file in StampV4 format: 
        this format is fully compatible with the Stamp4 MD code (read/write)
        )EOF";
    }

    inline void execute () override final
    {
      std::string file_name = data_file_path( *file );
      read_stamp_v4( ldbg, *mpi,file_name,*enlarge_bounds,*bounds_mode,*grid,*domain
                   , *timestep,*dt,*physical_time,*species
                   , *enable_xform_scale , *pbc_adjust_xform , *version , *molrig42 );
    }

    inline void yaml_initialize(const YAML::Node& node) override final
    {
      YAML::Node tmp;
      if( node.IsScalar() )
      {
        tmp["file"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize( tmp );
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(read_stamp_v4)
  {
    OperatorNodeFactory::instance()->register_factory( "read_stamp_v4", make_grid_variant_operator< ReadStampV4 > );
  }

}
