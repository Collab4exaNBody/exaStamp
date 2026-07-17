/*
   Licensed to the Apache Software Foundation (ASF) under one
   or more contributor license agreements.  See the NOTICE file
   distributed with this work for additional information
   regarding copyright ownership.  The ASF licenses this file
   to you under the Apache License, Version 2.0 (the
   "License"); you may not use this file except in compliance
   with the License.  You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <mpi.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
#include <onika/file_utils.h>
#include <exanb/core/domain.h>
#include <exanb/core/check_particles_inside_cell.h>
#include <exanb/core/simple_block_rcb.h>
#include <exanb/core/particle_type_id.h>

#include <ctime>
#include <string>
#include <numeric>

// rsa_mpi stuff
#include <rsa_data_storage.hxx>
#include <rsa_random.hxx>
#include <rsa_domain.hxx>
#include <rsa_decoration.hxx>
#include <operator_algorithm.hxx>
#include <radius_generator.hxx>

namespace exaStamp
{
  // accepts either a raw integer type id (type: 2) or a species name (type: Ta),
  // resolved against particle_type_map at execute() time.
  struct RSAParticleType
  {
    int id = 0;
    std::string name;
    bool is_id = true;
  };

  // one entry of a multi-species RSA packing: [radius, nb_particles, type]
  struct RSAParticleSpec
  {
    double radius = 0.0;
    long nb_particles = 0;
    RSAParticleType type;
  };
}

namespace YAML
{
  template<> struct convert< exaStamp::RSAParticleType >
  {
    static bool decode(const Node& node, exaStamp::RSAParticleType& v)
    {
      if( ! node.IsScalar() ) return false;
      try
      {
        v.id = node.as<int>();
        v.is_id = true;
      }
      catch( ... )
      {
        v.name = node.as<std::string>();
        v.is_id = false;
      }
      return true;
    }
  };

  template<> struct convert< exaStamp::RSAParticleSpec >
  {
    static bool decode(const Node& node, exaStamp::RSAParticleSpec& v)
    {
      if( ! node.IsSequence() || node.size() != 3 ) return false;
      v.radius = node[0].as<double>();
      v.nb_particles = node[1].as<long>();
      v.type = node[2].as<exaStamp::RSAParticleType>();
      return true;
    }
  };
}

namespace exaStamp
{

  using namespace exanb;

  template<
    class GridT
    >

  class InitRSA : public OperatorNode
  {
    ADD_SLOT( MPI_Comm           , mpi, INPUT, MPI_COMM_WORLD);
    ADD_SLOT( GridT              , grid, INPUT_OUTPUT);
    ADD_SLOT( Domain             , domain, INPUT_OUTPUT);
    ADD_SLOT( double             , enlarge_bounds, INPUT, 0.0);
    ADD_SLOT( std::vector<bool>  , periodicity, INPUT, OPTIONAL, DocString{"if set, overrides domain's periodicity stored in file with this value"});
    ADD_SLOT( bool               , expandable, INPUT, OPTIONAL, DocString{"if set, override domain expandability stored in file"});
    ADD_SLOT( AABB               , bounds, INPUT, REQUIRED, DocString{"if set, override domain's bounds, filtering out particle outside of overriden bounds"});
    ADD_SLOT( RSAParticleType    , type, INPUT, RSAParticleType{}, DocString{"single-species particle type, either an integer id (type: 2) or a species name (type: Ta), resolved via particle_type_map; ignored when 'species' is set"});
    ADD_SLOT( ParticleTypeMap    , particle_type_map, INPUT, OPTIONAL, DocString{"required only when a type is given as a species name"});
    ADD_SLOT( bool               , pbc_adjust_xform, INPUT, true);
    ADD_SLOT( double             , radius, INPUT, OPTIONAL, DocString{"single-species radius; ignored when 'species' is set"});
    ADD_SLOT( long                , nb_particles, INPUT, OPTIONAL, DocString{"single-species exact number of particles to place (RSA may place fewer if the domain can't fit that many); if unset (and 'species' isn't set either), packs the domain as densely as possible; ignored when 'species' is set"});
    ADD_SLOT( std::vector<RSAParticleSpec>, species, INPUT, OPTIONAL, DocString{"multi-species RSA packing: list of [radius, nb_particles, type] entries, e.g. species: [[0.5, 100, Ta], [0.25, 200, Cu]]; overrides radius/type/nb_particles when set"});
    ADD_SLOT( long                , seed, INPUT, 0, DocString{"seed for the RSA pseudo-random generator"});
    ADD_SLOT( bool                , verbose, INPUT, false, DocString{"enable rsa_mpi per-round diagnostics (miss rate, shot counts)"});
    ADD_SLOT( double             , rcut_max, INPUT_OUTPUT, 0.0);
    
 public:
    
  inline void execute() override final
    {
      //-------------------------------------------------------------------------------------------
      using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_id, field::_type>;
      using ParticleTuple = decltype(grid->cells()[0][0]);
      assert(grid->number_of_particles() == 0);

      // MPI Initialization
      int rank = 0, np = 1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      if( ! species.has_value() && ! radius.has_value() )
      {
        fatal_error() << "init_rsa: either 'radius' or 'species' must be set" << std::endl;
      }

      auto resolve_type = [&]( const RSAParticleType& t ) -> int
      {
        if( t.is_id ) return t.id;
        if( ! particle_type_map.has_value() )
        {
          fatal_error() << "init_rsa: type '"<<t.name<<"' given as a species name but no particle_type_map is available (run a 'species' operator upstream)" << std::endl;
        }
        auto it = particle_type_map->find(t.name);
        if( it == particle_type_map->end() )
        {
          fatal_error() << "init_rsa: unknown species '"<<t.name<<"'" << std::endl;
        }
        return it->second;
      };

      AABB b = *bounds;
      constexpr int DIM = 3;
      constexpr int method = 1;
      constexpr int ghost_layer = 1;
      std::array<double, DIM> domain_inf = {b.bmin.x, b.bmin.y, b.bmin.z};
      std::array<double, DIM> domain_sup = {b.bmax.x, b.bmax.y, b.bmax.z};

      double r_max = species.has_value() ? 0.0 : *radius;
      if( species.has_value() )
      {
        for( const auto& sp : *species ) r_max = std::max(r_max, sp.radius);
      }
      rsa_domain<DIM> rsa_domain(domain_inf, domain_sup, ghost_layer, r_max);

      int ParticleType = species.has_value() ? 0 : resolve_type(*type);

      size_t rsa_seed = size_t(*seed);
      if( species.has_value() )
      {
        // multi-species: each entry places its own exact (upper-bound) number of spheres
        std::vector<std::tuple<double,uint64_t,int>> cast_list;
        for( const auto& sp : *species )
        {
          cast_list.push_back( std::make_tuple(sp.radius, uint64_t(sp.nb_particles), resolve_type(sp.type)) );
        }
        sac_de_billes::RadiusGenerator<DIM> radius_generator(cast_list, 0.0);
        if( *verbose )
          algorithm::uniform_generate<DIM, method, true>(rsa_domain, radius_generator, 6000, 10, rsa_seed);
        else
          algorithm::uniform_generate<DIM, method, false>(rsa_domain, radius_generator, 6000, 10, rsa_seed);
      }
      else if( nb_particles.has_value() )
      {
        // pack an exact (upper-bound) number of spheres instead of filling the domain to volume fraction 1
        sac_de_billes::RadiusGenerator<DIM> radius_generator(
          std::vector<std::tuple<double,uint64_t,int>>{ { *radius, uint64_t(*nb_particles), ParticleType } }, 0.0 );
        if( *verbose )
          algorithm::uniform_generate<DIM, method, true>(rsa_domain, radius_generator, 6000, 10, rsa_seed);
        else
          algorithm::uniform_generate<DIM, method, false>(rsa_domain, radius_generator, 6000, 10, rsa_seed);
      }
      else
      {
        if( *verbose )
          algorithm::uniform_generate<DIM, method, true>(rsa_domain, *radius, 6000, 10, rsa_seed);
        else
          algorithm::uniform_generate<DIM, method, false>(rsa_domain, *radius, 6000, 10, rsa_seed);
      }
      auto spheres = rsa_domain.extract_spheres();

      if (rank == 0) {
        /** FILE_BOUNDS sounds wrong in this context, but it works. */
        compute_domain_bounds(*domain, exanb::ReadBoundsSelectionMode::FILE_BOUNDS, *enlarge_bounds, b, b, *pbc_adjust_xform);
      }

      // compute indexes
      int ns = spheres.size();
      MPI_Exscan(MPI_IN_PLACE, &ns, 1, MPI_INT, MPI_SUM, *mpi);
      
      // send bounds and size_box values to all cores
      MPI_Bcast(&(*domain), sizeof(Domain), MPI_CHARACTER, 0, *mpi);
      assert(check_domain(*domain));
      grid->set_offset(IJK{0, 0, 0});
      grid->set_origin(domain->bounds().bmin);
      grid->set_cell_size(domain->cell_size());
      grid->set_dimension(domain->grid_dimension());

      // add particles
      std::vector<ParticleTupleIO> particle_data;
      ParticleTupleIO pt;
      particle_data.resize(spheres.size());
      for (size_t s = 0; s < spheres.size(); s++) {
        auto pos = spheres[s].center;
        auto id = ns + s;
        // in multi-species mode, each sphere's phase carries its resolved particle type
        int pt_type = species.has_value() ? int(spheres[s].phase) : ParticleType;
        pt = ParticleTupleIO(pos[0], pos[1], pos[2], id, pt_type);
        particle_data[s] = pt;
      }

      // Fill grid, particles will migrate accross mpi processed
      // using the operator migrate_cell_particles
      for (auto p : particle_data) {
        Vec3d r{p[field::rx], p[field::ry], p[field::rz]};
        IJK loc = domain_periodic_location(*domain, r);  // grid.locate_cell(r);
        assert(grid->contains(loc));
        assert(min_distance2_between(r, grid->cell_bounds(loc)) < grid->epsilon_cell_size2());
        p[field::rx] = r.x;
        p[field::ry] = r.y;
        p[field::rz] = r.z;
        ParticleTuple t = p;
        //        t[field::homothety] = 1.0;
        grid->cell(loc).push_back(t, grid->cell_allocator());
      }
      
      uint64_t n_particles = particle_data.size();
      uint64_t n;
      MPI_Reduce(&n_particles, &n, 1, MPI_UINT64_T, MPI_SUM, 0, *mpi);
      
      // Display information
      lout << "=================================" << std::endl;
      lout << "Particles        = " << n << std::endl;
      lout << "Domain XForm     = " << domain->xform() << std::endl;
      lout << "Domain bounds    = " << domain->bounds() << std::endl;
      lout << "Domain size      = " << bounds_size(domain->bounds()) << std::endl;
      lout << "Real size        = "
           << bounds_size(domain->bounds()) * Vec3d{domain->xform().m11, domain->xform().m22, domain->xform().m33}
        << std::endl;
      lout << "Cell size        = " << domain->cell_size() << std::endl;
      lout << "Grid dimensions  = " << domain->grid_dimension() << " (" << grid_cell_count(domain->grid_dimension())
           << " cells)" << std::endl;
      lout << "=================================" << std::endl;
      
      grid->rebuild_particle_offsets();
      *rcut_max = std::max(*rcut_max, r_max);
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(init_rsa)
  {
    OperatorNodeFactory::instance()->register_factory("init_rsa", make_grid_variant_operator<InitRSA>);
  }
  
}
