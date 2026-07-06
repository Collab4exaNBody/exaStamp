/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <chrono>
#include <ctime>
#include <mpi.h>
#include <string>
#include <numeric>

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
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/check_particles_inside_cell.h>

namespace exaStamp
{
  using namespace exanb;

  // Read XYZ files in extended XYZ format (OVITO/ASE compatible).
  // Format:
  //   Line 1 : number of particles
  //   Line 2 : key=value pairs, notably:
  //              Lattice="a1 a2 a3 b1 b2 b3 c1 c2 c3"   (row-major 3x3 cell matrix)
  //              Properties=species:S:1:pos:R:3[:velo:R:3:...]
  //   Next lines : columns as described by Properties=

  // -------------------------------------------------------------------------
  // Helper: extract the value of a key from an extended XYZ comment line.
  // Handles both quoted values (key="...") and unquoted tokens (key=val).
  // Returns true and sets `value` on success.
  // -------------------------------------------------------------------------
  static bool extract_key_value( const std::string& line,
                                  const std::string& key,
                                  std::string&       value )
  {
    // Look for key=
    std::string search = key + "=";
    std::string::size_type pos = line.find(search);
    if( pos == std::string::npos ) return false;

    std::string::size_type start = pos + search.size();
    if( start >= line.size() ) return false;

    if( line[start] == '"' )
    {
      // Quoted value: find closing quote
      std::string::size_type end = line.find('"', start + 1);
      if( end == std::string::npos ) return false;
      value = line.substr(start + 1, end - start - 1);
    }
    else
    {
      // Unquoted: read until whitespace or end
      std::string::size_type end = line.find_first_of(" \t\r\n", start);
      value = line.substr(start, end - start);
    }
    return true;
  }

  // -------------------------------------------------------------------------
  // Helper: parse the Properties= string and determine column indices for
  // position and (optionally) velocity.
  //
  // Properties format: name:type:count:name:type:count:...
  //   type: S (string), R (real), I (integer)
  //
  // Returns false if pos columns are not found.
  // vel_col is set to -1 if velocities are absent.
  // -------------------------------------------------------------------------
  static bool parse_properties( const std::string& props,
                                  int& pos_col,
                                  int& vel_col )
  {
    pos_col = -1;
    vel_col = -1;

    // Split by ':'
    std::vector<std::string> tokens;
    {
      std::stringstream ss(props);
      std::string tok;
      while( std::getline(ss, tok, ':') )
        tokens.push_back(tok);
    }

    // Groups of 3: name, type, count
    int col = 0;
    for( size_t i = 0; i + 2 < tokens.size(); i += 3 )
    {
      const std::string& name  = tokens[i];
      // const std::string& type  = tokens[i+1];  // S/R/I — not needed here
      int count = std::stoi(tokens[i+2]);

      // Species column is always 1 wide and comes first in practice,
      // but we track it anyway so col is correct.
      if( name == "species" || name == "type" )
      {
        // species column: col already correct, just advance
      }
      else if( name == "pos" )
      {
        pos_col = col;
      }
      else if( name == "velo" || name == "vel" || name == "velocities" )
      {
        vel_col = col;
      }

      col += count;
    }

    return ( pos_col >= 0 );
  }

  template<typename GridT>
  class ReadXYZwXFormNode : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi             , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( std::string     , filename        , INPUT , REQUIRED );
    ADD_SLOT( Domain          , domain          , INPUT_OUTPUT );
    ADD_SLOT( GridT           , grid            , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species         , INPUT ); // optional
    ADD_SLOT( ReadBoundsSelectionMode, bounds_mode , INPUT , ReadBoundsSelectionMode::FILE_BOUNDS );
    ADD_SLOT( bool            , read_velocities , INPUT , false );

  public:
    inline void execute () override final
    {
      std::string file_name = onika::data_file_path( *filename );
      Domain& domain = *(this->domain);
      GridT& grid    = *(this->grid);

      // Print basename for log
      std::string basename;
      std::string::size_type p = file_name.rfind("/");
      if( p != std::string::npos ) basename = file_name.substr(p+1);
      else basename = file_name;
      lout << "======== " << basename << " ========" << std::endl;

      using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_id, field::_type>;
      using ParticleTuple   = decltype( grid.cells()[0][0] );

      assert( grid.number_of_particles() == 0 );

      int rank = 0, np = 1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      std::string line;
      uint64_t n_particles = 0;

      std::map<std::string,unsigned int> typeMap;
      unsigned int nextTypeId = 0;
      if( species.has_value() )
      {
        for( size_t i = 0; i < species->size(); i++ )
          typeMap[species->at(i).m_name] = i;
        nextTypeId = species->size();
      }

      std::vector<ParticleTupleIO> particle_data;

      // Velocity storage (parallel array to particle_data, rank-0 only)
      std::vector<Vec3d> velocity_data;

      bool uniform_scale = false;
      Mat3d domain_xform = make_identity_matrix();

      if( rank == 0 )
      {
        std::ifstream file;
        file.open(file_name, std::ifstream::in);
        if( !file.is_open() )
        {
          lerr << "Error in reading xyz: file " << file_name << " not found!" << std::endl;
          std::abort();
        }

        // ----------------------------------------------------------------
        // Line 1: number of atoms
        // ----------------------------------------------------------------
        ssize_t n_atoms = -1;
        std::getline(file, line);
        std::stringstream(line) >> n_atoms;

        // ----------------------------------------------------------------
        // Line 2: extended XYZ comment — parse Lattice= and Properties=
        // ----------------------------------------------------------------
        std::getline(file, line);

        // --- Parse Lattice="..." ---
        std::string lattice_str;
        if( !extract_key_value(line, "Lattice", lattice_str) )
        {
          lerr << "Error: could not find 'Lattice=' key in XYZ comment line." << std::endl;
          lerr << "  Comment line was: " << line << std::endl;
          std::abort();
        }

        Mat3d H = make_identity_matrix();
        {
          std::stringstream ss(lattice_str);
          // Lattice is stored row-major: a1 a2 a3 b1 b2 b3 c1 c2 c3
          ss >> H.m11 >> H.m12 >> H.m13
             >> H.m21 >> H.m22 >> H.m23
             >> H.m31 >> H.m32 >> H.m33;
          if( ss.fail() )
          {
            lerr << "Error: could not parse 9 values from Lattice string: \"" << lattice_str << "\"" << std::endl;
            std::abort();
          }
        }
        lout << "H = " << H << std::endl;

        Vec3d a = Vec3d{H.m11, H.m12, H.m13};
        Vec3d b = Vec3d{H.m21, H.m22, H.m23};
        Vec3d c = Vec3d{H.m31, H.m32, H.m33};

        double box_size_x = norm(a);
        double box_size_y = norm(b);
        double box_size_z = norm(c);

        lout << "a = " << box_size_x << std::endl;
        lout << "b = " << box_size_y << std::endl;
        lout << "c = " << box_size_z << std::endl;

        // --- Parse Properties=... ---
        std::string props_str;
        int pos_col = -1, vel_col = -1;
        if( extract_key_value(line, "Properties", props_str) )
        {
          if( !parse_properties(props_str, pos_col, vel_col) )
          {
            lerr << "Error: could not find 'pos' columns in Properties string: \"" << props_str << "\"" << std::endl;
            std::abort();
          }
        }
        else
        {
          // No Properties key: assume default layout "species:S:1:pos:R:3"
          lout << "Warning: no 'Properties=' key found; assuming layout: species x y z" << std::endl;
          pos_col = 1;
          vel_col = -1;
        }

        // Validate velocity availability
        if( *read_velocities && vel_col < 0 )
        {
          lerr << "Error: 'read_velocities' is true but no velocity column (velo/vel/velocities) "
               << "found in Properties string: \"" << props_str << "\"" << std::endl;
          std::abort();
        }

        if( *read_velocities )
          lout << "Reading velocities from column " << vel_col << std::endl;

        // ----------------------------------------------------------------
        // Read particle data
        // ----------------------------------------------------------------
        Mat3d Ht = transpose(H);

        while( std::getline(file, line) )
        {
          if( line.empty() ) continue;

          // Tokenise the line
          std::vector<std::string> toks;
          {
            std::stringstream ss(line);
            std::string tok;
            while( ss >> tok ) toks.push_back(tok);
          }
          if( toks.empty() ) continue;

          // Column 0 is always the species/type label
          std::string type = toks[0];

          // Position columns
          if( pos_col + 2 >= (int)toks.size() )
          {
            lerr << "Error: not enough columns for position on line: " << line << std::endl;
            std::abort();
          }
          double x = std::stod(toks[pos_col    ]);
          double y = std::stod(toks[pos_col + 1]);
          double z = std::stod(toks[pos_col + 2]);

          // Fold into fractional coords and back to Cartesian
          Vec3d r{x, y, z};
          r = inverse(Ht) * r;

          if( r.x < 0. ) r.x += 1.;
          if( r.y < 0. ) r.y += 1.;
          if( r.z < 0. ) r.z += 1.;
          if( r.x >= 1. ) r.x -= 1.;
          if( r.y >= 1. ) r.y -= 1.;
          if( r.z >= 1. ) r.z -= 1.;

          x = box_size_x * r.x;
          y = box_size_y * r.y;
          z = box_size_z * r.z;

          // Type mapping
          if( typeMap.find(type) == typeMap.end() )
          {
            typeMap[type] = nextTypeId;
            ++nextTypeId;
          }

          particle_data.push_back( ParticleTupleIO(x, y, z, n_particles++, typeMap[type]) );

          // Velocity columns (optional)
          if( *read_velocities )
          {
            if( vel_col + 2 >= (int)toks.size() )
            {
              lerr << "Error: not enough columns for velocity on line: " << line << std::endl;
              std::abort();
            }
            double vx = std::stod(toks[vel_col    ]);
            double vy = std::stod(toks[vel_col + 1]);
            double vz = std::stod(toks[vel_col + 2]);
            velocity_data.push_back( Vec3d{vx, vy, vz} );
          }
        }

        // ----------------------------------------------------------------
        // Domain bounds / xform
        // ----------------------------------------------------------------
        if( !domain.xform_is_identity() )
        {
          lerr << "needs initial XForm, resetting XForm" << std::endl;
          domain.set_xform( make_identity_matrix() );
        }

        AABB file_bounds = { {0., 0., 0.}, {box_size_x, box_size_y, box_size_z} };
        compute_domain_bounds(domain, *bounds_mode, 0.0, file_bounds, file_bounds, true);

        if( !domain.xform_is_identity() )
        {
          Mat3d inv_xform = domain.inv_xform();
          for( auto& pp : particle_data )
          {
            Vec3d rr = inv_xform * Vec3d{ pp[field::rx], pp[field::ry], pp[field::rz] };
            pp[field::rx] = rr.x;
            pp[field::ry] = rr.y;
            pp[field::rz] = rr.z;
          }
        }

        const Mat3d D    = diag_matrix(Vec3d{box_size_x, box_size_y, box_size_z});
        const Mat3d Hbis = Ht * inverse(D);
        domain_xform = Hbis * domain.xform();

        uniform_scale = is_uniform_scale(domain_xform);
        lout << "Uniform scale    = " << std::boolalpha << uniform_scale << std::endl;
        if( uniform_scale )
        {
          domain.set_xform( make_identity_matrix() );
          domain.set_cell_size( domain.cell_size() * domain_xform.m11 );
          domain.set_bounds( { domain.origin() * domain_xform.m11, domain.extent() * domain_xform.m11 } );
        }
        else
        {
          domain.set_xform(domain_xform);
        }

        lout << "Particles        = " << particle_data.size()          << std::endl;
        lout << "Domain XForm     = " << domain.xform()                << std::endl;
        lout << "Domain bounds    = " << domain.bounds()               << std::endl;
        lout << "Domain size      = " << bounds_size(domain.bounds())  << std::endl;
        lout << "Real size        = " << bounds_size(domain.bounds())
                                        * Vec3d{domain.xform().m11,
                                                domain.xform().m22,
                                                domain.xform().m33}    << std::endl;
        lout << "Cell size        = " << domain.cell_size()            << std::endl;
        lout << "Grid dimensions  = " << domain.grid_dimension()
             << " (" << grid_cell_count(domain.grid_dimension()) << " cells)" << std::endl;
      } // rank == 0

      // ----------------------------------------------------------------
      // Broadcast domain to all ranks
      // ----------------------------------------------------------------
      MPI_Bcast(&domain, sizeof(Domain), MPI_CHARACTER, 0, *mpi);
      assert( check_domain(domain) );

      grid.set_offset( IJK{0,0,0} );
      grid.set_origin( domain.bounds().bmin );
      grid.set_cell_size( domain.cell_size() );
      grid.set_dimension( domain.grid_dimension() );

      if( rank == 0 )
      {
        for( size_t i = 0; i < particle_data.size(); i++ )
        {
          auto& pp = particle_data[i];
          Vec3d r = { pp[field::rx], pp[field::ry], pp[field::rz] };
          if( uniform_scale ) r = r * domain_xform.m11;

          IJK loc = domain_periodic_location(domain, r);
          assert( grid.contains(loc) );
          assert( min_distance2_between(r, grid.cell_bounds(loc)) < grid.epsilon_cell_size2() );

          pp[field::rx] = r.x;
          pp[field::ry] = r.y;
          pp[field::rz] = r.z;

          ParticleTuple t = pp;
          grid.cell(loc).push_back(t);

          // Store velocity into the particle if requested
          // (assumes field::vx/vy/vz exist in the grid's field set)
          if( *read_velocities )
          {
            const Vec3d& v = velocity_data[i];
            auto& cell = grid.cell(loc);
            size_t idx  = cell.size() - 1;
            cell[idx][field::vx] = v.x;
            cell[idx][field::vy] = v.y;
            cell[idx][field::vz] = v.z;
          }
        }
      }

      lout << "============================" << std::endl;

      grid.rebuild_particle_offsets();

#     ifndef NDEBUG
      bool particles_inside_cell = check_particles_inside_cell(grid);
      assert( particles_inside_cell );
#     endif
    }

  };

  // === register factories ===
  __attribute__((constructor)) static void register_factories()
  {
    OperatorNodeFactory::instance()->register_factory(
      "read_xyz_file_with_xform",
      make_grid_variant_operator< ReadXYZwXFormNode >
    );
  }

} // namespace exaStamp
