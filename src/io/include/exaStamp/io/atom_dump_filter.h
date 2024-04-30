#pragma once

#include <exanb/fields.h>
#include <exanb/core/field_set_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/soatl/field_tuple.h>
#include <exanb/io/sim_dump_io.h>
#include <iostream>

#include <exanb/core/basic_types_operators.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/quaternion_operators.h>
#include <exanb/core/quaternion_stream.h>
#include <exanb/core/string_utils.h>

namespace exaStamp
{
  using namespace exanb;

  struct NullOptionalHeaderIO
  {
    template<class WriteFuncT> inline size_t write_optional_header( WriteFuncT ) { return 0; }
    template<class ReadFuncT> inline size_t read_optional_header( ReadFuncT ) { return 0; }
  };

  template< class GridT , class DumpFieldSet , class LDBG, class OptionalHeaderIO=NullOptionalHeaderIO >
  struct AtomDumpFilter
  {
    using GridFieldSet = typename GridT::Fields;
    using TupleT = onika::soatl::FieldTupleFromFieldIds< DumpFieldSet > ; //FieldTuple<field_ids...>;
    static constexpr bool has_field_id = TupleT::has_field(field::id);
    static constexpr bool has_field_type = TupleT::has_field(field::type);
    static constexpr bool has_field_id_type = has_field_id && has_field_type;
    static constexpr size_t MAX_PARTICLE_SPECIES = 1ul << 16;
    using StorageType = TupleT; /*std::conditional_t<
      has_field_id_type , 
      onika::soatl::FieldTupleFromFieldIds< AddDefaultFields< MergeFieldSet< RemoveFields< DumpFieldSet , FieldSet<field::_id,field::_type> > , FieldSet<field::_id_type> > > > , 
      TupleT
      >;  */

    ParticleSpecies& particle_species;
    LDBG& ldbg;

    // read/writes optional header data (i.e. useful for molecule species)
    OptionalHeaderIO optional_header_io = {};

    // optional cell size scaling
    double scale_cell_size = 1.0;

    // optional modification of domain periodicity
    bool override_periodicity = false;
    bool periodic_x = false;
    bool periodic_y = false;
    bool periodic_z = false;

    // optional override of domain expandability
    bool override_expandable = false;
    bool expandable = false;
    
    // optionally keep additional species from user
    bool keep_species = true;
    
    // optionally override domain bounds
    bool override_domain_bounds = false;
    bool shrink_to_fit = false;
    AABB domain_bounds = { {0,0,0} , {0,0,0} };

    // optionally override domain mirror boundary flags
    bool override_mirroring = false;
    bool mirror_x_min = false;
    bool mirror_x_max = false;
    bool mirror_y_min = false;
    bool mirror_y_max = false;
    bool mirror_z_min = false;
    bool mirror_z_max = false;

    // stats
    StorageType m_tuple_min;
    StorageType m_tuple_max;
    StorageType m_tuple_avg;
    size_t m_n_tuples = 0;

    struct StatsUpdater
    {
      AtomDumpFilter& m_ref;
      template<class F , class T> inline void operator () ( F f, T value )
      {
        static constexpr bool has_sum_div = std::is_arithmetic_v< typename F::value_type > || std::is_same_v< typename F::value_type , Vec3d > || std::is_same_v< typename F::value_type , Mat3d > || std::is_same_v< typename F::value_type , Quaternion >;
        static constexpr bool has_min_max = std::is_arithmetic_v< typename F::value_type > || std::is_same_v< typename F::value_type , Vec3d > ;
        if constexpr ( has_sum_div )
        {
          m_ref.m_tuple_avg[f] = m_ref.m_tuple_avg[f] + value;
        }
        if constexpr ( has_min_max )
        {
          if( m_ref.m_n_tuples == 0 )
          {
            m_ref.m_tuple_max[f] = value;
            m_ref.m_tuple_min[f] = value;
          }
          else
          {
            m_ref.m_tuple_max[f] = std::max( value , m_ref.m_tuple_max[f] );
            m_ref.m_tuple_min[f] = std::min( value , m_ref.m_tuple_min[f] );
          }
        }
        ++ m_ref.m_n_tuples;
      }
    };

    template<class StreamT>
    struct StatsPrinter
    {
      AtomDumpFilter& m_ref;
      StreamT& out;
      template<class F , class T> inline void operator () ( F f, T value )
      {
        using FiledPrintType = std::conditional_t< std::is_integral_v<typename F::value_type> , int64_t , typename F::value_type >;
        static constexpr bool has_sum_div = std::is_arithmetic_v< typename F::value_type > || std::is_same_v< typename F::value_type , Vec3d > || std::is_same_v< typename F::value_type , Mat3d > || std::is_same_v< typename F::value_type , Quaternion >;
        if constexpr ( has_sum_div )
        {
          if( m_ref.m_n_tuples > 0 ) m_ref.m_tuple_avg[f] = m_ref.m_tuple_avg[f] / double(m_ref.m_n_tuples);
        }
        out << f.short_name() << ": " << std::scientific << std::setprecision(6) 
            << static_cast<FiledPrintType>( m_ref.m_tuple_avg[f] ) << " / "
            << static_cast<FiledPrintType>( m_ref.m_tuple_min[f] ) << " / "
            << static_cast<FiledPrintType>( m_ref.m_tuple_max[f] ) << std::endl;
       }
    };
    template<class StreamT> StatsPrinter<StreamT> stats_printer(StreamT & out) { return StatsPrinter<StreamT>{*this,out}; }

    inline void process_domain(Domain& domain, Mat3d& particle_read_xform)
    {
      particle_read_xform = make_identity_matrix();
      if( scale_cell_size != 1.0 )
      {
        const Vec3d dom_size = domain.bounds_size();
        const Mat3d dom_xform = domain.xform();

        double desired_cell_size = domain.cell_size() * scale_cell_size;
        IJK grid_dims = make_ijk( Vec3d{0.5,0.5,0.5} + ( dom_size / desired_cell_size ) ); // round to nearest
        domain.set_grid_dimension( grid_dims );
        domain.set_cell_size( desired_cell_size );
        // domain bounds should remain the same
        domain.set_bounds( { domain.origin() , domain.origin() + ( grid_dims * desired_cell_size ) } );
        const Vec3d dom_new_size = domain.bounds_size();

        particle_read_xform = diag_matrix( dom_new_size / dom_size );
        const Mat3d dom_new_xform = inverse( particle_read_xform ) * dom_xform;
        domain.set_xform( dom_new_xform );
      }
      if( override_periodicity )
      {
        domain.set_periodic_boundary( periodic_x , periodic_y , periodic_z );
      }
      if( override_expandable )
      {
        domain.set_expandable( expandable );
      }
      if( override_domain_bounds )
      {
        if( ! domain.xform_is_identity() )
        {
          if( ! is_diagonal( domain.xform() ) )
          {
            fatal_error() << "cannot force domain bounds on a domain with non diagonal transform matrix" << std::endl;
          }
          else
          {
            Vec3d diag = { domain.xform().m11 , domain.xform().m22 , domain.xform().m33 };
            domain_bounds.bmin = domain_bounds.bmin / diag;
            domain_bounds.bmax = domain_bounds.bmax / diag;
          }
        }
        if( domain.periodic_boundary_x() ) { domain_bounds.bmin.x = domain.bounds().bmin.x; domain_bounds.bmax.x = domain.bounds().bmax.x; }
        if( domain.periodic_boundary_y() ) { domain_bounds.bmin.y = domain.bounds().bmin.y; domain_bounds.bmax.y = domain.bounds().bmax.y; }
        if( domain.periodic_boundary_z() ) { domain_bounds.bmin.z = domain.bounds().bmin.z; domain_bounds.bmax.z = domain.bounds().bmax.z; }
        std::cout << "overriden domain bounds : " << domain_bounds << std::endl;
      }
      if( override_mirroring )
      {
        domain.set_mirror_x_min( mirror_x_min ); domain.set_mirror_x_max( mirror_x_max );
        domain.set_mirror_y_min( mirror_y_min ); domain.set_mirror_y_max( mirror_y_max );
        domain.set_mirror_z_min( mirror_z_min ); domain.set_mirror_z_max( mirror_z_max );
      }
    }

    inline bool particle_input_filter(const Vec3d& r)
    {
      return ( ! override_domain_bounds ) || ( is_inside(domain_bounds,r) );
    }


    inline void update_sats( const StorageType & stp )
    {
      stp.apply_fields( StatsUpdater{*this} );
    }

    inline void initialize_write()
    {
      m_n_tuples = 0;
    }
    inline void initialize_read()
    {
      m_n_tuples = 0;
    }
    
    inline void post_process_domain(Domain& domain)
    {
      if( override_domain_bounds && shrink_to_fit)
      {
        ldbg << "shrinking domain bounds and grid to fit user fixed bounds" << std::endl;
        const IJK cell_shift = make_ijk( floor( ( domain_bounds.bmin - domain.origin() ) / domain.cell_size() ) );
        ldbg << "domain grid shift = "<<cell_shift<<std::endl;
        const Vec3d new_origin = domain.origin() + cell_shift * domain.cell_size();
        ldbg << "new_origin = "<<new_origin<<std::endl;
        const IJK new_grid_dims = make_ijk( ceil( ( domain_bounds.bmax - new_origin ) / domain.cell_size() ) );
        ldbg << "new_grid_dims = "<<new_grid_dims<<std::endl;
        const Vec3d new_extent = new_origin + new_grid_dims * domain.cell_size();
        ldbg << "new_extent = "<<new_extent<<std::endl;
        domain.set_bounds( { new_origin , new_extent } );
        domain.set_grid_dimension( new_grid_dims );
      }
    }
    
    inline void finalize_read()
    {    
      ldbg<< "field : avg / min / max" << std::endl;
      m_tuple_avg.apply_fields( stats_printer(ldbg) );      
    }
    inline void finalize_write()
    {
      ldbg<< "field : avg / min / max" << std::endl;
      m_tuple_avg.apply_fields( stats_printer(ldbg) );      
    }

    template<class WriteFuncT>
    inline size_t write_optional_header( WriteFuncT write_func )
    {
      size_t n = 0;
      if( particle_species.size() >= MAX_PARTICLE_SPECIES )
      {
        fatal_error() << "Number of particle species too big ("<<particle_species.size()<<")" <<std::endl;
      }
      size_t n_species = particle_species.size() + MAX_PARTICLE_SPECIES * ParticleSpecie::MaxRigidMolAtoms;
      ldbg << "write modified n_species : "<<n_species<<" = "<<particle_species.size()<<" + "<<MAX_PARTICLE_SPECIES<<" * "<< ParticleSpecie::MaxRigidMolAtoms<<std::endl;
      n += write_func( n_species );
      for(const auto& sp : particle_species) { n += write_func( sp ); }
      n += optional_header_io.write_optional_header( write_func );
      return n;
    }

    template<class ReadFuncT>
    inline size_t read_optional_header( ReadFuncT read_func )
    {
      ParticleSpecies& species = particle_species;
      ParticleSpecies backup_species = species;
      size_t n = 0;
      size_t n_species = 0;
      size_t max_rigid_molecule = LEGACY_MAX_RIGID_MOLECULE_ATOMS;
      n += read_func( n_species );
      if( n_species >= MAX_PARTICLE_SPECIES )
      {
        max_rigid_molecule = n_species / MAX_PARTICLE_SPECIES;
        n_species = n_species % MAX_PARTICLE_SPECIES;
        ldbg << "decode n_species in file to n_species="<<n_species<<" , max_rigid_molecule="<<max_rigid_molecule<<std::endl;
      }
      species.resize( n_species );
      
      if( max_rigid_molecule == ParticleSpecie::MaxRigidMolAtoms )
      {
        for(auto& sp : species) { n += read_func( sp ); }
      }
      else if( max_rigid_molecule == LEGACY_MAX_RIGID_MOLECULE_ATOMS )
      {
        for(auto& sp : species)
        {
          ParticleSpecieTmpl<LEGACY_MAX_RIGID_MOLECULE_ATOMS> tmp;
          n += read_func( tmp );
          if( tmp.m_rigid_atom_count > ParticleSpecie::MaxRigidMolAtoms )
          {
            fatal_error() << "Rigid molecule "<<tmp.name()<<" has more atoms than supported : "<<tmp.m_rigid_atom_count<<" > "<<ParticleSpecie::MaxRigidMolAtoms<<std::endl;
          }
          sp.copy_from( tmp );
        }
      }
      else
      {
        fatal_error()<<"Rigid molecule max atoms in dump ("<<max_rigid_molecule<<") is not supported in this software version (compiled for "<<ParticleSpecie::MaxRigidMolAtoms<<")"<<std::endl;
      }
      if( keep_species )
      {
        auto find_species = [&species](const std::string& name) -> bool { for(const auto& sp:species) { if(sp.name()==name) return true; } return false; };
        for(const auto& sp:backup_species)
        {
          if( ! find_species( sp.name() ) )
          {
            species.push_back( sp );
          }
        }
      }
      n += optional_header_io.read_optional_header( read_func );
      return n;
    }

    inline StorageType encode( const TupleT & tp )
    {
      /*
      StorageType stp = tp;
      if constexpr ( has_field_id && has_field_type )
      {
        stp[ field::id_type ] = ( tp[ field::id ] << 8 ) & tp[ field::type ];
      }
      return stp;
      */
      update_sats( tp );
      return tp;
    }

    inline TupleT decode( const StorageType & stp )
    {
      update_sats( stp );
      /*
      TupleT tp = stp;
      if constexpr ( has_field_id_type )
      {
        tp[ field::id ] = stp[ field::id_type ] >> 8 ;
        tp[ field::type ] = stp[ field::id_type ] ;
      }
      return tp;
      */
      return stp;
    }

    static inline constexpr size_t optional_cell_data_size(size_t) { return 0; }
    static inline constexpr void write_optional_cell_data(uint8_t*, size_t) { }
    static inline constexpr void append_cell_particle( size_t, size_t ) {};
    static inline constexpr void read_optional_data_from_stream( const uint8_t* , size_t ) {}
    static inline constexpr const uint8_t* optional_cell_data_ptr(size_t) { return nullptr; }
  };

  template< class GridT , class DumpFieldSet , class LDBG >
  inline AtomDumpFilter<GridT,DumpFieldSet,LDBG> make_atom_dump_filter( GridT& grid , ParticleSpecies& species, LDBG& ldbg , DumpFieldSet )
  {
    return { species, ldbg };
  }

}

