

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/fields.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/log.h>
#include <exanb/core/unityConverterHelper.h>

#include <exanb/io/mpi_file_io.h>
#include <exaStamp/io/stampv4_io.h>
#include <exanb/mpi/all_value_equal.h>
#include <onika/oarray.h>

#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>

//#include "exanb/vector_utils.h"
#include <exanb/core/file_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/molecule/molecule_species.h>

#include <exanb/core/check_particles_inside_cell.h>
#include <exanb/core/parallel_random.h>


#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
#include <limits>
#include <memory>
#include <map>

#include <ctime>
#include <mpi.h>
#include <string>
#include <numeric>

namespace exaStamp
{
  using namespace exanb;

  // Read XYZ files.
  // This files must be ASCII files
  // first line : number of particles
  // second line : xform of simulation domain
  // next lines : typeFF X Y Z molid moltype connectivity[5]
  // example : Cno2 5.3 9.23 -3.5 1 0 6 2 7 -1 -1
  // Note : the typeFF have to match a specie name

  template<typename GridT>
  class ReadXYZwXFormMoleculesNode : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
    ADD_SLOT( std::string     , file         , INPUT , REQUIRED );
    ADD_SLOT( Domain          , domain       , INPUT_OUTPUT );
    ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species      , INPUT ); // optional. if no species given, type ids are allocated automatically
    ADD_SLOT(MoleculeSpeciesVector  , molecules     , INPUT_OUTPUT                                 , DocString{"Molecule descriptions"} );
    ADD_SLOT(bool                   , build_molecule_species , INPUT , false                       , DocString{"detect molecule species from connectivity and fill molecules slot"} );    
    ADD_SLOT( ReadBoundsSelectionMode, bounds_mode   , INPUT , ReadBoundsSelectionMode::FILE_BOUNDS );

  public:
    inline void execute () override final
    {
      // static const double coord_conv = UnityConverterHelper::convert(1.0, "ang"); // conversion factor between input unity to exaStamp internal unity

      //-------------------------------------------------------------------------------------------
      // Reading datas from YAML or previous input
      std::string file_name = data_file_path( *file );
      Domain& domain = *(this->domain);
      GridT& grid = *(this->grid);
      MoleculeSpeciesVector& molvec = *(this->molecules);
      ParticleSpecies& species = *(this->species);
      auto & molspecies = molvec.m_molecules;      

      //-------------------------------------------------------------------------------------------
      std::string basename;
      std::string::size_type p = file_name.rfind("/");
      if( p != std::string::npos ) basename = file_name.substr(p+1);
      else basename=file_name;      
      lout << "======== " << basename << " ========" << std::endl;
      //-------------------------------------------------------------------------------------------

      using IOMolecules = stampv4::IOMolecules;
      using IndexT = std::conditional_t< true , int32_t , int64_t >;      
      using IOAtomes = stampv4::IOAtomes<IndexT>;

      using MoleculeTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, 
						       field::_vx, field::_vy, field::_vz, 
						       field::_id, field::_type, 
						       field::_orient , field::_angmom,
						       field::_idmol, field::_cmol>;
      using MoleculeTuple = decltype( grid.cells()[0][0] );      
      assert( grid.number_of_particles() == 0 );

      // MPI Initialization
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      //useful later
      std::string line;

      uint64_t count = 0;      
      // Get max and min positions
      // Need to define the size of the box
      // NOTE : only one processor need to do that
      std::vector<std::string> type_vec;
      std::vector<std::string> typeFF_vec;
      std::vector<int> molid_vec;
      std::vector<int> moltype_vec;
      std::vector<Vec3d> pos_vec;
      std::vector<Vec3d> vel_vec;
      std::vector<std::array<uint64_t,4>> cmol_vec;
      std::vector<Quaternion> quat_vec;
      std::vector<Vec3d> angmom_vec;
      
      if(rank==0)
      {
        //Open xyz file
        std::ifstream file;
        file.open(file_name, std::ifstream::in);
        if(!file.is_open())
        {
          lerr << "Error in reading xyz : file "<< file_name << " not found !" << std::endl;
          std::abort();
        }

        // Read number of atoms
        getline(file,line);
        std::stringstream(line) >> count;

        double box_size_x=-1.0;
        double box_size_y=-1.0;
        double box_size_z=-1.0;
        Mat3d H = make_identity_matrix();	
        
        std::getline(file,line);

	unsigned quote;
	quote = line.find("\"");
	quote += 1;
	std::string secondPart = line.substr(quote);
	quote = secondPart.find("\"");
	std::string lineclean = secondPart.substr(0,quote);

        std::stringstream(line) >> H.m11 >> H.m12 >> H.m13 >> H.m21 >> H.m22 >> H.m23 >> H.m31 >> H.m32 >> H.m33;

	lout << "H = " << H << std::endl;

        Vec3d a = Vec3d{H.m11, H.m12, H.m13};
        Vec3d b = Vec3d{H.m21, H.m22, H.m23};
        Vec3d c = Vec3d{H.m31, H.m32, H.m33};

        box_size_x = norm(a);
        box_size_y = norm(b);
        box_size_z = norm(c);

	lout << "a = " << box_size_x << std::endl;        
	lout << "b = " << box_size_y << std::endl;        
	lout << "c = " << box_size_z << std::endl;        

	
	AABB file_bounds  = { { 0., 0., 0. } , {box_size_x,box_size_y,box_size_z} };
        compute_domain_bounds(domain,*bounds_mode,0.0, file_bounds ,file_bounds, true );

	Mat3d inv_xform = make_identity_matrix();
        if( !domain.xform_is_identity() )
	  {
	    inv_xform = domain.inv_xform();
	  }
	
        // if( ! domain.xform_is_identity() )
	//   {
	//     lerr << "needs initial XForm, resetting XForm"<<std::endl;
	//     domain.set_xform( make_identity_matrix() );
	//   }
	
	Mat3d Ht = transpose(H);

	type_vec.resize( count );
	typeFF_vec.resize( count );
	molid_vec.resize( count , -1);
	moltype_vec.resize( count , -1);
	pos_vec.resize( count , Vec3d{0.,0.,0.} );
	vel_vec.resize( count , Vec3d{0.,0.,0.} );
	cmol_vec.resize( count );
	quat_vec.resize( count );
	angmom_vec.resize( count );
	
        // read one line at a time
	int cnt=0;
        while (std::getline(file,line))
        {
          std::string type;
          std::string typeFF;	  
	  size_t moleculeid;
	  size_t moleculetype;
          double x=0.0, y=0.0, z=0.0;
	  size_t c0,c1,c2,c3;

          std::stringstream(line) >> type >> typeFF >> moleculeid >> moleculetype >> x >> y >> z >> c0 >> c1 >> c2 >> c3;

	  Vec3d r{x, y, z};	  
	  r = inverse(Ht) * r;
	  
	  if (r.x < 0.) r.x += 1.;
	  if (r.y < 0.) r.y += 1.;
	  if (r.z < 0.) r.z += 1.;

	  if (r.x >= 1.) r.x -= 1.;
	  if (r.y >= 1.) r.y -= 1.;
	  if (r.z >= 1.) r.z -= 1.;

          x = box_size_x * r.x;
          y = box_size_y * r.y;
          z = box_size_z * r.z;

	  r = inv_xform * Vec3d{ x, y, z };
	  
	  typeFF_vec[cnt]=typeFF;
	  type_vec[cnt]=type;
	  molid_vec[cnt]=moleculeid;
	  moltype_vec[cnt]=moleculetype;
	  pos_vec[cnt]=Vec3d{r.x, r.y, r.z};
	  vel_vec[cnt]=Vec3d{0., 0., 0.};
	  quat_vec[cnt]={0., 0., 0., 0.};
	  angmom_vec[cnt]=Vec3d{0., 0., 0.};
	  cmol_vec[cnt] = {c0, c1, c2, c3};
	  cnt++;
        }
	file.close();

	Mat3d D = diag_matrix(Vec3d{box_size_x, box_size_y, box_size_z});
	Mat3d Hbis = Ht * inverse(D);
        domain.set_xform( Hbis * domain.xform());

	lout << "Particles        = "<<count<<std::endl;
	lout << "Domain XForm     = "<<domain.xform()<<std::endl;
	lout << "Domain bounds    = "<<domain.bounds()<<std::endl;
	lout << "Domain size      = "<<bounds_size(domain.bounds()) <<std::endl;
	lout << "Real size        = "<<bounds_size(domain.bounds()) * Vec3d{domain.xform().m11,domain.xform().m22,domain.xform().m33} <<std::endl;
	lout << "Cell size        = "<<domain.cell_size()<<std::endl;
	lout << "Grid dimensions  = "<<domain.grid_dimension()<<" ("<<grid_cell_count(domain.grid_dimension())<<" cells)"<< std::endl;
	
      }
      
      //send bounds and size_box values to all cores
      MPI_Bcast(&domain, sizeof(Domain), MPI_CHARACTER, 0, *mpi);
      assert( check_domain(domain) );

      grid.set_offset( IJK{0,0,0} );
      grid.set_origin( domain.bounds().bmin );
      grid.set_cell_size( domain.cell_size() );
      grid.set_dimension( domain.grid_dimension() );

      if (rank == 0) {

	std::unique_ptr<IOMolecules[]> molecules = nullptr;	
	molecules = std::make_unique<IOMolecules[]>(count);
	std::unique_ptr<IOAtomes[]> atomes = nullptr;
	atomes = std::make_unique<IOAtomes[]>(count);
	size_t cnt=0;
	for (size_t i=0; i<count; i++)
        {	  
	  atomes[cnt].numeroAtome = cnt+1;
	  strcpy(atomes[cnt].typeAtome, type_vec[cnt].c_str());
	  atomes[cnt].Position[0] = pos_vec[cnt].x;
	  atomes[cnt].Position[1] = pos_vec[cnt].y;
	  atomes[cnt].Position[2] = pos_vec[cnt].z;
	  atomes[cnt].Vitesse[0] = vel_vec[cnt].x;
      	  atomes[cnt].Vitesse[1] = vel_vec[cnt].y;
      	  atomes[cnt].Vitesse[2] = vel_vec[cnt].z;
      	  strcpy(molecules[cnt].typeAtomeFF, type_vec[cnt].c_str());
      	  molecules[cnt].numeroMolecule = molid_vec[cnt];
      	  molecules[cnt].typeMolecule = moltype_vec[cnt];
	  molecules[cnt].Connectivite[0] = cmol_vec[cnt][0];
      	  molecules[cnt].Connectivite[1] = cmol_vec[cnt][1];
      	  molecules[cnt].Connectivite[2] = cmol_vec[cnt][2];
      	  molecules[cnt].Connectivite[3] = cmol_vec[cnt][3];
      	  molecules[cnt].Connectivite[4] = -1;
      	  cnt++;	  
        }

      	// for (size_t i = 0;i<count; i++)
      	//   {
      	//     lout
      	//       << "\t\tatom id "       << std::setw( 6) << atomes[i].numeroAtome 
      	//       << " - typeFF "         << std::setw( 6) << molecules[i].typeAtomeFF 
      	//       << " - molecule id "    << std::setw( 6) << molecules[i].numeroMolecule 
      	//       << " - molecule type "  << std::setw( 6) << molecules[i].typeMolecule 
      	//       << " - connectivity ["  << std::setw( 6) << molecules[i].Connectivite[0] << " ; "
      	//       << std::setw( 6) << molecules[i].Connectivite[1] << " ; "
      	//       << std::setw( 6) << molecules[i].Connectivite[2] << " ; "
      	//       << std::setw( 6) << molecules[i].Connectivite[3] << " ; "
      	//       << std::setw( 6) << molecules[i].Connectivite[4] << "] "
      	//       << std::endl;
      	//   }
      	// lout << std::endl;

      	// for (size_t i = 0;i<count; i++)
      	//   {
      	//     lout
      	//       << "\t\tatom id "       << std::setw( 6) << atomes[i].numeroAtome 
      	//       << " - x "         << std::setw( 6) << atomes[i].Position[0] 
      	//       << " - y "    << std::setw( 6) << atomes[i].Position[1]
      	//       << " - z "  << std::setw( 6) << atomes[i].Position[2]
      	//       << std::endl;
      	//   }
      	// lout << std::endl;	
	
	// End of pushing info to atomes and molecule
	
      	auto type_name_to_index = [&species]( const char* typeAtome ) -> int 
      				  {
      				    int nSpecies = species.size();
      				    int at_type = 0;
      				    for(at_type = 0; at_type < nSpecies ; ++at_type)
      				      {
      					if( species[at_type].name() == typeAtome ) return at_type;
      				      }
      				    lerr << "Atom of type " << typeAtome << " is unknown, the "<<nSpecies<<" available types are :" << std::endl;
      				    for(at_type = 0; at_type < nSpecies ; ++at_type) lerr << "\t" << species[at_type].name()<<std::endl << std::flush;
      				    std::abort();
      				    return -1;
      				  };
	
      	// Start stuff like in read_stamp_v4 for molecules structure
	
      	//==========================================================================================
      	//==========================================================================================
      	//==========================================================================================
	
      	// ------------------------------------------------
      	// --------- comptage des types de molecule -------
      	std::vector<int> at_mol_place( count , 0 );

      	std::vector<int64_t> indices(count);
      	at_mol_place.assign( count , -1 );
      	for (size_t i = 0;i<count; i++) indices[i]=i;

      	std::stable_sort( indices.begin(), indices.end() ,
      			  [&molecules,&atomes](int64_t a, int64_t b)->bool
      			  {
      			    return molecules[a].numeroMolecule < molecules[b].numeroMolecule
      								 || ( molecules[a].numeroMolecule == molecules[b].numeroMolecule && atomes[a].numeroAtome < atomes[b].numeroAtome );
      			  } );
      
      	std::set< std::set< onika::oarray_t<int64_t,6> > > conn_types;
      	std::map< std::set< onika::oarray_t<int64_t,6> > , uint64_t > molecule_type_map;
      	std::vector< onika::oarray_t<int64_t,4> > curmol_conn_types;
      	std::vector<int64_t> curmol_at_ids;
      	std::vector<int> curmol_at_types;
      	std::vector<int> curmol_at_indices;
	
      	int64_t lastMolecule = -1;
      	int64_t molnum = -1;
      	uint64_t molecule_type_id = 0;
      	uint64_t molecule_instance_id = 0;

      	for (size_t idx=0;idx<=count;idx++)
      	  {
      	    int64_t i = -1;
      	    molnum = -1;
      	    if( idx < count )
      	      {
      		i = indices[idx];
      		molnum = molecules[i].numeroMolecule;
      	      }

      	    if( molnum != lastMolecule )
      	      {
      	    	if( lastMolecule != -1 && ! curmol_conn_types.empty() )
      	    	  {
      	    	    uint64_t moltype = 0;
      	    	    int64_t min_id = std::numeric_limits<int64_t>::max();
      	    	    int64_t max_id = -1;
      	    	    for(auto id:curmol_at_ids)
      	    	      {
      	    	    	min_id = std::min( min_id , id );
      	    	    	max_id = std::max( max_id , id );
      	    	      }
      	    	    size_t natoms = max_id-min_id+1;
      	    	    if( natoms != curmol_at_ids.size() )
      	    	      {
      	    	    	lerr<<"curmol_at_ids = [";
      	    	    	for(const auto& id:curmol_at_ids) lerr<<" "<<id;
      	    	    	lerr<<" ] , min_id="<<min_id<<" , max_id="<<max_id<<" , natoms="<<natoms <<std::endl;
      	    	    	fatal_error() << "atom ids not contiguous" << std::endl;
      	    	      }
            
      	    	    std::set< onika::oarray_t<int64_t,6> > molgraph;
      	    	    for(unsigned int ag=0;ag<natoms;ag++)
      	    	      {
      	    	    	auto ct = curmol_conn_types[ag];
      	    	    	for(auto& x:ct) if(x!=-1) x -= min_id;
      	    	    	onika::oarray_t<int64_t,6> gr = { curmol_at_ids[ag]-min_id , curmol_at_types[ag] , ct[0] , ct[1] , ct[2] , ct[3] };
      	    	    	molgraph.insert( gr );
      	    	      }
      	    	    if( conn_types.find(molgraph) == conn_types.end() )
      	    	      {
      	    	    	//lout << "insert molgraph from molid "<< lastMolecule << " with type id "<<molecule_type_id <<std::endl;
      	    	    	conn_types.insert( molgraph );
      	    	    	moltype = molecule_type_id++;
      	    	    	molecule_type_map[ molgraph ] = moltype;
      	    	      }
      	    	    else
      	    	      {
      	    	    	moltype = molecule_type_map[ molgraph ];
      	    	      }
      	    	    for(auto molatidx:curmol_at_indices)
      	    	      {
      	    	    	molecules[molatidx].typeMolecule = moltype;
      	    	    	molecules[molatidx].numeroMolecule = molecule_instance_id;
      	    	    	at_mol_place[molatidx] = atomes[molatidx].numeroAtome - min_id;
      	    	      }
      	    	    ++ molecule_instance_id;
      	    	  }
      	    	curmol_conn_types.clear();
      	    	curmol_at_ids.clear();
      	    	curmol_at_types.clear();
      	    	curmol_at_indices.clear();
      	    	lastMolecule = molnum;
      	      }
      	    if( i != -1 )
      	      {
      	    	//atom connectivity : note, ids are INT in stamp, not UINT.
      		auto cmol = onika::oarray_t<int64_t,4> {
      							molecules[i].Connectivite[0] ,
      							molecules[i].Connectivite[1] ,
      							molecules[i].Connectivite[2] ,
      							molecules[i].Connectivite[3] };		
      		curmol_conn_types.push_back( cmol );
      	    	curmol_at_ids.push_back( atomes[i].numeroAtome );
      	    	curmol_at_types.push_back( type_name_to_index(atomes[i].typeAtome) );
      	    	curmol_at_indices.push_back( i );
      	      }
      	  }

      	//lout << "different molecule graphs = "<<conn_types.size()<<std::endl;
      	molspecies.resize( conn_types.size() );
      	//int ngraphs = 0;
      	for(const auto& molgraph:conn_types)
      	  {
      	    int ngraphs = molecule_type_map[ molgraph ];
      	    molspecies[ngraphs].m_nb_atoms = molgraph.size();
      	    int a = 0;
      	    for(const auto& cmol:molgraph)
      	      {
      		molspecies[ngraphs].m_atom_type[a] = cmol[1];
      		for(int k=0;k<4;k++) molspecies[ngraphs].m_atom_connectivity[a][k] = cmol[k+2];
      		++a;
      	      }
        
      	    if( molspecies[ngraphs].name().empty() )
      	      {
      		std::ostringstream oss;
      		oss << "mol_"<<ngraphs;
      		molspecies[ngraphs].set_name( oss.str() );
      	      }
      	    molspecies[ngraphs].update_connectivity();
      	    molspecies[ngraphs].print( ldbg , species );
      	    //++ ngraphs;
      	  }

      	//==========================================================================================
      	//==========================================================================================
      	//==========================================================================================

      	static const double vel_conv = UnityConverterHelper::convert(1.0, "m/s");
      	Vec3d SumVelocity = {0.,0.,0.};

      	for (size_t i = 0;i<count; i++)
      	  {
      	    //atom absolute ID
      	    int64_t at_id  = atomes[i].numeroAtome;
      	    if( at_id < 0 )
      	      {
      		fatal_error() << "read_xyz_file_with_xform_molecules : Warning: bad id ("<<at_id<<") found" << std::endl;
      	      }

      	    //atom type ID
      	    int at_type_i = type_name_to_index(atomes[i].typeAtome);
      	    if( at_type_i == -1 )
      	      {
      	    	lerr << "\t\tatom id="<< std::setw( 6) << atomes[i].numeroAtome<<" index="<<i
      	    	     << " type="<< std::setw( 6) << atomes[i].typeAtome << "("<< at_type_i <<")"
      	    	     << " pos=("<< std::setw(14) << std::setprecision(6) << atomes[i].Position[0] << "," 
      	    	     << std::setw(14) << std::setprecision(6) << atomes[i].Position[1] << "," 
      	    	     << std::setw(14) << std::setprecision(6) << atomes[i].Position[2] <<")"
      	    	     << " vel=("<< std::setw(14) << std::setprecision(6) << atomes[i].Vitesse[0] << ","
      	    	     << std::setw(14) << std::setprecision(6) << atomes[i].Vitesse[1] << ","
      	    	     << std::setw(14) << std::setprecision(6) << atomes[i].Vitesse[2] << ")"
      	    	     << std::endl ;
      	    	std::abort();
      	      }
      	    uint8_t at_type = at_type_i;

      	    //atom position
      	    Vec3d r = { atomes[i].Position[0] , atomes[i].Position[1] , atomes[i].Position[2] };
      	    IJK loc = domain_periodic_location( domain, r );// - grid.offset();

      	    if( ! grid.contains(loc) )
      	      {
      	    	lerr<<"Domain = "<<domain<<std::endl;
      	    	lerr<<"Domain size = "<<domain.bounds_size()<<std::endl;
      	    	lerr<<"particle #"<<at_id<<", ro="<<r<<", r="<<r<<"<< in cell "<<loc<<" not in grid : offset="<<grid.offset()<<std::endl<<std::flush;
      	    	std::abort();
      	      }
	    
      	    if( ! is_inside(grid.cell_bounds(loc),r) )
      	      {
      	    	lerr<<"particle #"<<at_id<<", ro="<<r<<", r="<<r<<"<< in cell "<<loc<<" not inside cell bounds ="<<grid.cell_bounds(loc)<<std::endl<<std::flush;
      	    	std::abort();
      	      }

      	    //atom velocity
      	    Vec3d v = {atomes[i].Vitesse[0] * vel_conv, atomes[i].Vitesse[1] * vel_conv, atomes[i].Vitesse[2] * vel_conv};
      	    SumVelocity += v;

      	    uint64_t mol_id = 0;
      	    uint64_t molnum = 0;
      	    unsigned int moltype = 0;
      	    unsigned int molplace = 0;
      	    std::array<uint64_t,4> cmol = { 0, 0, 0, 0 };
	    
      	    //molecular ID : note, ids are INT in stamp, not UINT.
      	    molnum = static_cast<uint64_t>( molecules[i].numeroMolecule);
      	    moltype = static_cast<unsigned int>( molecules[i].typeMolecule);
      	    assert( at_mol_place[i] >= 0 );
      	    molplace = at_mol_place[i];
      	    mol_id = make_molecule_id( molnum, molplace, moltype );

      	    //atom connectivity : note, ids are INT in stamp, not UINT.
      	    cmol = std::array<uint64_t,4>{static_cast<uint64_t>(molecules[i].Connectivite[0]),
      	    				  static_cast<uint64_t>(molecules[i].Connectivite[1]),
      	    				  static_cast<uint64_t>(molecules[i].Connectivite[2]),
      	    				  static_cast<uint64_t>(molecules[i].Connectivite[3])};
	    
      	    Quaternion orient = { 0., 0., 0., 0. };
      	    Vec3d angmom = { 0., 0., 0. };
      	    MoleculeTuple t = MoleculeTupleIO(r.x, r.y, r.z, v.x, v.y, v.z, at_id, at_type, orient, angmom, mol_id, cmol);

      	    grid.cell(loc).push_back( t );
	    
      	  }
      	SumVelocity /= count;
      	lout << "Average velocity  = "<<SumVelocity<<std::endl;
	
      	//==========================================================================================
      	//==========================================================================================
      	//==========================================================================================

      }

      grid.set_offset( IJK{0,0,0} );
      grid.set_origin( domain.bounds().bmin );
      grid.set_cell_size( domain.cell_size() );
      grid.set_dimension( domain.grid_dimension() );
     
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
    OperatorNodeFactory::instance()->register_factory("read_xyz_file_with_xform_molecules", make_grid_variant_operator< ReadXYZwXFormMoleculesNode >);
  }

}
