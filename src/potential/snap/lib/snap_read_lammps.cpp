#include <exaStamp/potential/snap/snap_config.h>
#include <exaStamp/potential/snap/snap_read_lammps.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <exanb/core/physics_constants.h>
#include <exanb/core/log.h>

namespace SnapExt
{


void snap_read_lammps(const std::string& paramFileName, const std::string& coefFileName, SnapConfig& config , bool conv_units )
{ 
//  std::cout << "read parameters in "<<paramFileName << std::endl;

  config = SnapConfig{}; // set default values

  std::map< std::string , double > values;
  std::map< std::string , double >::iterator it = values.end();
  
  std::string line;
  std::ifstream params( paramFileName );
  if( !params.good() ) { exanb::fatal_error()<<"can't open file '"<<paramFileName<<"' for reading"<<std::endl; }
  while( params.good() )
  {
    std::getline(params,line);
    if( !line.empty() && line.find('#')!=0 )
    {
        std::string key;
        double value=0.;
        std::istringstream(line) >> key >> value;
        values[key] = value;
    }
  }
  params.close();
  
# define SNAP_SET_KEY_VALUE(name) if( (it=values.find(#name)) != values.end() ) config.set_##name( it->second );
  SNAP_SET_KEY_VALUE(rfac0)
  SNAP_SET_KEY_VALUE(rmin0)
  SNAP_SET_KEY_VALUE(bzeroflag)
  SNAP_SET_KEY_VALUE(rcutfac)
  SNAP_SET_KEY_VALUE(twojmax)
  SNAP_SET_KEY_VALUE(nelements)
  SNAP_SET_KEY_VALUE(switchflag)
  SNAP_SET_KEY_VALUE(chemflag)
  SNAP_SET_KEY_VALUE(bnormflag)
  SNAP_SET_KEY_VALUE(wselfallflag)
  SNAP_SET_KEY_VALUE(switchinnerflag)
  SNAP_SET_KEY_VALUE(quadraticflag)
  
# undef SNAP_SET_KEY_VALUE

//  std::cout << "read coefs in "<<coefFileName << std::endl;

  std::ifstream coefs( coefFileName );
  if( ! coefs.good() )
  {
    exanb::fatal_error() << "cannot read file '"<<coefFileName<<"' , aborting"<<std::endl;
  }
  int n_skip_lines = 0;
  std::getline(coefs,line);
  while( !coefs.eof() && (line.find('#')==0 || line.empty()) && n_skip_lines < 100 )
  {
//    std::cout << "line "<<n_skip_lines<<" : eof="<< std::boolalpha << coefs.eof() << ", find('#')="<<line.find('#')<< ", skip line '"<<line<<"'"<< std::endl;
    std::getline(coefs,line);
    ++ n_skip_lines;
  } 
  if( n_skip_lines >= 100 )
  {
    exanb::fatal_error() << "too many lines skipped in file '"<<coefFileName<<"' , aborting"<<std::endl;
  }
  
  size_t n_materials=0;
  size_t coefs_per_material=0;
  
  std::istringstream(line) >> n_materials >> coefs_per_material;

  static const double conv_energy_inv =  1e-4 * exanb::legacy_constant::elementaryCharge / exanb::legacy_constant::atomicMass;

  for(size_t m=0;m<n_materials;m++)
    {
//      std::cout << "Material #" << m << std::endl;
      std::getline(coefs,line);
//      std::cout << line << std::endl;
      if( line.find('#') != 0 )
	{
	  std::string name;
	  double radelem=0., weight=1.;
	  std::istringstream(line) >> name >> radelem >> weight;
//	  std::cout << "name radelem weight :" << name << " " << radelem << " " << weight << std::endl;
	  SnapMaterial mat;
	  mat.set_name(name);
	  mat.set_radelem(radelem);
	  mat.set_weight(weight);
	  mat.resize_coefficients( coefs_per_material );
	  for(size_t c=0;c<coefs_per_material;c++)
	    {
	      double coef = 0.;
	      double coef_converted = 0.;
	      std::getline(coefs,line);
//	      std::cout << "\tcoef "<<c<<" : line='"<<line<<"'"<<std::endl;
	      std::istringstream(line) >> coef;	  
	      coef_converted = coef;
	      if( conv_units ) coef_converted *= conv_energy_inv;
	      mat.set_coefficient(c,coef_converted);
	    }
	  config.materials().push_back(mat);      
	}
    }

  // config.materials().clear();

  // for(size_t m=0;m<n_materials;m++)
  // {
  //   std::getline(coefs,line);
  //   if( line.find('#') != 0 )
  //   {
  //     std::string name;
  //     double radelem=0., weight=1.;
  //     std::istringstream(line) >> name >> radelem >> weight;
  //     std::cout << "name radelem weight :" << name << " " << radelem << " " << weight << std::endl;
  //     SnapMaterial mat;
  //     mat.set_name(name);
  //     mat.set_radelem(radelem);
  //     mat.set_weight(weight);
  //     mat.resize_coefficients( coefs_per_material );
  //     for(size_t c=0;c<coefs_per_material;c++)
  //     {
  //       double coef = 0.;
  //       double coef_converted = 0.;	
  // 	coefs >> coef;
  //       coef_converted = coef;
  // 	      if( conv_units ) coef_converted *= conv_energy_inv;
  //       mat.set_coefficient(c,coef_converted);
  //     }
  //     config.materials().push_back(mat);
  //   }
  // }
  
}


}


