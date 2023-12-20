#include <exaStamp/potential/snaplegacy/snap_legacy_config.h>
#include <exaStamp/potential/snaplegacy/snap_legacy_read_lammps.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

void snap_legacy_read_lammps(const std::string& paramFileName, const std::string& coefFileName, SnapConfig& config )
{ 
//  std::cout << "snap_legacy_read_lammps( "<<paramFileName<<" , "<<coefFileName<<" )" std::endl;

  std::map< std::string , double > values;
  values["rfac0"] = 1.0;
  values["rmin0"] = 0.0;
  values["rcutfac"] = 1.0;
  values["twojmax"] = 6;
  
  std::string line;
  std::ifstream params( paramFileName );
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
  
  config.set_rfac0( values["rfac0"] );
  config.set_rmin0( values["rmin0"] );
  config.set_rcutfac( values["rcutfac"] );
  config.set_twojmax( values["twojmax"] );

  std::ifstream coefs( coefFileName );
  do
  {
    line = "";
    std::getline(coefs,line);
    //std::cout << "skip line "<<line<<std::endl;
  } while( !coefs.eof() && (line.find('#')==0 || line.empty()) );
 
  size_t n_materials=0;
  size_t coefs_per_material=0;
  std::istringstream(line) >> n_materials >> coefs_per_material;
  
  config.materials().clear();
  for(size_t m=0;m<n_materials;m++)
  {
    std::getline(coefs,line);
    if( line.find('#') != 0 )
    {
      std::string name;
      double radelem=0., weight=1.;
      std::istringstream(line) >> name >> radelem >> weight;
      SnapMaterial mat;
      mat.set_name(name);
      mat.set_radelem(radelem);
      mat.set_weight(weight);
      mat.resize_coefficients( coefs_per_material );
      double atomicMass = 1.66053904020e-27;  ///< Dalton atomic mass unit in kg
      double elementaryCharge = 1.6021892e-19;  ///< Elementary charge in Coulomb      
      double conv_energy_inv =  1e-4 * elementaryCharge / atomicMass;      
      for(size_t c=0;c<coefs_per_material;c++)
      {
        double coef = 0.;
        double coef_converted = 0.;	
        coefs >> coef;
	coef_converted = coef * conv_energy_inv;	
        mat.set_coefficient(c,coef_converted);
      }
      config.materials().push_back(mat);
    }
  }
  
}


