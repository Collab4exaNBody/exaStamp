#pragma once

#include <onika/log.h>

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>

#include <mpi.h>

namespace exaStamp
{

  template<class GridT>
  static inline void snap_check_bispectrum(MPI_Comm comm, const GridT& grid, const std::string& file_name, long ncoeff, const double* bispectrum, double max_l2_error = 1.e-12 )
  {
    std::map< long , std::vector<double> > id_bispectrum_map;
    bool write_check_file = false;
    if( std::ifstream(file_name).good() )
    {
      YAML::Node database = yaml_load_file_abort_on_except(file_name);
      for(auto it : database)
      {
        long id = it.first.as<long>();
        auto str_values = it.second.as< std::vector<std::string> >();
        std::vector<double> values;
        for(const auto & s:str_values) values.push_back( std::strtod( s.c_str() , NULL ) );
        id_bispectrum_map[id] = values;
      }
      write_check_file = false;
    }
    else
    {
      for(int i=1;i<=9;i++) id_bispectrum_map[i*111].clear();
      write_check_file = true;
    }
    
    std::map< long , long > bispectrum_idx_map;
    {
      long idx = 0;
      for(const auto & p : id_bispectrum_map) bispectrum_idx_map[ p.first ] = idx++;
    }
    std::vector<double> bispectrum_computed_values( id_bispectrum_map.size() * ncoeff , 0.0 );
    
    const size_t n_cells = grid.number_of_cells();
    const auto * cell_particle_offset = grid.cell_particle_offset_data();
    for(size_t ci=0;ci<n_cells;ci++)
    {
      if( ! grid.is_ghost_cell(ci) )
      {
        const auto& cell = grid.cell(ci);
        size_t n_particles = cell.size();
        for(size_t pi=0;pi<n_particles;pi++)
        {
          auto p_id = cell[field::id][pi];
          auto it = id_bispectrum_map.find( p_id );
          if( it != id_bispectrum_map.end() )
          {
            long idx = bispectrum_idx_map[ p_id ];
            for(int i=0;i<ncoeff;i++) bispectrum_computed_values[ idx * ncoeff + i ] = bispectrum[ ncoeff * ( cell_particle_offset[ci] + pi ) + i ];
          }
        }
      }
    }
    int rank = 0;
    MPI_Comm_rank(comm,&rank);
    MPI_Allreduce(MPI_IN_PLACE, bispectrum_computed_values.data(), bispectrum_computed_values.size(), MPI_DOUBLE, MPI_SUM, comm );
    if( rank == 0 )
    {
      if( write_check_file )
      {
        exanb::lout << "Write Bispectrum reference file '"<<file_name<<"'"<<std::endl;
        std::ofstream fout(file_name);
        for(const auto & p : bispectrum_idx_map)
        {
          const long idx = p.second;
          const char* s=" ";
          fout << p.first << ": [";
          for(int i=0;i<ncoeff;i++) { fout<<s << std::hexfloat << bispectrum_computed_values[idx*ncoeff+i] ; s=" , "; }
          fout << " ]" << std::endl;
        }
      }
      else
      {
        exanb::lout << "Check Bispectrum data with reference file '"<<file_name<<"'"<<std::endl;
        double max_l2norm = 0.0;
        for(const auto & p : id_bispectrum_map)
        {
          const long idx = bispectrum_idx_map[p.first];
          double l2norm = 0.0;
          for(int i=0;i<ncoeff;i++)
          {
            double d = bispectrum_computed_values[idx*ncoeff+i] - p.second[i];
            l2norm += d*d;
          }
          l2norm = sqrt(l2norm);
          max_l2norm = std::max( max_l2norm , l2norm );
        }
        if( max_l2norm > max_l2_error )
        {
          exanb::fatal_error() << "max L2 norm = "<<max_l2norm<<" , max error = "<<max_l2_error <<std::endl;
        }
      }
    }
  }


}


