#pragma once

#include <exanb/core/log.h>

namespace LAMMPS_NS
{
  struct ErrorLogWrapper
  {
    inline void warning( const char* file, int lineno , const std::string& str )
    {
#     pragma omp critical(dbg_mesg)
      {
        ::exanb::lerr << file << " : at line "<<lineno << " : " << str << std::endl;
      }
    }
    
    template<class... T>
    inline void one( const char* file, int lineno , const std::string& str , const T&... args )
    {
#     pragma omp critical(dbg_mesg)
      {
        ::exanb::lerr << file << " : at line "<<lineno << " : " << str << std::endl;
        ( ... , ( ::exanb::lerr << args << std::endl ) ) ;
      }
    }
  };

}

