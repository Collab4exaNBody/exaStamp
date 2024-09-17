#pragma once

#include <onika/cuda/cuda.h>

namespace exaStamp
{

  template<class SizeT>
  struct SnapXSTemporaryComplexArray
  {
    static inline constexpr int Size = 0;
  
    SizeT m_array_size = {};
    double * __restrict__ m_ptr = nullptr;
    
    ONIKA_HOST_DEVICE_FUNC
    inline void init( SizeT sz )
    {
      m_array_size = sz;
    }
    
    ONIKA_HOST_DEVICE_FUNC
    inline void alloc_array()
    {
      if( m_ptr == nullptr ) m_ptr = new double [ m_array_size * 2 ];
    }
    
    ONIKA_HOST_DEVICE_FUNC
    inline double * __restrict__ r()
    {
      alloc_array();
      return m_ptr;
    }

    ONIKA_HOST_DEVICE_FUNC
    inline double * __restrict__ i()
    {
      alloc_array();
      return m_ptr + m_array_size;
    }
    
    ONIKA_HOST_DEVICE_FUNC
    inline ~SnapXSTemporaryComplexArray()
    {
      if( m_ptr != nullptr ) delete [] m_ptr;
    }
  };

  template<int _SZ>
  struct SnapXSTemporaryComplexArray< onika::IntConst<_SZ> >
  {
    static inline constexpr int Size = _SZ;
    double m_storage[ Size * 2 ];

    ONIKA_HOST_DEVICE_FUNC
    inline void init( int sz )
    {
      assert( Size == sz );
    }

    ONIKA_HOST_DEVICE_FUNC
    inline void alloc_array() {}

    ONIKA_HOST_DEVICE_FUNC
    inline double * __restrict__ r()
    {
      return m_storage;
    }

    ONIKA_HOST_DEVICE_FUNC
    inline double * __restrict__ i()
    {
      return m_storage + Size;
    }
  };

  template<class SnapConfT>
  struct SnapXSForceExtStorage
  {  
    SnapXSTemporaryComplexArray< std::remove_cv_t< std::remove_reference_t< decltype( SnapConfT{}.idxu_max * SnapConfT{}.nelements ) > > > m_UTot_array = {};
    SnapXSTemporaryComplexArray< std::remove_cv_t< std::remove_reference_t< decltype( SnapConfT{}.idxu_max_alt * SnapConfT{}.nelements ) > > > m_Y_array    = {};
    SnapXSTemporaryComplexArray<int> m_U_array    = {}; // dynamically allocated only if needed for non specific cases (spedcialized implementations do not use these arrays)
    SnapXSTemporaryComplexArray<int> m_DU_array   = {};

    ONIKA_HOST_DEVICE_FUNC
    inline void init( const SnapConfT& snaconf )
    {
      m_U_array.init( snaconf.idxu_max );
      m_UTot_array.init( snaconf.idxu_max * snaconf.nelements );
      m_DU_array.init( snaconf.idxu_max * onika::IntConst<3>{} );
      m_Y_array.init( snaconf.idxu_max_alt * snaconf.nelements );
    }
  };

  template<class SnapConfT>
  struct SnapBSExtStorage
  {  
    SnapXSTemporaryComplexArray< std::remove_cv_t< std::remove_reference_t< decltype( SnapConfT{}.idxu_max                         ) > > > m_U_array    = {};
    SnapXSTemporaryComplexArray< std::remove_cv_t< std::remove_reference_t< decltype( SnapConfT{}.idxu_max * SnapConfT{}.nelements ) > > > m_UTot_array = {};
    SnapXSTemporaryComplexArray< std::remove_cv_t< std::remove_reference_t< decltype( SnapConfT{}.idxz_max * SnapConfT{}.ndoubles  ) > > > m_Z_array    = {};

    ONIKA_HOST_DEVICE_FUNC
    inline void init( const SnapConfT& snaconf )
    {
      m_U_array.init( snaconf.idxu_max );
      m_UTot_array.init( snaconf.idxu_max * snaconf.nelements );
      m_Z_array.init( snaconf.idxz_max * snaconf.ndoubles );
    }
  };

}

