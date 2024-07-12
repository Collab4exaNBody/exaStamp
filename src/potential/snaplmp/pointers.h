#pragma once

#include <memory>
#include <cstdint>
#include <map>
#include <iostream>
#include <algorithm>
#include <onika/cuda/cuda_math.h>
#include <onika/memory/allocator.h>

#define FLERR __FILE__,__LINE__
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

namespace LAMMPS_NS
{
  using bigint = long long;

  struct Size3 { size_t i=0,j=0,k=0; };

  class ArrayDescriptor
  {
  public:
    virtual ~ArrayDescriptor() {}
    virtual size_t count() const =0;
    virtual size_t bytes() const =0;
    virtual Size3 dim() const =0;
    virtual void print() const =0;
  };

  template<class T>  
  class ArrayDescriptorImpl : public ArrayDescriptor
  {
    size_t m_idim = 0;
    size_t m_jdim = 0;
    size_t m_kdim = 0;
    T * m_iptr = nullptr; // pointer to flat array valid for ptr[idx] with idx>=0 && idx<(dim[0]*dim[1]*dim[2])
    T ** m_jptr = nullptr; // arrays of kdim*jdim pointers to i rows
    T *** m_kptr = nullptr; // array of kdim pointers to jdim pointers
    char m_name[32] = { '\0' };

  public:
    inline size_t count() const override final { return m_idim * m_jdim * m_kdim; }
    inline size_t bytes() const override final { return count() * sizeof(T); }
    inline Size3 dim() const override final { return {  m_idim , m_jdim , m_kdim }; }
        
    inline T * ptr(size_t i=0, size_t j=0, size_t k=0)
    {
      assert( i < m_idim && j < m_jdim && k < m_kdim );
      return m_iptr + ( ( k * m_jdim + j ) * m_idim + i );
    }
    inline T ** jptr(size_t j=0, size_t k=0)
    {
      assert( j < m_jdim && k < m_kdim );
      return m_jptr + ( k * m_jdim + j );
    }
    inline T *** kptr(size_t k=0)
    {
      assert( k < m_kdim );
      return m_kptr + k;
    }

    inline ArrayDescriptorImpl( const char* s , size_t idim=1, size_t jdim=1, size_t kdim=1 )
      : m_idim(idim)
      , m_jdim(jdim)
      , m_kdim(kdim)
    {
      std::strncpy(m_name,s,32); m_name[31]='\0';
      const size_t elements = m_idim * m_jdim * m_kdim;
      const size_t data_bytes = ( elements * sizeof(T) + 63ull ) & ( ~ 63ull );
      const size_t pointer_bytes = ( m_kdim * sizeof(T**) ) + ( m_kdim * m_jdim * sizeof(T*) );
      posix_memalign( reinterpret_cast<void**>(&m_iptr), 64, data_bytes + pointer_bytes );
      m_jptr = reinterpret_cast<T**>( reinterpret_cast<uint8_t*>(m_iptr) + data_bytes );
      m_kptr = reinterpret_cast<T***>( m_jptr + ( m_kdim * m_jdim ) );

      for(size_t j=0;j<(m_kdim*m_jdim);j++)
      {
        m_jptr[j] = m_iptr + j*m_idim;
      }
      for(size_t k=0;k<m_kdim;k++)
      {
        m_kptr[k] = m_jptr + k*m_jdim;
      }
      for(size_t k=0;k<m_kdim;k++) for(size_t j=0;j<m_jdim;j++) for(size_t i=0;i<m_idim;i++)
      {
        assert( ptr(i,j,k) == & m_kptr[k][j][i] );
      }
    }

    inline void print() const override final
    {
      const size_t elements = m_idim * m_jdim * m_kdim;
      const size_t data_bytes = ( elements * sizeof(T) + 63ull ) & ( ~ 63ull );
      const size_t pointer_bytes = ( m_kdim * sizeof(T**) ) + ( m_kdim * m_jdim * sizeof(T*) );
      std::cout << "alloc "<<m_name<<" , size=";
      if(m_kdim>1) std::cout<<m_kdim<<"x";
      if(m_jdim>1) std::cout<<m_jdim<<"x";
      std::cout <<m_idim<<"="<<count()<<" , data_bytes="<<data_bytes<<" , pointer_bytes="<<pointer_bytes <<std::endl;
    }
    
    inline ~ArrayDescriptorImpl()
    {
      if( m_iptr != nullptr ) free( m_iptr );
      m_iptr = nullptr;
    }
  };

  struct Memory
  {
    std::map< void* , std::shared_ptr<ArrayDescriptor> > m_allocs;

    inline void destroy(void* p)
    {
      m_allocs.erase( p );
    }
    
    template<class T>
    inline void create(T* &p, size_t idim, const char* name )
    {
      m_allocs.erase( p );
      std::shared_ptr< ArrayDescriptorImpl<T> > array = std::make_shared< ArrayDescriptorImpl<T> >( name, idim );
      p = array->ptr();
      m_allocs[p] = array;
    }

    template<class T>
    inline void create(T** &p, size_t jdim, size_t idim, const char* name )
    {
      m_allocs.erase( p );
      std::shared_ptr< ArrayDescriptorImpl<T> > array = std::make_shared< ArrayDescriptorImpl<T> >( name, idim, jdim );
      p = array->jptr();
      m_allocs[p] = array;
    }

    template<class T>
    inline void create(T*** &p, size_t kdim, size_t jdim, size_t idim, const char* name )
    {
      m_allocs.erase( p );
      std::shared_ptr< ArrayDescriptorImpl<T> > array = std::make_shared< ArrayDescriptorImpl<T> >( name, idim, jdim, kdim );
      p = array->kptr();
      m_allocs[p] = array;
    }
    
    inline void print()
    {
      std::vector< std::shared_ptr<ArrayDescriptor> > desc;
      for(const auto& ap : m_allocs) { desc.push_back( ap.second ); }
      std::sort( desc.begin() , desc.end() , []( const std::shared_ptr<ArrayDescriptor>& a, const std::shared_ptr<ArrayDescriptor>& b ) -> bool { return a->bytes() > b->bytes(); } );
      for(const auto& array : desc) array->print();
    }
  };


}

